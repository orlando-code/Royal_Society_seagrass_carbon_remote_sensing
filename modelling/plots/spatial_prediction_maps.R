# Spatial prediction maps for all final models (XGB, GAM, GPR, LR)
#
# Loads final model RDS from output/<cv_regime>/final_models/, builds a
# prediction grid from rasters (same approach as gpr_predictions.R), predicts
# with each model, and saves prediction (and SE where available) maps to
# output/<cv_regime>/predictions/.
#
# Requires: fit_final_models.R has been run (XGB_final.rds, GAM_final.rds, GPR_final.rds, LR_final.rds)
# Outputs:  output/<cv_regime>/predictions/{xgb,gam,gpr,lm}_prediction_map.png,
#            {gpr,gam}_se_map.png (SE maps are available for models that return SE)
#
# Usage: source from the modelling driver, or run standalone.
#
# Optional globals (.GlobalEnv):
#   prediction_map_use_log_fill   — if TRUE, mean maps use log10 colour scale (SE maps stay linear).
#   prediction_map_log_fill_epsilon — floor for log scale lower limit (default 1e-6).
#   case_study_grid_cache_version — change string to invalidate case-study grid caches after edits.
#   case_study_regions — named list of list(name, lon, lat); default order: 2 right, rest bottom
#                        (West Mediterranean last = bottom-right).

if (!exists("seagrass_init_repo", mode = "function", inherits = TRUE)) {
  init_path <- file.path("modelling", "R", "init_repo.R")
  if (!file.exists(init_path)) {
    ff <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
    if (!length(ff)) stop("Run from repo root or with: Rscript /path/to/this_script.R", call. = FALSE)
    script_path <- normalizePath(sub("^--file=", "", ff[[1]]), winslash = "/", mustWork = FALSE)
    init_path <- normalizePath(file.path(dirname(script_path), "..", "R", "init_repo.R"), winslash = "/", mustWork = FALSE)
  }
  if (!file.exists(init_path)) stop("Missing bootstrap helper: modelling/R/init_repo.R", call. = FALSE)
  sys.source(init_path, envir = .GlobalEnv)
}
project_root <- seagrass_init_repo(
  include_helpers = FALSE,
  require_core_inputs = FALSE,
  check_renv = FALSE
)
set.seed(42)

source("modelling/R/helpers.R")
source("modelling/R/extract_covariates_from_rasters.R")
source("modelling/R/assign_region_from_latlon.R")
source("modelling/pipeline_config.R")
load_packages(c("here", "mgcv", "tidyverse", "ggplot2", "GauPro", "xgboost", "maps", "ncdf4"))

cfg <- get_pipeline_config()
apply_pipeline_defaults(
  cfg,
  c(
    "cv_output_dir", "target_var", "exclude_regions", "model_list", "n_lons", "n_lats", "dpi", "show_titles", "prediction_map_use_log_fill", "prediction_map_log_fill_epsilon"
  ),
  envir = .GlobalEnv
)

cv_out    <- get0("cv_output_dir", envir = .GlobalEnv, ifnotfound = "output")
out_dir   <- file.path(cv_out, "predictions")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
final_dir <- file.path(cv_out, "final_models")

target_var <- "median_carbon_density_100cm"
exclude_regions <- get0("exclude_regions", envir = .GlobalEnv, ifnotfound = character(0))
model_list <- get0("model_list", envir = .GlobalEnv, ifnotfound = c("XGB", "GAM", "GPR", "LR"))
n_lons <- as.integer(get0("n_lons", envir = .GlobalEnv, ifnotfound = 1000L))
n_lats <- as.integer(get0("n_lats", envir = .GlobalEnv, ifnotfound = 1000L))
dpi <- get0("dpi", envir = .GlobalEnv, ifnotfound = 150)
show_titles <- get0("show_titles", envir = .GlobalEnv, ifnotfound = TRUE)
prediction_map_use_log_fill <- isTRUE(get0("prediction_map_use_log_fill", envir = .GlobalEnv, ifnotfound = FALSE))
prediction_map_log_fill_epsilon <- as.numeric(get0("prediction_map_log_fill_epsilon", envir = .GlobalEnv, ifnotfound = 1e-6))

# -----------------------------------------------------------------------------
# Data and grid extent
# -----------------------------------------------------------------------------
cat("\n========================================\n")
cat("SPATIAL PREDICTION MAPS (XGB, GAM, GPR, LR)\n")
cat("========================================\n\n")

dat <- readr::read_rds("data/all_extracted_new.rds")
if (length(exclude_regions) > 0L) {
  if (!"region" %in% names(dat)) dat <- assign_region_from_latlon(dat)
  dat <- dat[is.na(dat$region) | !dat$region %in% exclude_regions, ]
}
dat <- process_rs_covariates(dat)

lon_min <- min(dat$longitude, na.rm = TRUE)
lon_max <- max(dat$longitude, na.rm = TRUE)
lat_min <- min(dat$latitude, na.rm = TRUE)
lat_max <- max(dat$latitude, na.rm = TRUE)
cat("Data extent: lon [", lon_min, ",", lon_max, "] lat [", lat_min, ",", lat_max, "]\n")
cat("Grid resolution:", n_lons, "x", n_lats, "\n\n")

# -----------------------------------------------------------------------------
# Prediction grid from rasters (same as gpr_predictions.R)
# -----------------------------------------------------------------------------
bathy_covariate <- grep("bathy|gebco|static_depth", raster_covariates, ignore.case = TRUE, value = TRUE)
reference_raster <- raster_covariates[1]
if (length(bathy_covariate) > 0L) reference_raster <- bathy_covariate[1]

prediction_grid <- create_prediction_grid_from_rasters(
  lon_range = c(lon_min, lon_max),
  lat_range = c(lat_min, lat_max),
  covariates = raster_covariates,
  n_lon = n_lons,
  n_lat = n_lats,
  reference_raster = reference_raster,
  cache_path = sprintf("output/cache/prediction_grid_cache_lons%s_lats%s.rds", n_lons, n_lats),
  use_cache = TRUE,
  show_progress = TRUE,
  fill_coastal = TRUE
)
cat("Prediction grid:", nrow(prediction_grid), "cells\n")

# Bathymetry filter
if (length(bathy_covariate) > 0L && bathy_covariate[1] %in% names(prediction_grid)) {
  bathy_col <- bathy_covariate[1]
  bathy_vals <- dat[[bathy_col]]
  # clip to 30m

  # Ensure column exists and subset correctly by name
  prediction_grid <- prediction_grid[prediction_grid[[bathy_col]] >= -30 & prediction_grid[[bathy_col]] <= 0, , drop = FALSE]


  # if (!all(is.na(bathy_vals))) {
  #   bathy_95perc <- quantile(bathy_vals, 0.01, na.rm = TRUE)
  #   bathy_max <- 1 * bathy_95perc
  #   prediction_grid <- prediction_grid[
  #     prediction_grid[[bathy_col]] >= bathy_max & prediction_grid[[bathy_col]] <= 0,
  #   ]
  #   cat("Filtered to bathymetry [", bathy_max, ", 0] m:", nrow(prediction_grid), "cells\n")
  # }
}
prediction_grid <- process_rs_covariates(prediction_grid)

# Categorical placeholder for grid (models may expect seagrass_species)
if (!"seagrass_species" %in% names(prediction_grid)) {
  prediction_grid$seagrass_species <- "Unspecified"
}

# -----------------------------------------------------------------------------
# Load final models and predict on grid
# -----------------------------------------------------------------------------
lon_range <- c(lon_min, lon_max)
lat_range <- c(lat_min, lat_max)
world <- ggplot2::map_data("world")
obs_for_map <- dat[, c("longitude", "latitude", target_var), drop = FALSE]
obs_for_map <- obs_for_map[complete.cases(obs_for_map), ]

# Color scale limits for all panels in a composite figure.
compute_map_limits <- function(v, use_log = FALSE, log_eps = 1e-6) {
  v <- v[is.finite(v)]
  if (length(v) == 0L) return(c(0, 1))
  upper <- if (any(v > 0.5, na.rm = TRUE)) 0.25 else max(v, na.rm = TRUE)
  lower <- min(v, na.rm = TRUE)
  if (!is.finite(lower)) lower <- 0
  if (lower > upper) lower <- upper
  if (isTRUE(use_log)) {
    if (upper <= 0) return(c(log_eps, max(log_eps * 10, 1)))
    lower <- max(lower, log_eps)
    if (lower >= upper) lower <- upper / 10
    if (lower <= 0) lower <- log_eps
    c(lower, upper)
  } else {
    c(lower, upper)
  }
}

# Single fill scale for rasters (observation points must not use fill = value or ggplot adds a second colourbar).
scale_fill_prediction <- function(fill_limits, fill_name, use_log = FALSE) {
  if (isTRUE(use_log)) {
    ggplot2::scale_fill_viridis_c(
      option = "turbo",
      name = fill_name,
      na.value = "transparent",
      limits = fill_limits,
      oob = scales::squish,
      trans = "log10"
    )
  } else {
    ggplot2::scale_fill_viridis_c(
      option = "turbo",
      name = fill_name,
      na.value = "transparent",
      limits = fill_limits,
      oob = scales::squish
    )
  }
}



make_spatial_panel <- function(grid, fill_col, lon_lim, lat_lim, world, fill_limits,
                               fill_name = "Predicted carbon density", panel_title = NULL,
                               obs_data = NULL, obs_value_col = target_var,
                               show_obs_points = TRUE, show_legend = TRUE,
                               case_panel = FALSE, use_log_fill = FALSE) {
  p <- ggplot2::ggplot() +
    ggplot2::geom_raster(
      data = grid,
      ggplot2::aes(x = longitude, y = latitude, fill = .data[[fill_col]])
    ) +
    scale_fill_prediction(fill_limits, fill_name, use_log = use_log_fill) +
    ggplot2::geom_polygon(
      data = world,
      ggplot2::aes(x = long, y = lat, group = group),
      fill = "#eeeeee", colour = "#a5a5a5"
    )

  # Do not map observations to fill — that creates a second fill scale / duplicate colourbar.
  if (isTRUE(show_obs_points) && !is.null(obs_data) && nrow(obs_data) > 0L) {
    p <- p + ggplot2::geom_point(
      data = obs_data,
      ggplot2::aes(x = longitude, y = latitude),
      shape = 21, size = 1.15, fill = "white", colour = "black", stroke = 0.45, na.rm = TRUE
    )
  }

  p <- p +
    ggplot2::coord_sf(xlim = lon_lim, ylim = lat_lim) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = NULL, y = NULL, title = panel_title) +
    ggplot2::theme(axis.title = ggplot2::element_blank())
  if (isTRUE(case_panel)) {
    p <- p + ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank()
    )
  }
  if (!show_legend) {
    p <- p + ggplot2::theme(legend.position = "none")
  } else {
    p <- p + ggplot2::theme(legend.position = "bottom", legend.key.width = ggplot2::unit(2, "cm"))
  }
  p
}

km_to_degree_counts <- function(lon_lim, lat_lim, cell_km = 1.0) {
  mid_lat <- mean(lat_lim)
  km_per_deg_lat <- 111.32
  km_per_deg_lon <- 111.32 * cos(mid_lat * pi / 180)
  km_per_deg_lon <- max(km_per_deg_lon, 1e-6)
  lon_km <- abs(diff(lon_lim)) * km_per_deg_lon
  lat_km <- abs(diff(lat_lim)) * km_per_deg_lat
  # Target ~cell_km spacing along each axis (regular seq length.out = n_*).
  n_lon <- as.integer(max(64L, ceiling(lon_km / cell_km) + 2L))
  n_lat <- as.integer(max(64L, ceiling(lat_km / cell_km) + 2L))
  list(n_lon = n_lon, n_lat = n_lat)
}

# Back-compatible helper for standalone mean/se maps.
plot_prediction_map <- function(grid, mean_col, se_col = NULL,
                                mean_title = "Predictions", se_title = "Standard error",
                                mean_name = "Predicted carbon density",
                                se_name = "Standard error",
                                lon_range, lat_range, obs_data = NULL, value_col = target_var,
                                use_log_fill_mean = FALSE) {
  if (nrow(grid) == 0L || !mean_col %in% names(grid)) return(list(p_mean = NULL, p_se = NULL))
  mean_limits <- compute_map_limits(
    grid[[mean_col]],
    use_log = use_log_fill_mean,
    log_eps = prediction_map_log_fill_epsilon
  )
  p_mean <- make_spatial_panel(
    grid = grid, fill_col = mean_col, lon_lim = lon_range, lat_lim = lat_range,
    world = world, fill_limits = mean_limits, fill_name = mean_name,
    panel_title = mean_title, obs_data = obs_data, obs_value_col = value_col,
    show_obs_points = TRUE, show_legend = TRUE,
    case_panel = FALSE, use_log_fill = use_log_fill_mean
  )

  p_se <- NULL
  if (!is.null(se_col) && se_col %in% names(grid)) {
    se_limits <- compute_map_limits(grid[[se_col]], use_log = FALSE)
    p_se <- make_spatial_panel(
      grid = grid, fill_col = se_col, lon_lim = lon_range, lat_lim = lat_range,
      world = world, fill_limits = se_limits, fill_name = se_name,
      panel_title = se_title, obs_data = obs_data, obs_value_col = value_col,
      show_obs_points = TRUE, show_legend = TRUE,
      case_panel = FALSE, use_log_fill = FALSE
    )
  }
  list(p_mean = p_mean, p_se = p_se)
}

write_prediction_grid_netcdf <- function(grid_df, prediction_slices, out_path, model_name, target_var) {
  if (!requireNamespace("ncdf4", quietly = TRUE)) {
    warning("Package 'ncdf4' is not available; skipping NetCDF export for ", model_name)
    return(invisible(NULL))
  }
  req_cols <- c("longitude", "latitude")
  if (!all(req_cols %in% names(grid_df))) {
    warning("Missing required columns for NetCDF export: ", paste(setdiff(req_cols, names(grid_df)), collapse = ", "))
    return(invisible(NULL))
  }
  if (length(prediction_slices) == 0L) {
    warning("No prediction slices provided for NetCDF export: ", model_name)
    return(invisible(NULL))
  }

  lon_vals <- sort(unique(grid_df$longitude))
  lat_vals <- sort(unique(grid_df$latitude))
  if (length(lon_vals) == 0L || length(lat_vals) == 0L) {
    warning("No lon/lat values available for NetCDF export: ", model_name)
    return(invisible(NULL))
  }

  lon_idx <- match(grid_df$longitude, lon_vals)
  lat_idx <- match(grid_df$latitude, lat_vals)
  valid_idx <- is.finite(lon_idx) & is.finite(lat_idx)

  species_vals <- names(prediction_slices)
  if (is.null(species_vals) || any(!nzchar(species_vals))) {
    species_vals <- paste0("species_", seq_along(prediction_slices))
  }
  # Replace all spaces in species_vals with underscores
  species_vals <- gsub(" ", "_", species_vals, fixed = TRUE)
  cat("  prediction_species_index mapping:\n")
  for (i in seq_along(species_vals)) {
    cat("    ", i, " -> ", species_vals[[i]], "\n", sep = "")
  }
  n_species <- length(prediction_slices)

  pred_arr <- array(NA_real_, dim = c(length(lon_vals), length(lat_vals), n_species))
  has_any_se <- any(vapply(prediction_slices, function(s) !is.null(s$se), logical(1)))
  se_arr <- if (has_any_se) array(NA_real_, dim = c(length(lon_vals), length(lat_vals), n_species)) else NULL

  for (i in seq_along(prediction_slices)) {
    slice <- prediction_slices[[i]]
    pred_vals <- as.numeric(slice$mean)
    if (length(pred_vals) != nrow(grid_df)) {
      warning("Prediction length mismatch for species '", species_vals[[i]], "'; skipping slice.")
      next
    }
    pred_mat <- matrix(NA_real_, nrow = length(lon_vals), ncol = length(lat_vals))
    pred_mat[cbind(lon_idx[valid_idx], lat_idx[valid_idx])] <- pred_vals[valid_idx]
    pred_arr[, , i] <- pred_mat

    if (has_any_se) {
      if (!is.null(slice$se) && length(slice$se) == nrow(grid_df)) {
        se_vals <- as.numeric(slice$se)
        se_mat <- matrix(NA_real_, nrow = length(lon_vals), ncol = length(lat_vals))
        se_mat[cbind(lon_idx[valid_idx], lat_idx[valid_idx])] <- se_vals[valid_idx]
        se_arr[, , i] <- se_mat
      } else {
        se_arr[, , i] <- NA_real_
      }
    }
  }

  fill_value <- -9999
  pred_arr[is.na(pred_arr)] <- fill_value
  if (has_any_se) se_arr[is.na(se_arr)] <- fill_value

  lon_dim <- ncdf4::ncdim_def("lon", "degrees_east", vals = lon_vals)
  lat_dim <- ncdf4::ncdim_def("lat", "degrees_north", vals = lat_vals)
  sp_dim <- ncdf4::ncdim_def("prediction_species_index", "", vals = seq_len(n_species), create_dimvar = TRUE)
  var_crs <- ncdf4::ncvar_def(
    name = "crs",
    units = "",
    dim = list(),
    missval = NULL,
    longname = "Coordinate reference system",
    prec = "integer"
  )
  var_pred <- ncdf4::ncvar_def(
    name = "prediction",
    units = if (grepl("carbon", target_var, ignore.case = TRUE)) "g C cm-3" else "",
    dim = list(lon_dim, lat_dim, sp_dim),
    missval = fill_value,
    longname = paste0(model_name, " predicted ", target_var),
    prec = "float"
  )
  var_defs <- list(var_crs, var_pred)
  if (has_any_se) {
    var_se <- ncdf4::ncvar_def(
      name = "prediction_se",
      units = "",
      dim = list(lon_dim, lat_dim, sp_dim),
      missval = fill_value,
      longname = paste0(model_name, " prediction standard error for ", target_var),
      prec = "float"
    )
    var_defs <- c(var_defs, list(var_se))
  }

  # Overwrite the file if it exists
  if (file.exists(out_path)) {
    tryCatch({
      file.remove(out_path)
    }, warning = function(w) {}, error = function(e) {})
  }

  nc <- ncdf4::nc_create(out_path, vars = var_defs, force_v4 = TRUE)
  on.exit(ncdf4::nc_close(nc), add = TRUE)
  ncdf4::ncvar_put(nc, var_crs, 0L)
  ncdf4::ncvar_put(nc, var_pred, pred_arr)
  if (has_any_se) ncdf4::ncvar_put(nc, var_se, se_arr)

  ncdf4::ncatt_put(nc, 0, "title", paste0("Spatial prediction grid for ", model_name))
  ncdf4::ncatt_put(nc, 0, "source_model_file", paste0(model_name, "_final.rds"))
  ncdf4::ncatt_put(nc, 0, "target_variable", target_var)
  ncdf4::ncatt_put(nc, 0, "Conventions", "CF-1.8")
  ncdf4::ncatt_put(nc, 0, "prediction_species_count", n_species)
  ncdf4::ncatt_put(nc, 0, "prediction_species_values", paste(species_vals, collapse = "|"))
  ncdf4::ncatt_put(nc, 0, "history", paste(format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ"), "created by spatial_prediction_maps.R"))
  ncdf4::ncatt_put(nc, 0, "geospatial_lon_min", min(lon_vals))
  ncdf4::ncatt_put(nc, 0, "geospatial_lon_max", max(lon_vals))
  ncdf4::ncatt_put(nc, 0, "geospatial_lat_min", min(lat_vals))
  ncdf4::ncatt_put(nc, 0, "geospatial_lat_max", max(lat_vals))
  ncdf4::ncatt_put(nc, 0, "geospatial_lon_units", "degrees_east")
  ncdf4::ncatt_put(nc, 0, "geospatial_lat_units", "degrees_north")
  ncdf4::ncatt_put(nc, "lon", "standard_name", "longitude")
  ncdf4::ncatt_put(nc, "lon", "long_name", "longitude")
  ncdf4::ncatt_put(nc, "lon", "axis", "X")
  ncdf4::ncatt_put(nc, "lat", "standard_name", "latitude")
  ncdf4::ncatt_put(nc, "lat", "long_name", "latitude")
  ncdf4::ncatt_put(nc, "lat", "axis", "Y")
  ncdf4::ncatt_put(nc, "prediction_species_index", "long_name", "Index of prediction species")
  ncdf4::ncatt_put(nc, "crs", "grid_mapping_name", "latitude_longitude")
  ncdf4::ncatt_put(nc, "crs", "epsg_code", "EPSG:4326")
  ncdf4::ncatt_put(nc, "crs", "crs_wkt", "GEOGCRS[\"WGS 84\",DATUM[\"World Geodetic System 1984\",ELLIPSOID[\"WGS 84\",6378137,298.257223563]],PRIMEM[\"Greenwich\",0],CS[ellipsoidal,2],AXIS[\"longitude\",east],AXIS[\"latitude\",north],ANGLEUNIT[\"degree\",0.0174532925199433],ID[\"EPSG\",4326]]")
  ncdf4::ncatt_put(nc, "crs", "semi_major_axis", 6378137.0)
  ncdf4::ncatt_put(nc, "crs", "inverse_flattening", 298.257223563)
  ncdf4::ncatt_put(nc, "crs", "longitude_of_prime_meridian", 0.0)
  ncdf4::ncatt_put(nc, "crs", "spatial_ref", "EPSG:4326")
  ncdf4::ncatt_put(nc, "prediction", "grid_mapping", "crs")
  ncdf4::ncatt_put(nc, "prediction", "coordinates", "lon lat prediction_species_index")
  if (has_any_se) {
    ncdf4::ncatt_put(nc, "prediction_se", "grid_mapping", "crs")
    ncdf4::ncatt_put(nc, "prediction_se", "coordinates", "lon lat prediction_species_index")
  }

  invisible(out_path)
}


for (model_name in model_list) {
  rds_path <- file.path(final_dir, paste0(model_name, "_final.rds"))
  if (!file.exists(rds_path)) {
    cat("Skip", model_name, ": no final model at", rds_path, "\n")
    next
  }
  obj <- readRDS(rds_path)
  pvars <- obj$predictor_vars %||% obj$best_covariate_set
  if (is.null(pvars) || length(pvars) == 0L) {
    cat("Skip", model_name, ": no predictor_vars in RDS\n")
    next
  }
  # Ensure grid has all predictors (add defaults for categoricals missing from rasters)
  missing_pv <- setdiff(pvars, names(prediction_grid))
  for (v in missing_pv) {
    if (v == "seagrass_species") {
      prediction_grid[[v]] <- "Unspecified"
    } else if (!is.null(obj$encoding) && v %in% names(obj$encoding)) {
      lvls <- obj$encoding[[v]]
      if (length(lvls) > 0L) prediction_grid[[v]] <- lvls[1] else prediction_grid[[v]] <- NA_character_
    } else {
      prediction_grid[[v]] <- NA
    }
  }
  mean_col <- paste0(tolower(model_name), "_mean")
  se_col <- paste0(tolower(model_name), "_se")
  prediction_grid[[mean_col]] <- NA_real_
  prediction_grid[[se_col]] <- NA_real_

  species_levels <- "all"
  if ("seagrass_species" %in% pvars) {
    if (!is.null(obj$encoding) && "seagrass_species" %in% names(obj$encoding) && length(obj$encoding$seagrass_species) > 0L) {
      species_levels <- unique(as.character(obj$encoding$seagrass_species))
    } else if ("seagrass_species" %in% names(dat)) {
      species_levels <- sort(unique(as.character(dat$seagrass_species)))
    } else {
      species_levels <- "Unspecified"
    }
  }
  species_levels <- species_levels[!is.na(species_levels) & nzchar(species_levels)]
  if (length(species_levels) == 0L) species_levels <- "Unspecified"
  cat("  prediction_species:", paste(species_levels, collapse = ", "), "\n")

  prediction_slices <- list()
  default_species <- if ("Unspecified" %in% species_levels) "Unspecified" else species_levels[[1]]
  default_grid_sub <- NULL
  default_res <- NULL

  for (sp in species_levels) {
    grid_work <- prediction_grid
    if ("seagrass_species" %in% pvars) {
      grid_work$seagrass_species <- sp
    }
    pred_ok_sp <- complete.cases(grid_work[, pvars, drop = FALSE])
    if (sum(pred_ok_sp) == 0L) {
      cat("  Skip species", sp, ": no complete-cases on grid for predictors\n")
      next
    }
    grid_sub_sp <- grid_work[pred_ok_sp, , drop = FALSE]
    res_sp <- tryCatch(
      predict_model(obj, grid_sub_sp, se = (model_name %in% c("GPR", "GAM"))),
      error = function(e) { cat("  Predict error (species=", sp, "): ", conditionMessage(e), "\n", sep = ""); NULL }
    )
    if (is.null(res_sp) || is.null(res_sp$mean)) next

    prediction_slices[[sp]] <- list(
      mean = as.numeric(res_sp$mean),
      se = if (!is.null(res_sp$se)) as.numeric(res_sp$se) else NULL,
      pred_ok = pred_ok_sp
    )

    if (identical(sp, default_species)) {
      prediction_grid[[mean_col]] <- NA_real_
      prediction_grid[[mean_col]][pred_ok_sp] <- as.numeric(res_sp$mean)
      if (!is.null(res_sp$se)) {
        prediction_grid[[se_col]] <- NA_real_
        prediction_grid[[se_col]][pred_ok_sp] <- as.numeric(res_sp$se)
      }
      default_grid_sub <- grid_sub_sp
      default_res <- res_sp
    }
  }
  if (length(prediction_slices) == 0L) {
    cat("  No valid predictions for", model_name, "\n")
    next
  }
  if (is.null(default_res)) {
    first_sp <- names(prediction_slices)[1]
    sl <- prediction_slices[[first_sp]]
    grid_work <- prediction_grid
    if ("seagrass_species" %in% pvars) grid_work$seagrass_species <- first_sp
    prediction_grid[[mean_col]] <- NA_real_
    prediction_grid[[mean_col]][sl$pred_ok] <- sl$mean
    if (!is.null(sl$se)) {
      prediction_grid[[se_col]] <- NA_real_
      prediction_grid[[se_col]][sl$pred_ok] <- sl$se
    }
    default_species <- first_sp
  }

  map_data <- prediction_grid[!is.na(prediction_grid[[mean_col]]) & complete.cases(prediction_grid[, c("longitude", "latitude")]), ]
  if (nrow(map_data) == 0L) {
    cat("  No valid predictions for", model_name, "\n")
    next
  }

  plots <- plot_prediction_map(
    map_data, mean_col, se_col,
    mean_title = if (show_titles) paste(model_name, "predictions (species:", default_species, ")") else NULL,
    se_title = if (show_titles) paste(model_name, "prediction uncertainty (SE, species:", default_species, ")") else NULL,
    lon_range = lon_range, lat_range = lat_range, obs_data = obs_for_map, value_col = target_var,
    use_log_fill_mean = prediction_map_use_log_fill
  )
  if (!is.null(plots$p_mean)) {
    f_mean <- file.path(out_dir, paste0(tolower(model_name), "_prediction_map.png"))
    ggsave(f_mean, plots$p_mean, width = 10, height = 10, dpi = dpi)
    cat("  Saved", basename(f_mean), "\n")
  }
  if (!is.null(plots$p_se)) {
    f_se <- file.path(out_dir, paste0(tolower(model_name), "_se_map.png"))
    ggsave(f_se, plots$p_se, width = 10, height = 10, dpi = dpi)
    cat("  Saved", basename(f_se), "\n")
  }

  f_nc <- file.path(out_dir, paste0(tolower(model_name), "_prediction_grid.nc"))
  netcdf_grid <- prediction_grid[complete.cases(prediction_grid[, c("longitude", "latitude")]), c("longitude", "latitude"), drop = FALSE]
  netcdf_slices <- lapply(prediction_slices, function(sl) {
    mean_full <- rep(NA_real_, nrow(netcdf_grid))
    se_full <- rep(NA_real_, nrow(netcdf_grid))
    idx <- which(sl$pred_ok)
    key <- paste(netcdf_grid$longitude, netcdf_grid$latitude)
    key_sp <- paste(prediction_grid$longitude[idx], prediction_grid$latitude[idx])
    map_idx <- match(key_sp, key)
    ok <- is.finite(map_idx)
    mean_full[map_idx[ok]] <- sl$mean[ok]
    if (!is.null(sl$se)) se_full[map_idx[ok]] <- sl$se[ok]
    list(mean = mean_full, se = if (all(is.na(se_full))) NULL else se_full)
  })
  write_prediction_grid_netcdf(
    grid_df = netcdf_grid,
    prediction_slices = netcdf_slices,
    out_path = f_nc,
    model_name = model_name,
    target_var = target_var
  )
  if (file.exists(f_nc)) cat("  Saved", basename(f_nc), "\n")
}

cat("\nSpatial prediction maps complete. Outputs in", out_dir, "\n")
