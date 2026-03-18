# Spatial prediction maps for all final models (XGB, GAM, GPR)
#
# Loads final model RDS from output/final_models/, builds a prediction grid from
# rasters (same approach as gpr_predictions.R), predicts with each model, and
# saves prediction (and SE where available) maps to output/predictions/.
#
# Requires: fit_final_models.R has been run (XGB_final.rds, GAM_final.rds, GPR_final.rds)
# Outputs:  output/predictions/{xgb,gam,gpr}_prediction_map.png, gpr_se_map.png
#
# Usage: source from run_paper.R (step 7), or run standalone.

setwd(here::here())
set.seed(42)

source("modelling/R/helpers.R")
source("modelling/R/extract_covariates_from_rasters.R")
source("modelling/R/assign_region_from_latlon.R")
load_packages(c("here", "mgcv", "tidyverse", "ggplot2", "GauPro", "xgboost", "maps"))

cv_out    <- get0("cv_output_dir", envir = .GlobalEnv, ifnotfound = "output")
out_dir   <- file.path(cv_out, "predictions")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
final_dir <- file.path(cv_out, "final_models")

target_var <- "median_carbon_density_100cm"
exclude_regions <- get0("exclude_regions", envir = .GlobalEnv, ifnotfound = character(0))
model_list <- get0("model_list", envir = .GlobalEnv, ifnotfound = c("XGB", "GAM", "GPR"))
n_lons <- as.integer(get0("n_lons", envir = .GlobalEnv, ifnotfound = 500L))
n_lats <- as.integer(get0("n_lats", envir = .GlobalEnv, ifnotfound = 500L))
dpi <- get0("dpi", envir = .GlobalEnv, ifnotfound = 150)
show_titles <- get0("show_titles", envir = .GlobalEnv, ifnotfound = TRUE)

# -----------------------------------------------------------------------------
# Data and grid extent
# -----------------------------------------------------------------------------
cat("\n========================================\n")
cat("SPATIAL PREDICTION MAPS (XGB, GAM, GPR)\n")
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
  if (!all(is.na(bathy_vals))) {
    bathy_95perc <- quantile(bathy_vals, 0.01, na.rm = TRUE)
    bathy_max <- 2 * bathy_95perc
    prediction_grid <- prediction_grid[
      prediction_grid[[bathy_col]] >= bathy_max & prediction_grid[[bathy_col]] <= 0,
    ]
    cat("Filtered to bathymetry [", bathy_max, ", 0] m:", nrow(prediction_grid), "cells\n")
  }
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

# Helper: plot mean (and optionally SE) map with given column names
plot_prediction_map <- function(grid, mean_col, se_col = NULL,
                                mean_title = "Predictions", se_title = "Standard error",
                                mean_name = "Predicted carbon density",
                                se_name = "Standard error",
                                lon_range, lat_range, obs_data = NULL, value_col = target_var) {
  if (nrow(grid) == 0L || !mean_col %in% names(grid)) return(list(p_mean = NULL, p_se = NULL))
  p_mean <- ggplot2::ggplot() +
    ggplot2::geom_raster(data = grid, ggplot2::aes(x = longitude, y = latitude, fill = .data[[mean_col]])) +
    ggplot2::scale_fill_viridis_c(option = "turbo", name = mean_name, na.value = "transparent") +
    ggplot2::geom_polygon(
      data = world,
      ggplot2::aes(x = long, y = lat, group = group),
      fill = "#eeeeee", colour = "#a5a5a5"
    )
  if (!is.null(obs_data) && nrow(obs_data) > 0L && value_col %in% names(obs_data)) {
    p_mean <- p_mean +
      ggplot2::geom_point(
        data = obs_data,
        ggplot2::aes(x = longitude, y = latitude, fill = .data[[value_col]]),
        shape = 21, size = 1.3, colour = "white", stroke = 0.5, show.legend = FALSE
      )
  }
  p_mean <- p_mean +
    ggplot2::coord_sf(xlim = lon_range, ylim = lat_range) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom", legend.key.width = ggplot2::unit(2, "cm")) +
    ggplot2::labs(x = "Longitude", y = "Latitude", title = mean_title)

  p_se <- NULL
  if (!is.null(se_col) && se_col %in% names(grid)) {
    p_se <- ggplot2::ggplot() +
      ggplot2::geom_raster(data = grid, ggplot2::aes(x = longitude, y = latitude, fill = .data[[se_col]])) +
      ggplot2::scale_fill_viridis_c(option = "turbo", name = se_name, na.value = "transparent") +
      ggplot2::geom_polygon(
        data = world,
        ggplot2::aes(x = long, y = lat, group = group),
        fill = "#eeeeee", colour = "#a5a5a5"
      )
    if (!is.null(obs_data) && nrow(obs_data) > 0L) {
      p_se <- p_se +
        ggplot2::geom_point(
          data = obs_data, ggplot2::aes(x = longitude, y = latitude),
          shape = 21, size = 1.5, fill = "white", colour = "black", stroke = 1.2, show.legend = FALSE
        )
    }
    p_se <- p_se +
      ggplot2::coord_sf(xlim = lon_range, ylim = lat_range) +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "bottom", legend.key.width = ggplot2::unit(2, "cm")) +
      ggplot2::labs(x = "Longitude", y = "Latitude", title = se_title)
  }
  list(p_mean = p_mean, p_se = p_se)
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
  # Ensure grid has all predictors (add default for categoricals missing from rasters)
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
  pred_ok <- complete.cases(prediction_grid[, pvars, drop = FALSE])
  if (sum(pred_ok) == 0L) {
    cat("Skip", model_name, ": no complete-cases on grid for predictors\n")
    next
  }
  grid_sub <- prediction_grid[pred_ok, , drop = FALSE]
  res <- tryCatch(
    predict_model(obj, grid_sub, se = (model_name == "GPR")),
    error = function(e) { cat("  Predict error:", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(res)) next

  mean_col <- paste0(tolower(model_name), "_mean")
  se_col <- if (model_name == "GPR" && !is.null(res$se)) paste0(tolower(model_name), "_se") else NULL
  prediction_grid[[mean_col]] <- NA_real_
  prediction_grid[[mean_col]][pred_ok] <- res$mean
  if (!is.null(se_col)) {
    prediction_grid[[se_col]] <- NA_real_
    prediction_grid[[se_col]][pred_ok] <- res$se
  }

  map_data <- prediction_grid[!is.na(prediction_grid[[mean_col]]) & complete.cases(prediction_grid[, c("longitude", "latitude")]), ]
  if (nrow(map_data) == 0L) {
    cat("  No valid predictions for", model_name, "\n")
    next
  }
  plots <- plot_prediction_map(
    map_data, mean_col, se_col,
    mean_title = if (show_titles) paste(model_name, "predictions") else NULL,
    se_title = if (show_titles) paste(model_name, "prediction uncertainty (SE)") else NULL,
    lon_range = lon_range, lat_range = lat_range, obs_data = obs_for_map, value_col = target_var
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
}

cat("\nSpatial prediction maps complete. Outputs in", out_dir, "\n")
