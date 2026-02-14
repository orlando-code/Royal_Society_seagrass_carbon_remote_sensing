# Gaussian Process Regression (GPR) Predictions and Spatial Mapping
#
# This script fits a GPR model to seagrass carbon density data and generates
# spatial predictions with uncertainty estimates. Optionally includes species
# and region as categorical variables, and can include spatial smoothing by
# adding longitude/latitude as predictors.

# TODO: better cmaps for spatial predictions
# zorder landmass to top
# decrease length of colourbar on joint plot


## ================================ SETUP ================================
setwd(here::here())
set.seed(42)

source("modelling/R/helpers.R")
# Ensure extractor helpers are defined in the global environment so raster_covariates is visible

source("modelling/R/extract_covariates_from_rasters.R")
source("modelling/R/gpr_funs.R")
source("modelling/R/assign_region_from_latlon.R")
source("modelling/R/plot_config.R")
source("modelling/plots/plot_gpr_spatial_maps.R")
load_packages(c("here", "mgcv", "tidyverse", "ggplot2", "GauPro", "viridisLite", "maps", "patchwork"))

## ================================ CONFIGURATION ================================
# Target variable
target_var <- "median_carbon_density_100cm"

# Prediction grid resolution (recommend 1000 for full paper; 100 for quick test)
n_lons <- 1000
n_lats <- 1000
if (exists("dpi", envir = .GlobalEnv)) dpi <- get("dpi", envir = .GlobalEnv) else dpi <- 150
if (exists("show_titles", envir = .GlobalEnv)) show_titles <- get("show_titles", envir = .GlobalEnv) else show_titles <- TRUE
# Options
# - Species: keep grid predictions species-agnostic (include_species = FALSE).
# - Region: CV showed Region hurts model fit; keep FALSE so predictions do not use region.
# - Spatial smoothing (lat/lon): if TRUE, predictions will revert to the mean of nearby data rather than the global mean.
include_species <- TRUE # Keep grid predictions species-agnostic
include_region <- FALSE # Region hurt CV fit; omit for prediction
include_spatial_smoothing <- TRUE # Lat/lon as predictors → revert to nearby mean (not global)
exclude_regions <- c("Black Sea") # Black Sea contains anomalously-high carbon density values and this leaks to eastern Mediterranean values
# Partial dependence plots
create_partial_dependence_plots <- TRUE # Generate partial dependence plots
pd_variables <- NULL # NULL = plot all numeric predictors, or specify vector of variable names
pd_n_points <- 100 # Number of evaluation points per variable
pd_reference_method <- "conditional" # "median", "mean", "sample", or "conditional" (uses actual data bins)
pd_n_samples <- 1000 # If pd_reference_method = "sample", number of samples to average over
pd_use_conditional <- TRUE # If TRUE, use conditional partial dependence (only actual data combinations)

# Variable importance analysis
compute_permutation_importance <- TRUE # Compute permutation importance (shuffle variable, measure performance drop)
n_permutations <- 5 # Number of permutations per variable (for stability)

# Kernel: use best from tuning if available
gpr_kernel <- "matern52" # Overridden by gpr_best_config.rds if present
use_cache_fitted_model <- TRUE # If TRUE and cache exists, load fitted model; also saves model for PDP script
model_cache_path <- "figures/cv_pipeline_output/gpr_fitted_model.rds"

# PDPs: explicitly over training data variable ranges (optionally restrict to percentile range)
pd_range_quantiles <- c(0.01, 0.99) # e.g. c(0.01, 0.99) to avoid extreme tails
pd_show_training_range <- TRUE # If TRUE, add rug of training values on PDPs

# Output directory
out_dir <- "figures/predictions"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
# Optional: create and save region collage (spatial prediction panels per region)
create_region_collage <- FALSE # Set TRUE to save gpr_region_collage.png; or run via run_paper_figures.R

## ================================ DATA LOADING ================================
cat("\n========================================\n")
cat("LOADING DATA\n")
cat("========================================\n\n")

# Load data with raster-extracted covariates
raster_covariates_dat <- read_rds("data/all_extracted_new.rds")
raster_covariates_dat <- process_rs_covariates(raster_covariates_dat)
# Optionally infer region from lon/lat for diagnostics/stratified summaries
if (!"region" %in% names(raster_covariates_dat) && all(c("longitude", "latitude") %in% names(raster_covariates_dat))) {
  raster_covariates_dat <- assign_region_from_latlon(raster_covariates_dat)
}

# Get coordinate ranges
lon_min <- min(raster_covariates_dat$longitude, na.rm = TRUE)
lon_max <- max(raster_covariates_dat$longitude, na.rm = TRUE)
lat_min <- min(raster_covariates_dat$latitude, na.rm = TRUE)
lat_max <- max(raster_covariates_dat$latitude, na.rm = TRUE)

cat("Data range:\n")
cat("  Longitude:", lon_min, "to", lon_max, "\n")
cat("  Latitude:", lat_min, "to", lat_max, "\n")
cat("  Number of data points:", nrow(raster_covariates_dat), "\n\n")

# Load pruned variables: prefer generic model permutation pruning (model_permutation_pruning.R),
gpr_pruned_vars <- get_pruned_predictors_for_model("GPR", data_cols = names(raster_covariates_dat))
if (!is.null(gpr_pruned_vars) && length(gpr_pruned_vars) >= 2L) {
  cat("Loaded", length(gpr_pruned_vars), "GPR covariates from model permutation pruning (1 km spatial CV).\n")
} else {
  gpr_pruned_vars <- NULL
}
if (is.null(gpr_pruned_vars) || length(gpr_pruned_vars) < 2L) {
  pruned_file_gpr <- "figures/covariate_selection/pruned_variables_to_include_gpr.csv"
  pruned_file <- "figures/covariate_selection/pruned_variables_to_include.csv"
  if (file.exists(pruned_file_gpr)) {
    gpr_pruned_vars <- read_csv(pruned_file_gpr, show_col_types = FALSE)$variable
    cat("Loaded", length(gpr_pruned_vars), "GPR-pruned covariates (gpr_covariate_pruning fallback).\n")
  } else if (file.exists(pruned_file)) {
    gpr_pruned_vars <- read_csv(pruned_file, show_col_types = FALSE)$variable
    cat("Loaded", length(gpr_pruned_vars), "pruned covariates (GAM list fallback).\n")
  } else {
    stop("Pruned variables not found. Run model_permutation_pruning.R, or covariate_pruning.R / gpr_covariate_pruning.R.")
  }
}
gpr_pruned_vars <- gpr_pruned_vars[gpr_pruned_vars %in% names(raster_covariates_dat)]
cat("  ", paste(gpr_pruned_vars, collapse = ", "), "\n\n")


source("modelling/R/plot_config.R")
# Pre-compute histograms with variable-specific binwidth (~50 bins in each variable's range)
n_bins <- 10
hist_bars <- purrr::map_dfr(gpr_pruned_vars, function(var) {
  x <- raster_covariates_dat[[var]]
  x <- x[!is.na(x) & is.finite(x)]
  if (length(x) < 2) {
    return(data.frame(variable = character(), xmin = numeric(), xmax = numeric(), density = numeric(), count = integer()))
  }
  r <- range(x)
  if (diff(r) == 0) r[2] <- r[1] + 1
  breaks <- seq(r[1], r[2], length.out = n_bins + 1)
  h <- hist(x, breaks = breaks, plot = FALSE)
  bin_widths <- diff(h$breaks)
  # Proper density: proportion of values per bin divided by bin width, total sums to 1
  densities <- h$counts / (sum(h$counts) * bin_widths)
  # normalize densities to sum to 1
  densities <- densities / sum(densities)
  data.frame(
    variable = var,
    xmin = h$breaks[-length(h$breaks)],
    xmax = h$breaks[-1],
    density = densities,
    count = h$counts
  )
})

# Faceted plot: one histogram per variable, each with appropriate binwidth (geom_rect for variable bar widths)
p_histograms <- ggplot(hist_bars, aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = density)) +
  geom_rect(fill = "steelblue", color = "black", alpha = 0.7, linewidth = 0.2) +
  facet_wrap(
    ~variable,
    scales = "free_x",
    labeller = as_labeller(function(x) {
      labs <- VAR_LABELS[x]
      labs[is.na(labs)] <- x[is.na(labs)]
      labs
    }))+
  labs(x = "Value", y = "Normalised density") +
  theme_minimal()

print(p_histograms)
ggsave(file.path(out_dir, "covariate_distributions.png"), p_histograms, width = 14, height = 10, dpi = dpi)

cat("\n========================================\n")
cat("COVARIATE DISTRIBUTION ANALYSIS\n")
cat("========================================\n\n")

# Analyze each variable's distribution
sparse_vars <- character()
continuous_vars <- character()

for (var in gpr_pruned_vars) {
  if (var %in% names(raster_covariates_dat)) {
    values <- raster_covariates_dat[[var]]
    values <- values[!is.na(values)]

    if (length(values) > 0) {
      # Calculate coefficient of variation and gaps
      n_unique <- length(unique(values))
      n_total <- length(values)
      unique_ratio <- n_unique / n_total

      # Check for large gaps (sparsity)
      sorted_vals <- sort(unique(values))
      if (length(sorted_vals) > 1) {
        gaps <- diff(sorted_vals)
        max_gap <- max(gaps)
        median_gap <- median(gaps)
        gap_ratio <- max_gap / (median_gap + 1e-10) # Avoid division by zero

        # Classify as sparse if: many unique values relative to total, or large gaps
        if (gap_ratio > 10 || unique_ratio > 0.5) {
          sparse_vars <- c(sparse_vars, var)
          cat("  SPARSE: ", var, "\n")
          cat("    - Unique values: ", n_unique, " / ", n_total,
            " (", round(100 * unique_ratio, 1), "%)\n",
            sep = ""
          )
          cat("    - Max gap ratio: ", round(gap_ratio, 2), "\n")
        } else {
          continuous_vars <- c(continuous_vars, var)
        }
      }
    }
  }
}


# Check for species/region columns if requested
if (include_species && !"seagrass_species" %in% names(raster_covariates_dat)) {
  warning("seagrass_species not found in data. Setting include_species = FALSE")
  include_species <- FALSE
}
if (include_region && !"Region" %in% names(raster_covariates_dat)) {
  warning("Region not found in data. Setting include_region = FALSE")
  include_region <- FALSE
}

## ================================ PREDICTION GRID ================================
cat("========================================\n")
cat("CREATING PREDICTION GRID\n")
cat("========================================\n\n")

# (Re)build a fresh covariate configuration here
source("modelling/R/extract_covariates_from_rasters.R")
bathy_covariate <- grep("bathy|gebco|static_depth", raster_covariates, ignore.case = TRUE, value = TRUE)
# Use first pruned covariate as reference raster (or find one that exists)
reference_raster <- gpr_pruned_vars[1]
if (!reference_raster %in% raster_covariates) {
  reference_raster <- raster_covariates[1]
}
cat("Using", reference_raster, "as reference raster\n")
cat("Using", bathy_covariate, "as bathymetry data\n\n")

cat("Creating prediction grid from rasters...\n")
prediction_grid <- create_prediction_grid_from_rasters(
  lon_range = c(lon_min, lon_max),
  lat_range = c(lat_min, lat_max),
  covariates = raster_covariates,
  n_lon = n_lons,
  n_lat = n_lats,
  reference_raster = reference_raster, # TODO: does the reference raster matter? e.g. its resolution?
  cache_path = sprintf("data/prediction_grid_cache_lons%s_lats%s.rds", n_lons, n_lats),
  use_cache = TRUE, # Use cache if available
  show_progress = TRUE,
  fill_coastal = TRUE
)
cat("Prediction grid: ", nrow(prediction_grid), " cells\n\n")
colnames(prediction_grid)
# Filter to reasonable bathymetry range (if bathy column exists)
if (length(bathy_covariate) > 0) {
  # If bathy_col is a column name, using [[bathy_col]] should work;
  # however, if there's a case/whitespace/naming issue, using $ for exact match helps.
  # Try falling back to $ operator if the first try is NA, also check for data issues.
  bathy_vals <- raster_covariates_dat[[bathy_covariate]]
  if (all(is.na(bathy_vals))) {
    # Try $ if [[bathy_covariate]] failed -- rare, but possible due to special characters
    if (bathy_covariate %in% names(raster_covariates_dat)) {
      bathy_vals <- raster_covariates_dat[[bathy_covariate]]
    } else if (bathy_covariate %in% colnames(raster_covariates_dat)) {
      bathy_vals <- raster_covariates_dat[, bathy_covariate]
    } else {
      stop(
        sprintf(
          "Bathymetry column '%s' not found in raster_covariates_dat. Names are: %s",
          bathy_covariate, paste(names(raster_covariates_dat), collapse = ", ")
        )
      )
    }
  }
  bathy_95perc <- quantile(bathy_vals, 0.01, na.rm = TRUE)
  bathy_max <- 2 * bathy_95perc
  bathy_grid <- prediction_grid[
    prediction_grid[[bathy_covariate]] >= bathy_max &
      prediction_grid[[bathy_covariate]] <= 0,
  ]
  cat("Filtered to bathymetry >= ", bathy_max, "m: ", nrow(bathy_grid), " cells\n\n")
} else {
  bathy_grid <- prediction_grid
  cat("No bathymetry column found; using full prediction grid\n\n")
}
bathy_grid <- process_rs_covariates(bathy_grid)
# Tag prediction grid with inferred region (for plotting/zoomed collages)
if (!"region" %in% names(bathy_grid) && all(c("longitude", "latitude") %in% names(bathy_grid))) {
  bathy_grid <- assign_region_from_latlon(bathy_grid)
}

world <- ggplot2::map_data("world")

# plot bathy_grid bathymetry (gebco)
p_bathy <- ggplot(bathy_grid, aes(x = longitude, y = latitude, fill = gebco_2025_n61.0_s34.0_w_10.0_e35.0)) +
  geom_raster() +
  geom_polygon(
    data = world,
    aes(x = long, y = lat, group = group),
    fill = "#eeeeee",
    colour = "#a5a5a5"
  ) +
  coord_sf(xlim = c(lon_min, lon_max), ylim = c(lat_min, lat_max)) +
  scale_fill_viridis_c() +
  labs(fill = "Bathymetry (m)", x = "Longitude", y = "Latitude") +
  theme_minimal()
print(p_bathy)
ggsave(file.path(out_dir, "bathy_grid_.png"), p_bathy, width = 10, height = 10, dpi = dpi)

## ================================ GPR FIT (see gpr_funs.R) ================================
## ================================ BUILD PREDICTOR SET ================================
cat("========================================\n")
cat("BUILDING PREDICTOR SET\n")
cat("========================================\n\n")

predictor_vars <- gpr_pruned_vars

# Add categorical variables if requested
if (include_species) {
  predictor_vars <- c(predictor_vars, "seagrass_species")
  cat("  Added seagrass_species as categorical predictor\n")
}
if (include_region) {
  predictor_vars <- c(predictor_vars, "Region")
  cat("  Added Region as categorical predictor\n")
}
if (include_spatial_smoothing) {
  predictor_vars <- c(predictor_vars, "longitude", "latitude")
  cat("  Added longitude and latitude as spatial predictors\n")
}
if (length(exclude_regions) > 0) {
  predictor_vars <- c(predictor_vars, "Region")
  cat("  Excluded regions: ", paste(exclude_regions, collapse = ", "), "\n")
  bathy_grid <- bathy_grid[!bathy_grid$region %in% exclude_regions, ]
}

# Ensure all predictor variables exist in data
missing_vars <- setdiff(predictor_vars, names(raster_covariates_dat))
if (length(missing_vars) > 0) {
  cat("WARNING: The following predictor variables are missing from data and will be excluded:\n")
  cat("  ", paste(missing_vars, collapse = ", "), "\n\n")
  predictor_vars <- setdiff(predictor_vars, missing_vars)
}

# Also check that target_var and coords exist
if (!target_var %in% names(raster_covariates_dat)) {
  stop(
    "Target variable '", target_var, "' not found in data.\n",
    "Available columns include: ", paste(head(names(raster_covariates_dat), 10), collapse = ", "), "..."
  )
}
if (!all(c("longitude", "latitude") %in% names(raster_covariates_dat))) {
  stop("Coordinate columns (longitude, latitude) not found in data.")
}

# Ensure all predictor variables exist in prediction grid
numeric_vars <- predictor_vars[!predictor_vars %in% c("seagrass_species", "Region")]
missing_grid_vars <- setdiff(numeric_vars, names(bathy_grid))
if (length(missing_grid_vars) > 0) {
  warning(
    "The following predictor variables are missing from prediction grid:\n",
    paste(missing_grid_vars, collapse = ", "),
    "\nThese will be excluded from predictions."
  )
  predictor_vars <- setdiff(predictor_vars, missing_grid_vars)
}
# Add 'Unspecified' seagrass_species to prediction grid
bathy_grid$seagrass_species <- "Unspecified"

cat("Final predictor set (", length(predictor_vars), " variables):\n")
cat("  ", paste(predictor_vars, collapse = ", "), "\n\n")

## ================================ FIT GPR ================================
cat("========================================\n")
cat("FITTING GPR MODEL\n")
cat("========================================\n\n")

# Use best kernel from tuning if available (check figures then modelling output dirs)
best_config_path <- "figures/cv_pipeline_output/gpr_best_config.rds"
if (!file.exists(best_config_path)) best_config_path <- "modelling/cv_pipeline_output/gpr_best_config.rds"
if (file.exists(best_config_path)) {
  best_config_list <- readRDS(best_config_path)
  gpr_kernel <- best_config_list$kernel
  nug_min <- best_config_list$nug_min
  nug_max <- best_config_list$nug_max
  cat("Using kernel, nug_min, nug_max from tuning: ", gpr_kernel, ", ", nug_min, ", ", nug_max, "\n")
}

gpr_result <- NULL
if (use_cache_fitted_model && file.exists(model_cache_path)) {
  cached <- tryCatch(readRDS(model_cache_path), error = function(e) NULL)
  if (is_valid_gpr_cache(cached)) {
    gpr_result <- cached
    cat("Loaded fitted model from cache:", model_cache_path, "\n")
    cat("  (Cached model has", gpr_result$n_train, "training points.)\n")
    if (!identical(gpr_result$prediction_grid, bathy_grid)) {
      cat("  Re-predicting on current grid...\n")
      pv <- gpr_result$predictor_vars
      missing_pv <- setdiff(pv, names(bathy_grid))
      if (length(missing_pv) > 0L) {
        warning(
          "Cached model expects predictor(s) missing from grid: ",
          paste(missing_pv, collapse = ", "), ". Fitting new model instead."
        )
        gpr_result <- NULL
      } else {
        pred_complete <- complete.cases(bathy_grid[, pv, drop = FALSE])
        pred_X <- bathy_grid[pred_complete, pv, drop = FALSE]
        pred_X_num <- gpr_apply_onehot(pred_X, gpr_result$encoding, pv)
        XX_scaled <- gpr_apply_scale(pred_X_num, gpr_result$scale_params)
        XX_mat <- as.matrix(XX_scaled[, gpr_result$encoded_names, drop = FALSE])
        storage.mode(XX_mat) <- "double"
        pred_result <- gpr_result$model$pred(XX_mat, se.fit = TRUE, return_df = TRUE)
        bathy_grid$gpr_mean <- NA_real_
        bathy_grid$gpr_se <- NA_real_
        bathy_grid$gpr_mean[pred_complete] <- pred_result$mean
        bathy_grid$gpr_se[pred_complete] <- pred_result$se
        gpr_result$prediction_grid <- bathy_grid
        gpr_result$predictions <- pred_result$mean
        gpr_result$se <- pred_result$se
        gpr_result$n_pred <- sum(pred_complete)
      }
    }
  } else {
    warning(
      "Cached model at ", model_cache_path, " is invalid or empty (missing model/scale_params/predictor_vars). ",
      "Fitting new model."
    )
  }
}
if (is.null(gpr_result)) {
  cat("Fitting GPR model...\n")
  gpr_result <- fit_gaussian_process_regression(
    dat = raster_covariates_dat,
    value_var = target_var,
    coords = c("longitude", "latitude"),
    predictor_vars = predictor_vars,
    formula = NULL,
    prediction_grid = bathy_grid,
    include_spatial = include_spatial_smoothing,
    kernel = gpr_kernel,
    nug_min = nug_min,
    nug_max = nug_max,
  )
  if (use_cache_fitted_model) {
    tryCatch(
      {
        saveRDS(gpr_result, model_cache_path)
        cat("Saved fitted model to cache.\n")
      },
      error = function(e) warning("Could not save model cache: ", e$message)
    )
  }
}

if (is.null(gpr_result) || is.null(gpr_result$model)) {
  stop(
    "GPR model could not be fitted or loaded. ",
    "Check that data and predictor columns are valid."
  )
}

cat("GPR model fitted successfully!\n")
cat("  Training points: ", gpr_result$n_train, "\n")
cat("  Prediction points: ", gpr_result$n_pred, "\n\n")

## ================================ PLOTS ================================
cat("========================================\n")
cat("GENERATING PLOTS\n")
cat("========================================\n\n")

# Filter prediction grid to valid predictions
gpr_pred_map_data <- subset(
  gpr_result$prediction_grid,
  !is.na(gpr_mean) & !is.na(longitude) & !is.na(latitude)
)

cat("Valid predictions for plotting: ", nrow(gpr_pred_map_data), "\n\n")

p_mean_full <- NULL
p_se_full <- NULL
if (nrow(gpr_pred_map_data) == 0L) {
  cat("No valid grid predictions available; skipping spatial GPR plots.\n")
} else {
  # 1. Histogram of predictions
  save_path <- file.path(out_dir, "gpr_prediction_histogram.png")
  cat("1. Histogram of GPR predictions (saving to ", save_path, ")...\n")
  p_hist <- ggplot(data = gpr_pred_map_data, aes(x = gpr_mean)) +
    geom_histogram(binwidth = 0.01, fill = "steelblue", color = "black", alpha = 0.7) +
    labs(x = "GPR prediction", y = "Count", title = if (show_titles) "Histogram of GPR predictions" else NULL) +
    theme_minimal()
  print(p_hist)
  ggsave(save_path, p_hist, width = 10, height = 4, dpi = dpi)

  # 2 & 3. Spatial maps using shared plotting (geom_raster, land on top, obs points)
  lon_range <- c(lon_min, lon_max)
  lat_range <- c(lat_min, lat_max)
  obs_for_map <- raster_covariates_dat[, c("longitude", "latitude", target_var)]
  obs_for_map <- obs_for_map[complete.cases(obs_for_map), ]
  spatial_plots <- plot_gpr_spatial_maps(
    grid = gpr_pred_map_data,
    lon_range = lon_range,
    lat_range = lat_range,
    obs_data = obs_for_map,
    value_col = target_var,
    mean_name = if (target_var %in% names(VAR_LABELS)) VAR_LABELS[target_var] else target_var,
    mean_title = if (show_titles) "GPR Mean Predictions" else NULL,
    se_title = if (show_titles) "GPR Prediction Uncertainty (Standard Error)" else NULL,
    width = 10,
    height = 10,
    dpi = dpi
  )
  p_mean_full <- spatial_plots$p_mean
  p_se_full <- spatial_plots$p_se
  if (!is.null(p_mean_full)) {
    cat("2. Spatial map of GPR mean predictions (saving to ", file.path(out_dir, "gpr_prediction_map.png"), ")...\n")
    print(p_mean_full)
    ggsave(file.path(out_dir, "gpr_prediction_map.png"), p_mean_full, width = 10, height = 10, dpi = dpi)
  }
  if (!is.null(p_se_full)) {
    cat("3. Spatial map of GPR standard error (saving to ", file.path(out_dir, "gpr_se_map.png"), ")...\n")
    print(p_se_full)
    ggsave(file.path(out_dir, "gpr_se_map.png"), p_se_full, width = 10, height = 10, dpi = dpi)
  }
}

# 4. Zoomed maps (if data is concentrated in a smaller area)
data_lon_range <- range(raster_covariates_dat$longitude, na.rm = TRUE)
data_lat_range <- range(raster_covariates_dat$latitude, na.rm = TRUE)
data_lon_span <- diff(data_lon_range)
data_lat_span <- diff(data_lat_range)

if (!is.null(p_mean_full) && !is.null(p_se_full) &&
    (data_lon_span < 0.5 * (lon_max - lon_min) || data_lat_span < 0.5 * (lat_max - lat_min))) {
  save_path_predictions <- file.path(out_dir, "gpr_prediction_map_zoom.png")
  save_path_se <- file.path(out_dir, "gpr_se_map_zoom.png")
  cat("4. Creating zoomed maps (data-focused extent) (saving to ", save_path_predictions, " and ", save_path_se, ")...\n")

  # Add small buffer around data
  lon_buffer <- max(1, 0.1 * data_lon_span)
  lat_buffer <- max(1, 0.1 * data_lat_span)
  zoom_lon_min <- data_lon_range[1] - lon_buffer
  zoom_lon_max <- data_lon_range[2] + lon_buffer
  zoom_lat_min <- data_lat_range[1] - lat_buffer
  zoom_lat_max <- data_lat_range[2] + lat_buffer

  p_mean_zoom <- p_mean_full +
    coord_sf(xlim = c(zoom_lon_min, zoom_lon_max), ylim = c(zoom_lat_min, zoom_lat_max)) +
    labs(title = if (show_titles) "GPR Mean Predictions (Zoomed)" else NULL)
  print(p_mean_zoom)
  ggsave(save_path_predictions, p_mean_zoom, width = 10, height = 10, dpi = dpi)

  p_se_zoom <- p_se_full +
    coord_sf(xlim = c(zoom_lon_min, zoom_lon_max), ylim = c(zoom_lat_min, zoom_lat_max)) +
    labs(title = if (show_titles) "GPR Prediction Uncertainty (Zoomed)" else NULL)
  print(p_se_zoom)
  ggsave(save_path_se, p_se_zoom, width = 10, height = 10, dpi = dpi)
} else {
  cat("4. Skipping zoomed maps (data spans most of full extent)\n")
}

# 5. Combined plot: mean and SE side by side
save_path <- file.path(out_dir, "gpr_mean_and_se.png")
cat("5. Combined plot (mean and SE) (saving to ", save_path, ")...\n")
if (!is.null(p_mean_full) && !is.null(p_se_full) && nrow(gpr_pred_map_data) > 0L) {
  if (requireNamespace("patchwork", quietly = TRUE)) {
    p_combined <- p_mean_full + p_se_full +
      plot_layout(ncol = 2) +
      plot_annotation(title = if (show_titles) "GPR Predictions: Mean and Uncertainty" else NULL)
    print(p_combined)
    ggsave(save_path, p_combined, width = 20, height = 10, dpi = dpi)
  } else {
    cat("  (patchwork not available; skipping combined plot)\n")
  }
} else {
  cat("  (no valid spatial maps available; skipping combined plot)\n")
}

## ================================ PERMUTATION IMPORTANCE ================================
if (compute_permutation_importance) {
  cat("========================================\n")
  cat("COMPUTING PERMUTATION IMPORTANCE\n")
  cat("========================================\n\n")
  pv_imp <- gpr_result$predictor_vars
  train_data_imp <- raster_covariates_dat[, c(target_var, pv_imp), drop = FALSE]
  train_data_imp <- train_data_imp[complete.cases(train_data_imp), ]
  importance_results <- gpr_permutation_importance(
    gpr_result$model, train_data_imp, target_var, pv_imp,
    n_permutations = n_permutations, n_val = 500,
    scale_params = gpr_result$scale_params,
    encoding = gpr_result$encoding,
    encoded_names = gpr_result$encoded_names
  )
  cat("Permutation Importance Results:\n")
  cat("(Higher values = more important; negative values indicate variable may be harmful)\n\n")
  print(importance_results)
  cat("\n")
  write.csv(importance_results,
    file.path(out_dir, "gpr_permutation_importance.csv"),
    row.names = FALSE
  )
  cat("Saved to: gpr_permutation_importance.csv\n\n")
  importance_results$variable_label <- label_vars(importance_results$variable)
  p_importance <- ggplot(importance_results, aes(x = reorder(variable_label, rmse_increase), y = rmse_increase)) +
    geom_col(fill = "steelblue", alpha = 0.85) +
    coord_flip() +
    scale_y_sqrt() +
    labs(x = NULL, y = "RMSE increase (permutation importance; sqrt scale)") +
    theme_minimal(base_size = 12) +
    theme(axis.text.y = element_text(size = rel(0.9)))
  print(p_importance)
  ggsave(file.path(out_dir, "gpr_permutation_importance.png"), p_importance, width = 9, height = max(5, nrow(importance_results) * 0.35), dpi = dpi)
  cv_out <- "figures/cv_pipeline_output"
  if (dir.exists(cv_out)) {
    tryCatch(
      {
        ggsave(file.path(cv_out, "gpr_feature_importance_best_model.png"), p_importance, width = 9, height = max(5, nrow(importance_results) * 0.35), dpi = dpi)
        write.csv(importance_results[, setdiff(names(importance_results), "variable_label")], file.path(cv_out, "gpr_feature_importance_best_model.csv"), row.names = FALSE)
        cat("Also saved to", cv_out, "\n")
      },
      error = function(e) invisible(NULL)
    )
  }
  cat("INTERPRETATION: Permutation importance captures variable effects including interactions.\n\n")
}

## ================================ PARTIAL DEPENDENCE PLOTS ================================
# PDPs are over training data variable ranges (see gpr_funs.R create_partial_dependence).
if (create_partial_dependence_plots) {
  cat("========================================\n")
  cat("GENERATING PARTIAL DEPENDENCE PLOTS (training data ranges)\n")
  cat("========================================\n\n")

  partial_dependence_dir <- file.path(out_dir, "partial_dependence")
  if (!dir.exists(partial_dependence_dir)) {
    dir.create(partial_dependence_dir, recursive = TRUE)
  }
  # Use model's predictor set (includes lon/lat when spatial smoothing enabled)
  pv_pd <- gpr_result$predictor_vars
  # Exclude categorical variables and spatial coordinates from PDP (plot env effects only)
  numeric_predictors <- pv_pd[
    !pv_pd %in% c("seagrass_species", "Region") &
      !pv_pd %in% c("longitude", "latitude")
  ]

  if (is.null(pd_variables)) {
    # Plot all numeric predictors (limit to first 12 for readability)
    vars_to_plot <- head(numeric_predictors, 12)
    if (length(numeric_predictors) > 12) {
      cat("Note: Plotting first 12 of", length(numeric_predictors), "numeric predictors.\n")
      cat("Set pd_variables to specify which variables to plot.\n\n")
    }
  } else {
    # Plot specified variables
    vars_to_plot <- intersect(pd_variables, numeric_predictors)
    if (length(vars_to_plot) < length(pd_variables)) {
      missing <- setdiff(pd_variables, vars_to_plot)
      warning(
        "The following variables were not found or are not numeric: ",
        paste(missing, collapse = ", ")
      )
    }
  }

  if (length(vars_to_plot) == 0) {
    warning("No variables available for partial dependence plots.")
  } else {
    cat("Plotting partial dependence for", length(vars_to_plot), "variables:\n")
    cat("  ", paste(vars_to_plot, collapse = ", "), "\n\n")

    # Prepare training data (complete cases, full predictor set for model compatibility)
    train_data <- raster_covariates_dat[, c(target_var, pv_pd), drop = FALSE]
    train_data <- train_data[complete.cases(train_data), ]

    # Generate plots
    if (length(vars_to_plot) == 1) {
      # Single variable: individual plot
      cat("Creating partial dependence plot for", vars_to_plot[1], "...\n")
      p_pd <- plot_partial_dependence(
        gpr_model = gpr_result$model,
        dat = train_data,
        var_names = vars_to_plot,
        predictor_vars = pv_pd,
        n_points = pd_n_points,
        reference_method = pd_reference_method,
        n_samples = pd_n_samples,
        use_conditional = pd_use_conditional,
        target_var_name = target_var,
        range_quantiles = pd_range_quantiles,
        show_training_range = pd_show_training_range,
        scale_params = gpr_result$scale_params,
        encoding = gpr_result$encoding, encoded_names = gpr_result$encoded_names
      )
      print(p_pd)
      filename <- paste0("gpr_partial_dependence_", gsub("[^A-Za-z0-9]", "_", vars_to_plot[1]), ".png")
      ggsave(file.path(partial_dependence_dir, filename), p_pd, width = 8, height = 6, dpi = dpi)
      cat("  Saved:", filename, "\n")
    } else {
      # Multiple variables: faceted plot
      cat("Creating faceted partial dependence plots...\n")
      p_pd <- plot_partial_dependence(
        gpr_model = gpr_result$model,
        dat = train_data,
        var_names = vars_to_plot,
        predictor_vars = pv_pd,
        n_points = pd_n_points,
        reference_method = pd_reference_method,
        n_samples = pd_n_samples,
        use_conditional = pd_use_conditional,
        target_var_name = target_var,
        range_quantiles = pd_range_quantiles,
        show_training_range = pd_show_training_range,
        scale_params = gpr_result$scale_params,
        encoding = gpr_result$encoding, encoded_names = gpr_result$encoded_names
      )
      print(p_pd)
      ggsave(file.path(partial_dependence_dir, "gpr_partial_dependence_all.png"), p_pd,
        width = 14, height = ceiling(length(vars_to_plot) / 2) * 4, 
        dpi = dpi
      )
      cat("  Saved: gpr_partial_dependence_all.png\n")

      # Also create individual plots for each variable
      cat("\nCreating individual partial dependence plots...\n")
      for (v in vars_to_plot) {
        cat("  Plotting", v, "...\n")
        p_pd_ind <- plot_partial_dependence(
          gpr_model = gpr_result$model,
          dat = train_data,
          var_names = v,
          predictor_vars = pv_pd,
          n_points = pd_n_points,
          reference_method = pd_reference_method,
          n_samples = pd_n_samples,
          use_conditional = pd_use_conditional,
          target_var_name = target_var,
          range_quantiles = pd_range_quantiles,
          show_training_range = pd_show_training_range,
          scale_params = gpr_result$scale_params,
          encoding = gpr_result$encoding, encoded_names = gpr_result$encoded_names
        )
        filename <- paste0("gpr_partial_dependence_", gsub("[^A-Za-z0-9]", "_", v), ".png")
        ggsave(file.path(partial_dependence_dir, filename), p_pd_ind, width = 8, height = 6, dpi = dpi)
      }
      cat("  Individual plots saved.\n")
    }
    cat("\n")
  }
}

cat("\n========================================\n")
cat("GPR PREDICTIONS COMPLETE\n")
cat("========================================\n\n")
cat("Summary:\n")
cat("  Model: Gaussian Process Regression (kernel: ", gpr_kernel, ")\n")
cat("  Predictors: ", length(predictor_vars), " variables\n")
if (include_spatial_smoothing) {
  cat("  Spatial smoothing: Enabled (lon/lat as predictors → predictions revert to nearby mean, not global)\n")
} else {
  cat("  Spatial smoothing: Off (predictions revert to global mean in data-sparse areas)\n")
}
if (include_species) cat("  Species: Included as categorical\n")
if (include_region) cat("  Region: Included as categorical\n")
cat("  Training points: ", gpr_result$n_train, "\n")
cat("  Prediction points: ", gpr_result$n_pred, "\n")
cat("  Output directory: ", out_dir, "\n\n")


cat("Plots saved:\n")
cat("  - ", file.path(out_dir, "gpr_prediction_histogram.png"), "\n")
cat("  - ", file.path(out_dir, "gpr_prediction_map.png"), "\n")
cat("  - ", file.path(out_dir, "gpr_se_map.png"), "\n")
if (exists("p_mean_zoom")) {
  cat("  - ", file.path(out_dir, "gpr_prediction_map_zoom.png"), "\n")
  cat("  - ", file.path(out_dir, "gpr_se_map_zoom.png"), "\n")
}
if (exists("p_combined")) {
  cat("  - ", file.path(out_dir, "gpr_mean_and_se.png"), "\n")
}
if (create_partial_dependence_plots && exists("vars_to_plot") && length(vars_to_plot) > 0) {
  cat("  - ", file.path(out_dir, "gpr_partial_dependence_all.png"), " (faceted)\n")
  if (length(vars_to_plot) > 1) {
    cat("  - ", file.path(out_dir, "gpr_partial_dependence_*.png"), " (individual plots)\n")
  }
}
if (compute_permutation_importance) {
  cat("  - ", file.path(out_dir, "gpr_permutation_importance.png"), "\n")
  cat("  - ", file.path(out_dir, "gpr_permutation_importance.csv"), "\n")
}
# Save prediction grid for later use (e.g. region collages)
if (exists("gpr_result") && !is.null(gpr_result$prediction_grid) && nrow(gpr_result$prediction_grid) > 0L) {
  tryCatch(
    {
      saveRDS(gpr_result$prediction_grid, file.path(out_dir, "gpr_prediction_grid.rds"))
      cat("  - ", file.path(out_dir, "gpr_prediction_grid.rds"), " (for region collages)\n")
    },
    error = function(e) invisible(NULL)
  )
}
# Optional: create region collage (mean/SE panels per region)
if (exists("create_region_collage") && isTRUE(create_region_collage) && exists("gpr_result") && !is.null(gpr_result$prediction_grid)) {
  tryCatch(
    {
      source("modelling/plots/gpr_region_collages.R")
      collage <- create_gpr_region_collage(gpr_result$prediction_grid)
      ggsave(file.path(out_dir, "gpr_region_collage.png"), collage, width = 14, height = 10, dpi = dpi)
      cat("  - ", file.path(out_dir, "gpr_region_collage.png"), "\n")
    },
    error = function(e) warning("Region collage failed: ", e$message)
  )
}
