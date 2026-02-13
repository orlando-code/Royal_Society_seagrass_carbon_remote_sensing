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
# Ensure extractor helpers are defined in the global environment so COVARIATE_CONFIG is visible
sys.source("modelling/R/extract_covariates_from_rasters.R", envir = .GlobalEnv)
if (!exists("COVARIATE_CONFIG", inherits = TRUE)) {
  # Rebuild if not created during sourcing (defensive against environment quirks)
  if (!exists("RASTER_DIR", inherits = TRUE)) {
    RASTER_DIR <- "data/env_rasters"
  }
  COVARIATE_CONFIG <<- build_covariate_config_from_dir(RASTER_DIR)
}
source("modelling/R/gpr_funs.R")
source("modelling/R/assign_region_from_latlon.R")
source("modelling/R/plot_config.R")
load_packages(c("here", "mgcv", "tidyverse", "ggplot2", "GauPro", "viridisLite", "maps"))

## ================================ CONFIGURATION ================================
# Target variable
target_var <- "median_carbon_density_100cm"

# Prediction grid resolution (use 5000 for full paper; 100 for quick test)
n_lons <- 100
n_lats <- 100

# Options
# - Species: keep grid predictions species-agnostic (include_species = FALSE).
# - Region: CV showed Region hurts model fit; keep FALSE so predictions do not use region.
# - Spatial smoothing (lat/lon): keep TRUE so predictions revert to the mean of *nearby* data
#   rather than the global mean. Without lat/lon, the GPR would revert to the global dataset
#   mean in data-sparse or extrapolation areas. With lat/lon, we get smooth spatial variation
#   (local reversion). This is independent of CV: we can report best CV from env-only or
#   env+lat/lon, but for maps we want spatially informed predictions.
include_species <- FALSE            # Keep grid predictions species-agnostic
include_region <- FALSE             # Region hurt CV fit; omit for prediction
include_spatial_smoothing <- TRUE   # Lat/lon as predictors → revert to nearby mean (not global)

# Partial dependence plots
create_partial_dependence_plots <- TRUE  # Generate partial dependence plots
pd_variables <- NULL  # NULL = plot all numeric predictors, or specify vector of variable names
pd_n_points <- 100     # Number of evaluation points per variable
pd_reference_method <- "conditional"  # "median", "mean", "sample", or "conditional" (uses actual data bins)
pd_n_samples <- 1000   # If pd_reference_method = "sample", number of samples to average over
pd_use_conditional <- TRUE  # If TRUE, use conditional partial dependence (only actual data combinations)

# Variable importance analysis
compute_permutation_importance <- TRUE  # Compute permutation importance (shuffle variable, measure performance drop)
n_permutations <- 5  # Number of permutations per variable (for stability)

# Kernel: use best from tuning if available
gpr_kernel <- "matern52"  # Overridden by gpr_best_config.rds if present
use_cache_fitted_model <- FALSE  # If TRUE and cache exists, load fitted model instead of re-fitting
model_cache_path <- "figures/cv_pipeline_output/gpr_fitted_model.rds"

# PDPs: explicitly over training data variable ranges (optionally restrict to percentile range)
pd_range_quantiles <- c(0.01, 0.99)  # e.g. c(0.01, 0.99) to avoid extreme tails
pd_show_training_range <- TRUE  # If TRUE, add rug of training values on PDPs

# Output directory
out_dir <- "figures/interpolation"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
# Optional: create and save region collage (spatial prediction panels per region)
create_region_collage <- FALSE  # Set TRUE to save gpr_region_collage.png; or run via run_paper_figures.R

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
# then GPR-specific list, then GAM list
gam_pruned_vars <- get_pruned_predictors_for_model("GPR", data_cols = names(raster_covariates_dat))
if (!is.null(gam_pruned_vars) && length(gam_pruned_vars) >= 2L) {
  cat("Loaded", length(gam_pruned_vars), "GPR covariates from model permutation pruning (1 km spatial CV).\n")
} else {
  gam_pruned_vars <- NULL
}
if (is.null(gam_pruned_vars) || length(gam_pruned_vars) < 2L) {
  pruned_file_gpr <- "figures/covariate_selection/pruned_variables_to_include_gpr.csv"
  pruned_file <- "figures/covariate_selection/pruned_variables_to_include.csv"
  if (file.exists(pruned_file_gpr)) {
    gam_pruned_vars <- read_csv(pruned_file_gpr, show_col_types = FALSE)$variable
    cat("Loaded", length(gam_pruned_vars), "GPR-pruned covariates (gpr_covariate_pruning fallback).\n")
  } else if (file.exists(pruned_file)) {
    gam_pruned_vars <- read_csv(pruned_file, show_col_types = FALSE)$variable
    cat("Loaded", length(gam_pruned_vars), "pruned covariates (GAM list fallback).\n")
  } else {
    stop("Pruned variables not found. Run model_permutation_pruning.R, or covariate_pruning.R / gpr_covariate_pruning.R.")
  }
}
gam_pruned_vars <- gam_pruned_vars[gam_pruned_vars %in% names(raster_covariates_dat)]
cat("  ", paste(gam_pruned_vars, collapse = ", "), "\n\n")

# plot histogram of pruned variables in data on a facet plot
p_histograms <- gam_pruned_vars |>
  purrr::map_dfr(function(var) {
    data.frame(
      value = raster_covariates_dat[[var]],
      variable = var
    )
  }) |>
  ggplot(aes(x = value)) +
    geom_histogram(binwidth = 0.01, fill = "steelblue", color = "black", alpha = 0.7) +
    facet_wrap(~ variable, scales = "free") +
    labs(x = "Value", y = "Count") +
    theme_minimal()
print(p_histograms)
ggsave(file.path(out_dir, "covariate_distributions.png"), p_histograms, width = 14, height = 10)

cat("\n========================================\n")
cat("COVARIATE DISTRIBUTION ANALYSIS\n")
cat("========================================\n\n")
cat("SPARSE DISTRIBUTIONS DETECTED:\n")
cat("----------------------------------------\n")
cat("Many covariates show sparse, gapped, or discrete distributions:\n\n")

# Analyze each variable's distribution
sparse_vars <- character()
continuous_vars <- character()

for (var in gam_pruned_vars) {
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
        gap_ratio <- max_gap / (median_gap + 1e-10)  # Avoid division by zero
        
        # Classify as sparse if: many unique values relative to total, or large gaps
        if (gap_ratio > 10 || unique_ratio > 0.5) {
          sparse_vars <- c(sparse_vars, var)
          cat("  SPARSE: ", var, "\n")
          cat("    - Unique values: ", n_unique, " / ", n_total, 
              " (", round(100 * unique_ratio, 1), "%)\n", sep = "")
          cat("    - Max gap ratio: ", round(gap_ratio, 2), "\n")
        } else {
          continuous_vars <- c(continuous_vars, var)
        }
      }
    }
  }
}

cat("\nCONTINUOUS DISTRIBUTIONS:\n")
cat("----------------------------------------\n")
if (length(continuous_vars) > 0) {
  cat("  ", paste(continuous_vars, collapse = ", "), "\n")
} else {
  cat("  (None detected)\n")
}

cat("\nIMPLICATIONS FOR PARTIAL DEPENDENCE PLOTS:\n")
cat("----------------------------------------\n")
cat("1. SPARSE VARIABLES WILL SHOW JAGGED PDPs:\n")
cat("   - Large gaps = model uncertainty in those regions\n")
cat("   - GPR interpolates/extrapolates in gaps\n")
cat("   - Predictions can vary wildly between empty regions\n")
cat("   - This is EXPECTED and not necessarily overfitting\n\n")

cat("2. OVERFITTING CONCERNS:\n")
cat("   - If R² is good on spatial CV, model generalizes well\n")
cat("   - Jagged PDPs may reflect true conditional effects\n")
cat("   - Use permutation importance to verify variable importance\n")
cat("   - Consider smoother kernels (matern32 vs matern52)\n\n")

cat("3. RECOMMENDATIONS:\n")
cat("   - Focus interpretation on well-sampled ranges\n")
cat("   - Use smoothing (LOESS) for visualization\n")
cat("   - Acknowledge uncertainty in sparse regions\n")
cat("   - Run hyperparameter tuning (see gpr_hyperparameter_tuning.R)\n\n")

cat("4. TUNING FOR REGULARIZATION:\n")
cat("   - Test different kernels: matern52, matern32, gaussian\n")
cat("   - Matern32 is smoother (more regularization)\n")
cat("   - Use 1km spatial CV to select best kernel\n")
cat("   - See: modelling/gpr_hyperparameter_tuning.R\n\n")

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

# Get raster covariates from COVARIATE_CONFIG (built in extract_covariates_from_rasters.R)
RASTER_DIR <- "data/env_rasters"
# (Re)build a fresh covariate configuration here to avoid relying on
# side-effects from sourcing in different environments.
source("modelling/R/extract_covariates_from_rasters.R") # TODO: this is a bit of a hack to ensure build_covariate_config_from_dir is found 
COVARIATE_CONFIG <- build_covariate_config_from_dir(RASTER_DIR)
# raster_covariates <- c(bathy_covariate, gam_pruned_vars)
raster_covariates <- tolower(names(COVARIATE_CONFIG))
bathy_covariate <- grep("bathy|gebco|static_depth", raster_covariates, ignore.case = TRUE, value = TRUE)
# Use first pruned covariate as reference raster (or find one that exists)
reference_raster <- gam_pruned_vars[1]
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
  reference_raster = reference_raster,  # TODO: does the reference raster matter? e.g. its resolution?
  cache_path = sprintf("data/prediction_grid_cache_lons%s_lats%s.rds", n_lons, n_lats),
  use_cache = TRUE,  # Use cache if available
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
        sprintf("Bathymetry column '%s' not found in raster_covariates_dat. Names are: %s",
                bathy_covariate, paste(names(raster_covariates_dat), collapse = ", "))
      )
    }
  }
  bathy_95perc <- quantile(bathy_vals, 0.99, na.rm = TRUE)
  bathy_max <- 2 * bathy_95perc
  bathy_grid <- prediction_grid[
    prediction_grid[[bathy_covariate]] <= bathy_max & 
    prediction_grid[[bathy_covariate]] >= 0, 
  ]
  cat("Filtered to bathymetry <= ", bathy_max, "m: ", nrow(bathy_grid), " cells\n\n")
} else {
  bathy_grid <- prediction_grid
  cat("No bathymetry column found; using full prediction grid\n\n")
}
bathy_grid <- process_rs_covariates(bathy_grid)
# Tag prediction grid with inferred region (for plotting/zoomed collages)
if (!"region" %in% names(bathy_grid) && all(c("longitude", "latitude") %in% names(bathy_grid))) {
  bathy_grid <- assign_region_from_latlon(bathy_grid)
}

## ================================ GPR FIT (see gpr_funs.R) ================================
## ================================ BUILD PREDICTOR SET ================================
cat("========================================\n")
cat("BUILDING PREDICTOR SET\n")
cat("========================================\n\n")

predictor_vars <- gam_pruned_vars

# Add categorical variables if requested
if (include_species) {
  predictor_vars <- c(predictor_vars, "seagrass_species")
  cat("  Added seagrass_species as categorical predictor\n")
}
if (include_region) {
  predictor_vars <- c(predictor_vars, "Region")
  cat("  Added Region as categorical predictor\n")
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
  stop("Target variable '", target_var, "' not found in data.\n",
       "Available columns include: ", paste(head(names(raster_covariates_dat), 10), collapse = ", "), "...")
}
if (!all(c("longitude", "latitude") %in% names(raster_covariates_dat))) {
  stop("Coordinate columns (longitude, latitude) not found in data.")
}

# Ensure all predictor variables exist in prediction grid (except categoricals)
numeric_vars <- predictor_vars[!predictor_vars %in% c("seagrass_species", "Region")]
missing_grid_vars <- setdiff(numeric_vars, names(bathy_grid))
if (length(missing_grid_vars) > 0) {
  warning("The following predictor variables are missing from prediction grid:\n",
          paste(missing_grid_vars, collapse = ", "),
          "\nThese will be excluded from predictions.")
  predictor_vars <- setdiff(predictor_vars, missing_grid_vars)
}

cat("Final predictor set (", length(predictor_vars), " variables):\n")
cat("  ", paste(predictor_vars, collapse = ", "), "\n\n")

# Create formula
pruned_formula_reg <- as.formula(paste(target_var, "~", paste(predictor_vars, collapse = " + ")))

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
  cat("Using kernel from tuning: ", gpr_kernel, "\n")
}

gpr_result <- NULL
if (use_cache_fitted_model && file.exists(model_cache_path)) {
  cat("Loading fitted model from cache:", model_cache_path, "\n")
  gpr_result <- readRDS(model_cache_path)
  cat("  (Cached model has", gpr_result$n_train, "training points.)\n")
  if (!identical(gpr_result$prediction_grid, bathy_grid)) {
    cat("  Re-predicting on current grid...\n")
    pred_complete <- complete.cases(bathy_grid[, predictor_vars, drop = FALSE])
    XX <- as.matrix(bathy_grid[pred_complete, predictor_vars, drop = FALSE])
    if (!is.null(gpr_result$scale_params)) XX <- gpr_apply_scale(XX, gpr_result$scale_params)
    pred_result <- gpr_result$model$pred(XX, se.fit = TRUE, return_df = TRUE)
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
source("modelling/R/gpr_funs.R")
if (is.null(gpr_result)) {
  gpr_result <- fit_gaussian_process_regression(
    dat = raster_covariates_dat,
    value_var = target_var,
    coords = c("longitude", "latitude"),
    predictor_vars = predictor_vars,
    formula = pruned_formula_reg,
    prediction_grid = bathy_grid,
    include_spatial = include_spatial_smoothing,
    kernel = gpr_kernel
  )
  if (use_cache_fitted_model) {
    tryCatch({
      saveRDS(gpr_result, model_cache_path)
      cat("Saved fitted model to cache.\n")
    }, error = function(e) warning("Could not save model cache: ", e$message))
  }
}

cat("GPR model fitted successfully!\n")
cat("  Training points: ", gpr_result$n_train, "\n")
cat("  Prediction points: ", gpr_result$n_pred, "\n\n")

## ================================ PLOTS ================================
cat("========================================\n")
cat("GENERATING PLOTS\n")
cat("========================================\n\n")

# World map background
world <- map_data("world")

# Filter prediction grid to valid predictions
gpr_pred_map_data <- subset(
  gpr_result$prediction_grid,
  !is.na(gpr_mean) & !is.na(longitude) & !is.na(latitude)
)

cat("Valid predictions for plotting: ", nrow(gpr_pred_map_data), "\n\n")

if (nrow(gpr_pred_map_data) == 0L) {
  cat("No valid grid predictions available; skipping spatial GPR plots.\n")
} else {
  # 1. Histogram of predictions
  save_path <- file.path(out_dir, "gpr_prediction_histogram.png")
  cat("1. Histogram of GPR predictions (saving to ", save_path, ")...\n")
  p_hist <- ggplot(data = gpr_pred_map_data, aes(x = gpr_mean)) +
    geom_histogram(binwidth = 0.01, fill = "steelblue", color = "black", alpha = 0.7) +
    labs(x = "GPR prediction", y = "Count", title = "Histogram of GPR predictions") +
    scale_y_log10() +
    theme_minimal()
  print(p_hist)
  ggsave(save_path, p_hist, width = 10, height = 4)

  gpr_pred_map_data_sub <- gpr_pred_map_data[seq(1, nrow(gpr_pred_map_data), 100), ]
  # 2. Spatial map of mean predictions (full extent)
  save_path <- file.path(out_dir, "gpr_prediction_map.png")
  cat("2. Spatial map of GPR mean predictions (full extent) (saving to ", save_path, ")...\n")
  # subset to every 100th point
  p_mean_full <- ggplot() +
  geom_polygon(
    data = world,
    aes(x = long, y = lat, group = group),
    fill = "grey95",
    color = "grey60",
    linewidth = 0.2
  ) +
  geom_point(
    data = gpr_pred_map_data_sub,
    aes(x = longitude, y = latitude, color = gpr_mean),
    size = 0.1,
    alpha = 0.85
  ) +
  scale_color_gradientn(
    colors = viridisLite::turbo(256),
    name = target_var
  ) +
  coord_cartesian(
    xlim = c(lon_min, lon_max),
    ylim = c(lat_min, lat_max),
    expand = FALSE
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(2, "cm"),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid = element_line(color = "grey90")
  ) +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "GPR Mean Predictions"
  )
  print(p_mean_full)
  ggsave(save_path, p_mean_full, width = 10, height = 10)

  # 3. Spatial map of standard error (full extent)
  save_path <- file.path(out_dir, "gpr_se_map.png")
  cat("3. Spatial map of GPR standard error (full extent) (saving to ", save_path, ")...\n")
  p_se_full <- ggplot() +
  geom_polygon(
    data = world,
    aes(x = long, y = lat, group = group),
    fill = "grey95",
    color = "grey60",
    linewidth = 0.2
  ) +
  geom_point(
    data = gpr_pred_map_data_sub,
    aes(x = longitude, y = latitude, color = gpr_se),
    size = 0.1,
    alpha = 0.85
  ) +
  scale_color_gradientn(
    colors = viridisLite::viridis(256),
    name = "Standard Error"
  ) +
  coord_cartesian(
    xlim = c(lon_min, lon_max),
    ylim = c(lat_min, lat_max),
    expand = FALSE
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(2, "cm"),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid = element_line(color = "grey90")
  ) +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "GPR Prediction Uncertainty (Standard Error)"
  )
  print(p_se_full)
  ggsave(save_path, p_se_full, width = 10, height = 10)
}

# 4. Zoomed maps (if data is concentrated in a smaller area)
# Calculate data extent
data_lon_range <- range(raster_covariates_dat$longitude, na.rm = TRUE)
data_lat_range <- range(raster_covariates_dat$latitude, na.rm = TRUE)
data_lon_span <- diff(data_lon_range)
data_lat_span <- diff(data_lat_range)

# If data spans less than 50% of full extent, create zoomed plots
if (data_lon_span < 0.5 * (lon_max - lon_min) || data_lat_span < 0.5 * (lat_max - lat_min)) {
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
    coord_cartesian(
      xlim = c(zoom_lon_min, zoom_lon_max),
      ylim = c(zoom_lat_min, zoom_lat_max),
      expand = FALSE
    ) +
    labs(title = "GPR Mean Predictions (Zoomed)")
  print(p_mean_zoom)
  ggsave(save_path, p_mean_zoom, width = 10, height = 10)
  
  p_se_zoom <- p_se_full +
    coord_cartesian(
      xlim = c(zoom_lon_min, zoom_lon_max),
      ylim = c(zoom_lat_min, zoom_lat_max),
      expand = FALSE
    ) +
    labs(title = "GPR Prediction Uncertainty (Zoomed)")
  print(p_se_zoom)
  ggsave(save_path_se, p_se_zoom, width = 10, height = 10)
} else {
  cat("4. Skipping zoomed maps (data spans most of full extent)\n")
}

# 5. Combined plot: mean and SE side by side
save_path <- file.path(out_dir, "gpr_mean_and_se.png")
cat("5. Combined plot (mean and SE) (saving to ", save_path, ")...\n")
if (exists("p_mean_full") && exists("p_se_full") && nrow(gpr_pred_map_data) > 0L) {
  cat("5. Combined plot (mean and SE) (saving to ", save_path, ")...\n")
  if (requireNamespace("patchwork", quietly = TRUE)) {
    p_combined <- p_mean_full + p_se_full +
      plot_layout(ncol = 2) +
      plot_annotation(title = "GPR Predictions: Mean and Uncertainty")
    print(p_combined)
    ggsave(save_path, p_combined, width = 20, height = 10)
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
  train_data_imp <- raster_covariates_dat[, c(target_var, predictor_vars), drop = FALSE]
  train_data_imp <- train_data_imp[complete.cases(train_data_imp), ]
  importance_results <- gpr_permutation_importance(
    gpr_result$model, train_data_imp, target_var, predictor_vars,
    n_permutations = n_permutations, n_val = 500,
    scale_params = gpr_result$scale_params
  )
  cat("Permutation Importance Results:\n")
  cat("(Higher values = more important; negative values indicate variable may be harmful)\n\n")
  print(importance_results)
  cat("\n")
  write.csv(importance_results,
            file.path(out_dir, "gpr_permutation_importance.csv"),
            row.names = FALSE)
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
  ggsave(file.path(out_dir, "gpr_permutation_importance.png"), p_importance, width = 9, height = max(5, nrow(importance_results) * 0.35), dpi = 150)
  cv_out <- "figures/cv_pipeline_output"
  if (dir.exists(cv_out)) {
    tryCatch({
      ggsave(file.path(cv_out, "gpr_feature_importance_best_model.png"), p_importance, width = 9, height = max(5, nrow(importance_results) * 0.35), dpi = 150)
      write.csv(importance_results[, setdiff(names(importance_results), "variable_label")], file.path(cv_out, "gpr_feature_importance_best_model.csv"), row.names = FALSE)
      cat("Also saved to", cv_out, "\n")
    }, error = function(e) invisible(NULL))
  }
  cat("INTERPRETATION: Permutation importance captures variable effects including interactions.\n\n")
}

## ================================ PARTIAL DEPENDENCE PLOTS ================================
# PDPs are over training data variable ranges (see gpr_funs.R create_partial_dependence).
if (create_partial_dependence_plots) {
  cat("========================================\n")
  cat("GENERATING PARTIAL DEPENDENCE PLOTS (training data ranges)\n")
  cat("========================================\n\n")
  
  # Determine which variables to plot
  # Exclude categorical variables and spatial coordinates (unless spatial smoothing is enabled)
  numeric_predictors <- predictor_vars[
    !predictor_vars %in% c("seagrass_species", "Region") &
    (!predictor_vars %in% c("longitude", "latitude") | include_spatial_smoothing)
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
      warning("The following variables were not found or are not numeric: ",
              paste(missing, collapse = ", "))
    }
  }
  
  if (length(vars_to_plot) == 0) {
    warning("No variables available for partial dependence plots.")
  } else {
    cat("Plotting partial dependence for", length(vars_to_plot), "variables:\n")
    cat("  ", paste(vars_to_plot, collapse = ", "), "\n\n")
    
    # Prepare training data (complete cases only)
    train_data <- raster_covariates_dat[, c(target_var, predictor_vars), drop = FALSE]
    train_data <- train_data[complete.cases(train_data), ]
    
    # Generate plots
    if (length(vars_to_plot) == 1) {
      # Single variable: individual plot
      cat("Creating partial dependence plot for", vars_to_plot[1], "...\n")
      p_pd <- plot_partial_dependence(
        gpr_model = gpr_result$model,
        dat = train_data,
        var_names = vars_to_plot,
        predictor_vars = predictor_vars,
        n_points = pd_n_points,
        reference_method = pd_reference_method,
        n_samples = pd_n_samples,
        use_conditional = pd_use_conditional,
        target_var_name = target_var,
        range_quantiles = pd_range_quantiles,
        show_training_range = pd_show_training_range,
        scale_params = gpr_result$scale_params
      )
      print(p_pd)
      filename <- paste0("gpr_partial_dependence_", gsub("[^A-Za-z0-9]", "_", vars_to_plot[1]), ".png")
      ggsave(file.path(out_dir, filename), p_pd, width = 8, height = 6)
      cat("  Saved:", filename, "\n")
    } else {
      # Multiple variables: faceted plot
      cat("Creating faceted partial dependence plots...\n")
      p_pd <- plot_partial_dependence(
        gpr_model = gpr_result$model,
        dat = train_data,
        var_names = vars_to_plot,
        predictor_vars = predictor_vars,
        n_points = pd_n_points,
        reference_method = pd_reference_method,
        n_samples = pd_n_samples,
        use_conditional = pd_use_conditional,
        target_var_name = target_var,
        range_quantiles = pd_range_quantiles,
        show_training_range = pd_show_training_range,
        scale_params = gpr_result$scale_params
      )
      print(p_pd)
      ggsave(file.path(out_dir, "gpr_partial_dependence_all.png"), p_pd, 
             width = 14, height = ceiling(length(vars_to_plot) / 2) * 4)
      cat("  Saved: gpr_partial_dependence_all.png\n")
      
      # Also create individual plots for each variable
      cat("\nCreating individual partial dependence plots...\n")
      for (v in vars_to_plot) {
        cat("  Plotting", v, "...\n")
        p_pd_ind <- plot_partial_dependence(
          gpr_model = gpr_result$model,
          dat = train_data,
          var_names = v,
          predictor_vars = predictor_vars,
          n_points = pd_n_points,
          reference_method = pd_reference_method,
          n_samples = pd_n_samples,
          use_conditional = pd_use_conditional,
          target_var_name = target_var,
          range_quantiles = pd_range_quantiles,
          show_training_range = pd_show_training_range,
          scale_params = gpr_result$scale_params
        )
        filename <- paste0("gpr_partial_dependence_", gsub("[^A-Za-z0-9]", "_", v), ".png")
        ggsave(file.path(out_dir, filename), p_pd_ind, width = 8, height = 6)
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

if (create_partial_dependence_plots && compute_permutation_importance) {
  cat("========================================\n")
  cat("UNDERSTANDING FLAT PARTIAL DEPENDENCE\n")
  cat("========================================\n\n")
  cat("If your partial dependence plots are flat but the model performs well:\n\n")
  cat("1. STRONG INTERACTIONS:\n")
  cat("   - Variables may have important effects only when combined with others\n")
  cat("   - Partial dependence averages over other variables, masking interactions\n")
  cat("   - Check permutation importance - variables with high importance but flat PDP\n")
  cat("     likely have strong interactions\n\n")
  cat("2. SPATIAL STRUCTURE DOMINATES:\n")
  cat("   - If spatial coordinates (lon/lat) are included, they may capture most signal\n")
  cat("   - Environmental variables may add refinement but show little marginal effect\n")
  cat("   - This is common in spatial models\n\n")
  cat("3. CONDITIONAL EFFECTS:\n")
  cat("   - A variable's effect may depend strongly on values of other variables\n")
  cat("   - Averaging over reference values smooths out these conditional effects\n")
  cat("   - Try ICE plots or 2D partial dependence to see interactions\n\n")
  cat("4. NON-LINEAR COMBINATIONS:\n")
  cat("   - GPR can capture complex non-linear relationships\n")
  cat("   - Individual variables may have little effect alone but matter in combination\n")
  cat("   - Permutation importance captures this better than partial dependence\n\n")
  cat("RECOMMENDATION: Use permutation importance to identify truly important variables.\n")
  cat("Partial dependence is useful for understanding marginal effects, but may miss\n")
  cat("variables that are important through interactions.\n\n")
}

cat("Plots saved:\n")
cat("  - gpr_prediction_histogram.png\n")
cat("  - gpr_prediction_map.png\n")
cat("  - gpr_se_map.png\n")
if (exists("p_mean_zoom")) {
  cat("  - gpr_prediction_map_zoom.png\n")
  cat("  - gpr_se_map_zoom.png\n")
}
if (exists("p_combined")) {
  cat("  - gpr_mean_and_se.png\n")
}
if (create_partial_dependence_plots && exists("vars_to_plot") && length(vars_to_plot) > 0) {
  cat("  - gpr_partial_dependence_all.png (faceted)\n")
  if (length(vars_to_plot) > 1) {
    cat("  - gpr_partial_dependence_*.png (individual plots)\n")
  }
}
if (compute_permutation_importance) {
  cat("  - gpr_permutation_importance.png\n")
  cat("  - gpr_permutation_importance.csv\n")
}
# Save prediction grid for later use (e.g. region collages)
if (exists("gpr_result") && !is.null(gpr_result$prediction_grid) && nrow(gpr_result$prediction_grid) > 0L) {
  tryCatch({
    saveRDS(gpr_result$prediction_grid, file.path(out_dir, "gpr_prediction_grid.rds"))
    cat("  - gpr_prediction_grid.rds (for region collages)\n")
  }, error = function(e) invisible(NULL))
}
# Optional: create region collage (mean/SE panels per region)
if (exists("create_region_collage") && isTRUE(create_region_collage) && exists("gpr_result") && !is.null(gpr_result$prediction_grid)) {
  tryCatch({
    source("modelling/plots/gpr_region_collages.R")
    collage <- create_gpr_region_collage(gpr_result$prediction_grid)
    ggsave(file.path(out_dir, "gpr_region_collage.png"), collage, width = 14, height = 10, dpi = 150)
    cat("  - gpr_region_collage.png\n")
  }, error = function(e) warning("Region collage failed: ", e$message))
}