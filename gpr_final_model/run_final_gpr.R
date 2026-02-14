# =============================================================================
# Final GPR model: fit, train-test evaluation, and N by N grid predictions
#
# Run from project root:  source("gpr_final_model/run_final_gpr.R")
# Or:  Rscript gpr_final_model/run_final_gpr.R
#
# Produces:
#   - gpr_final_model.rds       Fitted model for load + predict_gpr()
#   - train_test_metrics.csv    R² and RMSE from random train-test split
#   - grid_predictions.rds       Mean and SE on N by N covariate grid (optional)
#   - grid_plot_mean.png, grid_plot_se.png (optional)
# =============================================================================

if (!requireNamespace("here", quietly = TRUE)) install.packages("here", quiet = TRUE)
setwd(here::here())
set.seed(42)

# -----------------------------------------------------------------------------
# Options
# -----------------------------------------------------------------------------
train_frac <- 0.8 # Fraction of data for training in train-test split
n_grid_side <- 1000L # Points per axis for N by N spatial grid

if (n_grid_side > 1000L) {
  stop("Your machine will likely run out of memory with a grid size of ", n_grid_side, " by ", n_grid_side, ". Please reduce the grid size.")
}

save_grid_plots <- TRUE
mask_to_bathymetry <- TRUE # Mask to bathymetry < 40m
bathymetry_threshold <- 40
out_dir <- "gpr_final_model"
figures_dir <- file.path(out_dir, "figures")
spatial_cache_dir <- file.path(out_dir, "spatial_grid_cache")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(spatial_cache_dir, recursive = TRUE, showWarnings = FALSE)
# -----------------------------------------------------------------------------
# Dependencies and helpers
# -----------------------------------------------------------------------------
source("modelling/R/helpers.R")
source("modelling/R/gpr_funs.R")
# For spatial grid: covariate extraction from rasters (builds raster_covariates)
source("modelling/R/extract_covariates_from_rasters.R")
source("modelling/plots/plot_gpr_spatial_maps.R")
load_packages(c("here", "readr", "dplyr", "GauPro", "ggplot2", "viridisLite"))


# -----------------------------------------------------------------------------
# Data
# -----------------------------------------------------------------------------
target_var <- "median_carbon_density_100cm"
if (file.exists("data/all_extracted_new.rds")) {
  dat_raw <- readRDS("data/all_extracted_new.rds")
  dat <- dat_raw[dat_raw$region != "Black Sea", ] # Remove Black Sea values (which drags GPR kernel weirdly into the Eastern Med)
  dat <- process_rs_covariates(dat)
} else {
  stop("Data file `data/all_extracted_new.rds` not found. Run main pipeline (`modelling/run_paper_figures.R`) first.")
}


# Best kernel from tuning. TODO: also pending further feature pruning refinement
best_config_path <- "figures/cv_pipeline_output/gpr_best_config.rds"
if (file.exists(best_config_path)) {
  best_config <- readRDS(best_config_path)
  gpr_kernel <- best_config$kernel
  predictor_vars <- best_config$features_used
  gpr_min_nug <- best_config$nug_min
  gpr_max_nug <- best_config$nug_max
} else {
  # Predictor set: from model permutation pruning CSV or fallback
  predictor_vars <- NULL
  if (is.null(predictor_vars) || length(predictor_vars) < 2L) {
    pruned_file <- "figures/covariate_selection/pruned_variables_to_include_gpr.csv" # TODO: this is a placeholder pending further feature pruning refinement
    if (file.exists(pruned_file)) {
      predictor_vars <- readr::read_csv(pruned_file, show_col_types = FALSE)$variable
      predictor_vars <- predictor_vars[predictor_vars %in% colnames(dat)]
    }
  }
  if (is.null(predictor_vars) || length(predictor_vars) < 2L) {
    stop("GPR predictor set not found. Run model_permutation_pruning.R (or the full pipeline) first.")
  }}

# Add spatial and species predictors
predictor_vars <- unique(c("seagrass_species", predictor_vars))
# Add lon/lat for spatial smoothing
coords <- c("longitude", "latitude")
predictor_vars <- unique(c(coords, predictor_vars))
predictor_vars <- predictor_vars[predictor_vars %in% colnames(dat) | predictor_vars == "seagrass_species"]

need_vars <- c(target_var, predictor_vars)
dat <- dat[, intersect(names(dat), need_vars)]
dat <- dat[complete.cases(dat), ]


cat("Kernel:", gpr_kernel, "\n")
cat("Nugget:", gpr_min_nug, "to", gpr_max_nug, "\n")
cat("Predictors:", length(predictor_vars), "\n\n")

# -----------------------------------------------------------------------------
# Train-test split evaluation
# -----------------------------------------------------------------------------
n <- nrow(dat)
idx_train <- sample.int(n, size = round(train_frac * n))
train_data <- dat[idx_train, ]
test_data <- dat[-idx_train, ]

cat("Train-test split: ", nrow(train_data), " train, ", nrow(test_data), " test\n")

fit_train <- fit_gaussian_process_regression(
  dat             = train_data,
  value_var       = target_var,
  coords          = coords,
  predictor_vars  = predictor_vars,
  formula         = NULL,
  prediction_grid = test_data,
  include_spatial = TRUE,
  kernel          = gpr_kernel,
  nug_min         = gpr_min_nug,
  nug_max         = gpr_max_nug,
)

test_pred <- predict_gpr(fit_train, test_data[, predictor_vars, drop = FALSE], se = TRUE)
#  N.B. this is an optimistic evaluation due to only random splitting (potentially very small spatial separation between train and test)
obs <- test_data[[target_var]]
r2 <- cor(obs, test_pred$mean)^2
rmse <- sqrt(mean((obs - test_pred$mean)^2))

cat("Test R²:", round(r2, 4), "\n")
cat("Test RMSE:", round(rmse, 6), "\n\n")

write.csv(
  data.frame(metric = c("r2", "rmse"), value = c(r2, rmse)),
  file.path(out_dir, "train_test_metrics.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# Fit on full data and build N by N spatial prediction grid
# -----------------------------------------------------------------------------
# N by N grid over the spatial range of the dataset, with environmental covariates
# extracted from rasters at each (lon, lat) -> spatial prediction map
env_covariates <- setdiff(predictor_vars, c("longitude", "latitude"))
env_covariates <- intersect(env_covariates, raster_covariates)
if (length(env_covariates) == 0L) {
  stop("No predictor variables found in raster_covariates. Check raster config and predictor set.")
}
lon_range <- range(dat$longitude, na.rm = TRUE)
lat_range <- range(dat$latitude, na.rm = TRUE)
grid_cache_path <- file.path(spatial_cache_dir, sprintf("spatial_grid_n%d.rds", n_grid_side))

cat("Building ", n_grid_side, " by ", n_grid_side, " spatial grid over data extent...\n", sep = "")
reference_raster <- env_covariates[1L]
spatial_grid <- create_prediction_grid_from_rasters(
  lon_range        = lon_range,
  lat_range        = lat_range,
  covariates       = env_covariates,
  n_lon            = n_grid_side,
  n_lat            = n_grid_side,
  reference_raster = reference_raster,
  cache_path       = grid_cache_path,
  use_cache        = FALSE,
  show_progress    = TRUE,
  fill_coastal     = TRUE
)
spatial_grid <- process_rs_covariates(spatial_grid)
# Restrict to predictor columns and complete cases
spatial_grid <- spatial_grid[, intersect(names(spatial_grid), predictor_vars), drop = FALSE]
spatial_grid <- spatial_grid[complete.cases(spatial_grid), ]
cat("Spatial grid: ", nrow(spatial_grid), " cells with complete covariates\n", sep = "")

# Add 'unspecified' seagrass_species to all points on spatial grid
spatial_grid$seagrass_species <- "Unspecified"

# Use the trained model (fit_train) to predict on the full spatial grid
gpr_pred_full <- predict_gpr(fit_train, spatial_grid[, predictor_vars, drop = FALSE], se = TRUE)

gpr_full <- list(
  model = fit_train$model,
  predictions = gpr_pred_full$mean,
  se = gpr_pred_full$se,
  prediction_grid = spatial_grid,
  n_train = nrow(train_data),
  n_pred = nrow(spatial_grid)
)

# Attach metadata for deployment
gpr_full$value_var <- target_var
gpr_full$predictor_vars <- predictor_vars
gpr_full$kernel <- gpr_kernel

saveRDS(gpr_full, file.path(out_dir, "gpr_final_model.rds"))
cat("Saved", file.path(out_dir, "gpr_final_model.rds"), "\n")

grid_out <- gpr_full$prediction_grid
if (!"gpr_mean" %in% names(grid_out)) {
  grid_out$gpr_mean <- NA_real_
  # Clip negative values in gpr_mean to 0
  grid_out$gpr_se <- NA_real_
  pred_ok <- complete.cases(grid_out[, predictor_vars, drop = FALSE])
  grid_out$gpr_mean[pred_ok] <- gpr_full$predictions
  grid_out$gpr_mean <- pmax(grid_out$gpr_mean, 0)
  grid_out$gpr_se[pred_ok] <- gpr_full$se
}
saveRDS(grid_out, file.path(out_dir, "grid_predictions.rds"))
cat("Saved", file.path(out_dir, "grid_predictions.rds"), "\n\n")

# Plot a histogram of the predictions. N.B. this includes points on land (result of GPR not handling NA values properly). This could be trimmed via a land mask.
ggplot(grid_out, aes(x = gpr_se)) +
  geom_histogram(binwidth = max(grid_out$gpr_se, na.rm = TRUE) / 100) +
  labs(title = "Histogram of GPR predictions")
ggsave(file.path(figures_dir, "grid_predictions_histogram.png"), width = 8, height = 4, dpi = 150)
cat("Saved grid_predictions_histogram.png\n")

# -----------------------------------------------------------------------------
# Spatial prediction maps (shared plotting)
# -----------------------------------------------------------------------------

if (save_grid_plots && nrow(grid_out) > 0L && all(c("longitude", "latitude", "gpr_mean", "gpr_se") %in% names(grid_out))) {
  grid_plot <- subset(grid_out, !is.na(gpr_mean) & !is.na(longitude) & !is.na(latitude))
  obs_for_map <- dat[, c("longitude", "latitude", target_var)]
  obs_for_map <- obs_for_map[complete.cases(obs_for_map), ]
  spatial_plots <- plot_gpr_spatial_maps(
    grid = grid_plot,
    lon_range = lon_range,
    lat_range = lat_range,
    obs_data = obs_for_map,
    value_col = target_var,
    mean_title = "GPR Mean Predictions",
    se_title = "GPR Prediction Uncertainty",
    width = 8,
    height = 6,
    dpi = 150
  )
  if (!is.null(spatial_plots$p_mean)) {
    ggsave(file.path(figures_dir, "grid_plot_mean.png"), spatial_plots$p_mean, width = 8, height = 6, dpi = 150)
    cat("Saved grid_plot_mean.png\n")
  }
  if (!is.null(spatial_plots$p_se)) {
    ggsave(file.path(figures_dir, "grid_plot_se.png"), spatial_plots$p_se, width = 8, height = 6, dpi = 150)
    cat("Saved grid_plot_se.png\n")
  }
}

cat("\nPipeline complete.\nUse readRDS('gpr_final_model/gpr_final_model.rds') and predict_gpr() for deployment.\n")
