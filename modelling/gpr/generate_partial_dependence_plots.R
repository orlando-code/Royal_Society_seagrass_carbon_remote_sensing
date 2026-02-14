# Generate GPR partial dependence plots (PDPs)
#
# Standalone script to create PDPs from a fitted GPR model. Uses gpr_result from
# the global environment (e.g. after sourcing gpr_predictions.R), or loads from
# cache (figures/cv_pipeline_output/gpr_fitted_model.rds) if available.
#
# Run from project root: source("modelling/gpr/generate_partial_dependence_plots.R")
# Also called by run_paper_figures.R after Step 7.

if (!requireNamespace("here", quietly = TRUE)) install.packages("here", quiet = TRUE)
setwd(here::here())
set.seed(42)

source("modelling/R/helpers.R")
source("modelling/R/extract_covariates_from_rasters.R")
source("modelling/R/gpr_funs.R")
source("modelling/R/assign_region_from_latlon.R")
source("modelling/R/plot_config.R")
load_packages(c("here", "tidyverse", "ggplot2", "GauPro"))

target_var <- "median_carbon_density_100cm"
out_dir <- "figures/interpolation"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# PDP options (match gpr_predictions.R)
pd_n_points <- 100
pd_reference_method <- "conditional"
pd_n_samples <- 1000
pd_use_conditional <- TRUE
pd_range_quantiles <- c(0.01, 0.99)
pd_show_training_range <- TRUE
max_vars_to_plot <- 12

cat("\n========================================\n")
cat("GENERATING PARTIAL DEPENDENCE PLOTS\n")
cat("========================================\n\n")

# Get gpr_result: from global env (after gpr_predictions) or from cache
gpr_result <- NULL
if (exists("gpr_result", envir = .GlobalEnv) && !is.null(get("gpr_result", envir = .GlobalEnv))) {
  gpr_result <- get("gpr_result", envir = .GlobalEnv)
  cat("Using gpr_result from environment (from gpr_predictions.R)\n")
}
if (is.null(gpr_result)) {
  model_cache_path <- "figures/cv_pipeline_output/gpr_fitted_model.rds"
  if (file.exists(model_cache_path)) {
    gpr_result <- tryCatch(readRDS(model_cache_path), error = function(e) NULL)
    if (!is.null(gpr_result) && "model" %in% names(gpr_result)) {
      cat("Loaded fitted GPR from cache:", model_cache_path, "\n")
    } else {
      gpr_result <- NULL
    }
  }
}
if (is.null(gpr_result)) {
  stop(
    "No fitted GPR model found. Run gpr_predictions.R first (or run_paper_figures.R Step 7).\n",
    "To cache the model for future PDP runs, set use_cache_fitted_model <- TRUE in gpr_predictions.R before running."
  )
}

predictor_vars <- gpr_result$predictor_vars

# Load training data (needed for PDP reference values)
raster_covariates_dat <- readr::read_rds("data/all_extracted_new.rds")
raster_covariates_dat <- process_rs_covariates(raster_covariates_dat)
train_data <- raster_covariates_dat[, c(target_var, predictor_vars), drop = FALSE]
train_data <- train_data[complete.cases(train_data), ]

# Determine which variables to plot (exclude categoricals and spatial coords)
numeric_predictors <- predictor_vars[
  !predictor_vars %in% c("seagrass_species", "Region") &
  !predictor_vars %in% c("longitude", "latitude")
]
vars_to_plot <- head(numeric_predictors, max_vars_to_plot)

if (length(vars_to_plot) == 0) {
  warning("No numeric predictors available for partial dependence plots.")
} else {
  cat("Plotting partial dependence for", length(vars_to_plot), "variables\n")

  if (length(vars_to_plot) == 1) {
    p_pd <- plot_partial_dependence(
      gpr_model = gpr_result$model, dat = train_data, var_names = vars_to_plot,
      predictor_vars = predictor_vars, n_points = pd_n_points,
      reference_method = pd_reference_method, n_samples = pd_n_samples,
      use_conditional = pd_use_conditional, target_var_name = target_var,
      range_quantiles = pd_range_quantiles, show_training_range = pd_show_training_range,
      scale_params = gpr_result$scale_params,
      encoding = gpr_result$encoding, encoded_names = gpr_result$encoded_names
    )
    filename <- paste0("gpr_partial_dependence_", gsub("[^A-Za-z0-9]", "_", vars_to_plot[1]), ".png")
    ggplot2::ggsave(file.path(out_dir, filename), p_pd, width = 8, height = 6)
    cat("  Saved:", filename, "\n")
  } else {
    p_pd <- plot_partial_dependence(
      gpr_model = gpr_result$model, dat = train_data, var_names = vars_to_plot,
      predictor_vars = predictor_vars, n_points = pd_n_points,
      reference_method = pd_reference_method, n_samples = pd_n_samples,
      use_conditional = pd_use_conditional, target_var_name = target_var,
      range_quantiles = pd_range_quantiles, show_training_range = pd_show_training_range,
      scale_params = gpr_result$scale_params,
      encoding = gpr_result$encoding, encoded_names = gpr_result$encoded_names
    )
    ggplot2::ggsave(file.path(out_dir, "gpr_partial_dependence_all.png"), p_pd,
      width = 14, height = ceiling(length(vars_to_plot) / 2) * 4
    )
    cat("  Saved: gpr_partial_dependence_all.png\n")
    # also save each plot individually
    for (v in vars_to_plot) {
      p_pd_ind <- plot_partial_dependence(
        gpr_model = gpr_result$model, dat = train_data, var_names = v,
        predictor_vars = predictor_vars, n_points = pd_n_points,
        reference_method = pd_reference_method, n_samples = pd_n_samples,
        use_conditional = pd_use_conditional, target_var_name = target_var,
        range_quantiles = pd_range_quantiles, show_training_range = pd_show_training_range,
        scale_params = gpr_result$scale_params,
        encoding = gpr_result$encoding, encoded_names = gpr_result$encoded_names
      )
      filename <- paste0("gpr_partial_dependence_", gsub("[^A-Za-z0-9]", "_", v), ".png")
      ggplot2::ggsave(file.path(out_dir, filename), p_pd_ind, width = 8, height = 6)
    }
    cat("  Saved", length(vars_to_plot), "individual plots\n")
  }
  cat("\nPartial dependence plots saved to", out_dir, "\n")
}
