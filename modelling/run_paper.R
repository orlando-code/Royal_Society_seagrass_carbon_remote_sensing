# =============================================================================
# Paper pipeline driver: exploration, variable choice, prediction, supplement
#
# Runs the full pipeline in order to regenerate all paper and supplement figures.
# Uses cached data where possible (spatial folds, prediction grids) in output/cache.
# All outputs go under output/: cache, covariate_selection, cv_pipeline,
# final_models, predictions, supplement.
#
# Usage: setwd(project_root); source("modelling/run_paper.R")
#        Or: Rscript modelling/run_paper.R
# =============================================================================

setwd(here::here())
set.seed(42)  # globally

cat("\n")
cat(paste(rep("=", 94), collapse = ""), "\n")
cat("PAPER PIPELINE: CONFIGURE -> PRUNING -> TUNING -> CV -> IMPORTANCE -> PREDICTION -> SUPPLEMENT\n")
cat(paste(rep("=", 94), collapse = ""), "\n\n")

# -----------------------------------------------------------------------------
# -1. Configuring pipeline
# -----------------------------------------------------------------------------
cat("\t\tStep -1: Configuring paper pipeline ...\n")

# Plotting configuration
dpi <- 150
show_titles <- TRUE
assign("dpi", dpi, envir = .GlobalEnv)
assign("show_titles", show_titles, envir = .GlobalEnv)

target_var <- "median_carbon_density_100cm"
log_transform_target <- TRUE
# Region exclusion: set to character(0) to include all regions, or e.g. c("Black Sea") to exclude
# Black Sea has anomalously high carbon density that leaks to eastern Mediterranean predictions
# exclude_regions <- character(0)
exclude_regions <- c("Black Sea")

# Pruning configuration (run before CV pipeline)
use_correlation_filter       <- TRUE
correlation_filter_threshold <- 0.8
permutation_max_vars         <- 15L   # max vars retained after permutation importance
n_permutations               <- 1L # increase for paper run
permutation_coverage         <- 0.99  # cumulative importance coverage threshold
use_shap_per_model           <- TRUE  # per-model covariate sets via iml SHAP
do_cv_on_defaults <- TRUE  # whether to run Step 2 (default CV) + Step 2b

model_list    <- c("GPR", "GAM", "XGB")
n_folds       <- 5L

# --- CV strategy for tuning, pruning, importance, and final-model fallback ---
# cv_type controls fold construction across the whole pipeline (except the
# diagnostic CV scan in cv_pipeline.R, which always runs all methods).
#
#   "random"           – naive i.i.d. row-level split (optimistic; leaks duplicates)
#   "location_grouped" – group by identical (lon, lat); prevents coordinate leakage
#   "pixel_grouped"    – group by identical covariate vector; prevents both
#                        coordinate AND coarse-raster-pixel leakage (recommended)
#   "spatial"          – spatial blocks via blockCV at cv_blocksize metres
cv_type       <- get0("cv_type", envir = .GlobalEnv, ifnotfound = "pixel_grouped")
cv_blocksize  <- get0("cv_blocksize", envir = .GlobalEnv, ifnotfound = 1000L)  # metres; only used when cv_type = "spatial"

# cv_blocksize_scan: block sizes evaluated in the diagnostic CV scan
# (cv_pipeline.R always runs random, location-grouped, and pixel-grouped
# in addition to these spatial block sizes)
cv_blocksize_scan <- c(1000L, 5000L, 10000L, 100000L)

# --- Output directory structure ---
# Regime-specific outputs go under output/<cv_regime>/ so results from
# different splitting strategies are kept separate. Shared resources
# (cache, supplement) stay at the top level.
cv_regime_name <- switch(cv_type,
  random           = "random",
  location_grouped = "location_grouped",
  pixel_grouped    = "pixel_grouped",
  spatial          = paste0("spatial_", cv_blocksize, "m"),
  cv_type  # fallback
)
cv_output_dir <- file.path("output", cv_regime_name)

# Assign all globals for sourced sub-scripts
for (nm in c("target_var", "log_transform_target",
             "exclude_regions", "model_list",
             "use_correlation_filter", "correlation_filter_threshold",
             "permutation_max_vars", "permutation_coverage", "n_permutations",
             "use_shap_per_model", "n_folds", "cv_type", "cv_blocksize", "cv_blocksize_scan",
             "do_cv_on_defaults", "cv_output_dir")) {
  assign(nm, get(nm), envir = .GlobalEnv)
}

# Create directories: shared (cache, supplement) + regime-specific subdirs
for (d in c("output", "output/cache", "output/supplement",
            file.path(cv_output_dir, c("covariate_selection", "cv_pipeline",
                                        "final_models", "predictions")))) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# Flag to distinguish pre- vs post-tuning CV in cv_pipeline.R (to distinguish the outputs of steps 2 and 4)
post_tuning_validation <- FALSE
assign("post_tuning_validation", post_tuning_validation, envir = .GlobalEnv)

source("modelling/R/extract_covariates_from_rasters.R")
if (exists("raster_covariates"))
  assign("raster_covariates", raster_covariates, envir = .GlobalEnv)

cat("\n")
cat("  target_var:                  ", target_var, "\n")
cat("  log_transform_target:        ", log_transform_target, "\n")
cat("  exclude_regions:             ",
    if (length(exclude_regions) == 0) "none" else paste(exclude_regions, collapse = ", "), "\n")
cat("  use_correlation_filter:      ", use_correlation_filter, "\n")
cat("  correlation_filter_threshold:", correlation_filter_threshold, "\n")
cat("  permutation_max_vars:        ", permutation_max_vars, "\n")
cat("  n_permutations:              ", n_permutations, "\n")
cat("  use_shap_per_model:          ", use_shap_per_model, "\n")
cat("  model_list:                  ", paste(model_list, collapse = ", "), "\n")
cat("  n_folds:                     ", n_folds, "\n")
cat("  do_cv_on_defaults:          ", do_cv_on_defaults, "\n")
cat("  cv_type:                     ", cv_type, "\n")
cat("  cv_blocksize:               ", cv_blocksize, "\n")
cat("  cv_blocksize_scan:          ", if (length(cv_blocksize_scan)) paste(cv_blocksize_scan, collapse = ", ") else "none", "\n")
cat("  cv_output_dir:              ", cv_output_dir, "\n")
cat("  dpi:                         ", dpi, "\n")
cat("\n")


# -----------------------------------------------------------------------------
# 0. Data (build if missing)
# -----------------------------------------------------------------------------
if (!file.exists("data/all_extracted_new.rds")) {
  cat("\t\tStep 0: Building all_extracted_new.rds ...\n")
  source("modelling/pipeline/build_all_extracted_new.R")
  cat("\n")
} else {
  cat("\t\tStep 0: Using existing data/all_extracted_new.rds\n\n")
}

# -----------------------------------------------------------------------------
# 1. Covariate pruning (correlation and/or nested selection)
# -----------------------------------------------------------------------------
cat("\t\tStep 1: Covariate pruning (correlation and/or nested selection)\n")
source("modelling/pipeline/covariate_pruning_pipeline.R")
cat("\n")

# -----------------------------------------------------------------------------
# 2. CV pipeline (folds, CV results, inline plots)
# -----------------------------------------------------------------------------
if (isTRUE(get0("do_cv_on_defaults", envir = .GlobalEnv, ifnotfound = TRUE))) {
  cat("\t\tStep 2: CV pipeline (spatial folds, comparison, results)\n")
  source("modelling/pipeline/cv_pipeline.R")
  cat("\n")

  # -----------------------------------------------------------------------------
  # 2b. Spatial / lat-lon / region effect (all models)
  # -----------------------------------------------------------------------------
  cat("\t\tStep 2b: Spatial and categorical effect (lat, lon, region) for all models (default hyperparams)\n")
  source("modelling/plots/spatial_categorical_effect_all_models.R")
  cat("\n")
} else {
  cat("\t\tStep 2: Skipping default CV pipeline (do_cv_on_defaults=FALSE)\n\n")
}

# -----------------------------------------------------------------------------
# 3. Hyperparameter tuning (all models; cv_type; progress)
# -----------------------------------------------------------------------------
cat("\t\tStep 3: Hyperparameter tuning (GPR, GAM, XGB)\n")
source("modelling/pipeline/hyperparameter_tuning_pipeline.R")
cat("\n")

# -----------------------------------------------------------------------------
# 3b. Permutation and SHAP importance (best config + best vars per model)
# -----------------------------------------------------------------------------
cat("\t\tStep 3b: Permutation importance (tuned models)\n")
source("modelling/pipeline/permutation_importance_final.R")
cat("\n")
cat("\t\tStep 3c: SHAP importance (tuned models)\n")
source("modelling/pipeline/shap_importance_final.R")
cat("\n")

# -----------------------------------------------------------------------------
# 4. Cross-fold validation (tuned models and pruned covariates)
# -----------------------------------------------------------------------------
post_tuning_validation <- TRUE  # now using tuned models
assign("post_tuning_validation", post_tuning_validation, envir = .GlobalEnv)
cat("\t\tStep 4: Cross-fold validation (tuned models)\n")
source("modelling/pipeline/cv_pipeline.R")
cat("\n")

# -----------------------------------------------------------------------------
# 4b. Spatial / lat-lon / region effect (tuned models)
# -----------------------------------------------------------------------------
categorical_use_tuned <- TRUE # use tuned models for spatial and categorical effect plots
assign("categorical_use_tuned", categorical_use_tuned, envir = .GlobalEnv)
cat("\t\tStep 4b: Spatial and categorical effect (lat, lon, region) for all models (tuned hyperparams)\n")
source("modelling/plots/spatial_categorical_effect_all_models.R")
cat("\n")

# -----------------------------------------------------------------------------
# 5. Fit and save final models (all data); XGB, GAM, GPR
# -----------------------------------------------------------------------------
cat("\t\tStep 5: Fit and save final models\n")
source("modelling/pipeline/fit_final_models.R")
cat("\n")

# -----------------------------------------------------------------------------
# 6. Partial dependence plots for all final models
# -----------------------------------------------------------------------------
cat("\t\tStep 6: Partial dependence plots (XGB, GAM, GPR)\n")
source("modelling/plots/partial_dependence_all_models.R")
cat("\n")


# -----------------------------------------------------------------------------
# 7. Spatial prediction maps for each model
# -----------------------------------------------------------------------------
cat("\t\tStep 7: Spatial prediction maps (XGB, GAM, GPR)\n")
source("modelling/plots/spatial_prediction_maps.R")
cat("\n")

# -----------------------------------------------------------------------------
# 8. Supplement (shared, not regime-specific)
# -----------------------------------------------------------------------------
cat("\t\tStep 8: Supplement figures\n")
source("modelling/plots/supplement.R")
cat("\n")

# -----------------------------------------------------------------------------
# Done
# -----------------------------------------------------------------------------
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("PAPER PIPELINE COMPLETE\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")
cat("Outputs (cv_regime = ", cv_regime_name, "):\n", sep = "")
cat("  output/cache/                           – Shared cached spatial folds and prediction grids\n")
cat("  output/supplement/                      – Shared supplement figures\n")
cat("  ", cv_output_dir, "/covariate_selection/ – Covariate pruning results\n", sep = "")
cat("  ", cv_output_dir, "/cv_pipeline/         – CV results, tuning, importance plots\n", sep = "")
cat("  ", cv_output_dir, "/final_models/        – Fitted model RDS (XGB, GAM, GPR)\n", sep = "")
cat("  ", cv_output_dir, "/predictions/         – Prediction maps, PDPs\n", sep = "")
cat("\n")
