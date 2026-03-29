# =============================================================================
# Robust (across-seed) pixel-grouped pipeline driver
#
# Seed-robust analogue of `modelling/run_paper.R`, focused on the
# pixel_grouped CV strategy with multiple seed instantiations. Evaluates model performance across
# multiple fold-seed instantiations to characterise split-variance in a
# small dataset where a single random split can be unrepresentative.
#
# Pipeline phases:
#   Phase I  – Warm-start (optional): baseline pruning + tuning
#   Phase II – Robust selection: multi-seed tuning, SHAP pruning, re-tuning
#   Phase III – Evaluation & diagnostics
#
# Multi-seed scripts live in `modelling/multiseed/`.
#
# Seed policy:
#   - Global set.seed(42) for reproducible non-seed-controlled randomness.
#   - robust_fold_seed_list (from pipeline_config): seeds for selection/tuning phase.
#     Each produces an independent fold instantiation.
#   - eval_fold_seed_list (from pipeline_config): seeds for final evaluation.
#     Separate from selection seeds to prevent optimistic bias
#     ("optimise on A, evaluate on B").
#   - Grid sampling in robust tuning uses set.seed(42).
#   - SHAP: set.seed(seed + fold_k) per (seed, fold) for determinism.
#
# Usage: setwd(project_root); source("<robust pipeline driver>")
#        Or: Rscript <robust pipeline driver>
# =============================================================================

setwd(here::here())
set.seed(42)

cat("\n")
cat(paste(rep("=", 94), collapse = ""), "\n")
cat("ROBUST PIXEL_GROUPED: CONFIG -> ROBUST-PRUNE -> ROBUST-TUNE -> ROBUST-EVAL -> FINAL-MODELS\n")
cat(paste(rep("=", 94), collapse = ""), "\n\n")

# -----------------------------------------------------------------------------
# -1. Configuration
# -----------------------------------------------------------------------------
cat("\t\tStep -1: Configuring robust pipeline ...\n")

source("modelling/config/pipeline_config.R")
cfg <- get_pipeline_config()

# Plotting
dpi <- cfg$dpi
show_titles <- cfg$show_titles

# Response
target_var           <- cfg$target_var
log_transform_target <- isTRUE(cfg$log_transform_target)

# Region exclusion (Black Sea has anomalously high carbon density)
exclude_regions <- cfg$exclude_regions

# Models
model_list <- cfg$model_list
include_seagrass_species <- isTRUE(cfg$include_seagrass_species)

# Baseline covariate pruning (used in warm-start and as fallback)
use_correlation_filter       <- isTRUE(cfg$use_correlation_filter)
correlation_filter_threshold <- cfg$correlation_filter_threshold
permutation_max_vars         <- as.integer(cfg$permutation_max_vars)
n_permutations               <- as.integer(cfg$n_permutations)
permutation_coverage         <- as.numeric(cfg$permutation_coverage)
use_shap_per_model           <- isTRUE(cfg$use_shap_per_model)
shap_selection_policy        <- cfg$shap_selection_policy

# CV configuration
n_folds      <- as.integer(cfg$n_folds)
cv_type      <- cfg$cv_type
cv_blocksize <- cfg$cv_blocksize  # unused for pixel_grouped; kept for interface consistency

# Output directory structure
cv_regime_name <- cfg$cv_regime_name
cv_output_dir  <- cfg$cv_output_dir
cv_type_label  <- cfg$cv_type_label

# Seed policy
robust_fold_seed_list <- as.integer(cfg$robust_fold_seed_list)   # seeds for selection/tuning
eval_fold_seed_list   <- as.integer(cfg$eval_fold_seed_list)     # seeds for held-out evaluation
run_output_dir <- file.path(
  cv_output_dir,
  "cv_pipeline",
  build_seeded_run_folder_name(
    cv_type_label = cv_type_label,
    folder_type = "evaluation",
    repeat_seed_list = eval_fold_seed_list,
    robust_seed_list = robust_fold_seed_list,
    include_seed_values = TRUE
  )
)

# SHAP pruning controls
shap_n_points       <- as.integer(cfg$shap_n_points)
shap_folds_per_seed <- as.integer(cfg$shap_folds_per_seed)
shap_max_gpr_train  <- as.integer(cfg$shap_max_gpr_train)
shap_selection_coverage_grid <- as.numeric(cfg$shap_selection_coverage_grid)
shap_selection_max_vars_grid <- as.integer(cfg$shap_selection_max_vars_grid)

# Robust pruning importance type
robust_pruned_importance_type <- cfg$robust_pruned_importance_type
robust_rmse_lambda <- as.numeric(cfg$robust_rmse_lambda)

# Pipeline toggles (control runtime for development vs publication)
do_warm_start          <- isTRUE(cfg$do_warm_start)
do_shap_refined_tuning <- isTRUE(cfg$do_shap_refined_tuning)
do_sensitivity         <- isTRUE(cfg$do_sensitivity)
# Optional nested robust seed-count sweep (very expensive):
# re-run robust tuning/pruning/eval at multiple robust seed counts.
# This is now executed via standalone script:
#   Rscript modelling/analysis/tuning_seed_sweep.R
do_tuning_seed_sweep <- isTRUE(cfg$do_tuning_seed_sweep)
tuning_seed_sweep_counts <- as.integer(cfg$tuning_seed_sweep_counts)
tuning_seed_pool <- as.integer(cfg$tuning_seed_pool)
tuning_sweep_eval_seed_list <- as.integer(cfg$tuning_sweep_eval_seed_list)
tuning_seed_sampling <- cfg$tuning_seed_sampling
do_tuning_seed_sweep_refined_tuning <- isTRUE(cfg$do_tuning_seed_sweep_refined_tuning)
do_diagnostics         <- isTRUE(cfg$do_diagnostics)
do_fit_final_models    <- isTRUE(cfg$do_fit_final_models)
do_supplement          <- isTRUE(cfg$do_supplement)
xgb_max_grid_candidates <- as.integer(cfg$xgb_max_grid_candidates)
n_lons <- as.integer(cfg$n_lons)
n_lats <- as.integer(cfg$n_lats)

# Assign all config to .GlobalEnv for sourced sub-scripts
config_vars <- c(
  "dpi", "show_titles",
  "target_var", "log_transform_target", "exclude_regions",
  "model_list",
  "include_seagrass_species",
  "use_correlation_filter", "correlation_filter_threshold",
  "permutation_max_vars", "n_permutations", "permutation_coverage",
  "use_shap_per_model", "shap_selection_policy",
  "n_folds", "cv_type", "cv_blocksize",
  "cv_regime_name", "cv_output_dir", "cv_type_label",
  "run_output_dir",
  "robust_fold_seed_list", "eval_fold_seed_list",
  "shap_n_points", "shap_folds_per_seed", "shap_max_gpr_train",
  "shap_selection_coverage_grid", "shap_selection_max_vars_grid",
  "robust_rmse_lambda",
  "do_tuning_seed_sweep", "tuning_seed_sweep_counts", "tuning_seed_pool",
  "tuning_sweep_eval_seed_list", "tuning_seed_sampling",
  "do_tuning_seed_sweep_refined_tuning",
  "robust_pruned_importance_type",
  "xgb_max_grid_candidates",
  "n_lons", "n_lats"
)
for (nm in config_vars) assign(nm, get(nm), envir = .GlobalEnv)

# Create output directories
for (d in c("output", "output/cache",
            file.path(cv_output_dir, c("covariate_selection", "cv_pipeline")))) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# Persist effective run configuration for reproducibility/auditing.
run_metadata_dir <- file.path(cv_output_dir, "run_metadata")
dir.create(run_metadata_dir, recursive = TRUE, showWarnings = FALSE)
run_cfg_vars <- unique(c(
  config_vars,
  "do_warm_start", "do_shap_refined_tuning", "do_sensitivity",
  "do_diagnostics", "do_fit_final_models", "do_supplement",
  "cv_output_dir"
))
effective_run_cfg <- setNames(lapply(run_cfg_vars, function(nm) get(nm, envir = .GlobalEnv)), run_cfg_vars)
saveRDS(effective_run_cfg, file.path(run_metadata_dir, "pipeline_config_effective.rds"))
dput(effective_run_cfg, file = file.path(run_metadata_dir, "pipeline_config_effective.dput"))
cat("  wrote run config to:         ", file.path(run_metadata_dir, "pipeline_config_effective.rds"), "\n")

source("modelling/R/extract_covariates_from_rasters.R")
if (exists("raster_covariates"))
  assign("raster_covariates", raster_covariates, envir = .GlobalEnv)

cat("\n")
cat("  target_var:                  ", target_var, "\n")
cat("  log_transform_target:        ", log_transform_target, "\n")
cat("  exclude_regions:             ",
    if (length(exclude_regions) == 0) "none" else paste(exclude_regions, collapse = ", "), "\n")
cat("  model_list:                  ", paste(model_list, collapse = ", "), "\n")
cat("  include_seagrass_species:    ", include_seagrass_species, "\n")
cat("  n_folds:                     ", n_folds, "\n")
cat("  cv_type (baseline):          ", cv_type, "\n")
cat("  cv_output_dir:               ", cv_output_dir, "\n")
cat("  run_output_dir:              ", run_output_dir, "\n")
cat("  robust_fold_seed_list:       ", paste(robust_fold_seed_list, collapse = ", "), "\n")
cat("  eval_fold_seed_list:         ", paste(eval_fold_seed_list, collapse = ", "), "\n")
cat("  shap_n_points:               ", shap_n_points, "\n")
cat("  shap_folds_per_seed:         ", shap_folds_per_seed, "\n")
cat("  shap_max_gpr_train:          ", shap_max_gpr_train, "\n")
cat("  shap_selection_coverage_grid:", paste(shap_selection_coverage_grid, collapse = ", "), "\n")
cat("  shap_selection_max_vars_grid:", paste(shap_selection_max_vars_grid, collapse = ", "), "\n")
cat("  robust_pruned_importance:    ", robust_pruned_importance_type, "\n")
cat("  robust_rmse_lambda:          ", robust_rmse_lambda, "\n")
cat("  do_warm_start:               ", do_warm_start, "\n")
cat("  do_shap_refined_tuning:      ", do_shap_refined_tuning, "\n")
cat("  do_sensitivity:              ", do_sensitivity, "\n")
cat("  do_tuning_seed_sweep:        ", do_tuning_seed_sweep, "\n")
cat("  tuning_seed_sweep_counts:    ", paste(tuning_seed_sweep_counts, collapse = ", "), "\n")
cat("  tuning_seed_pool:            ", paste(tuning_seed_pool, collapse = ", "), "\n")
cat("  tuning_sweep_eval_seed_list: ", paste(tuning_sweep_eval_seed_list, collapse = ", "), "\n")
cat("  tuning_seed_sampling:        ", tuning_seed_sampling, "\n")
cat("  sweep_refined_tuning:        ", do_tuning_seed_sweep_refined_tuning, "\n")
cat("  do_diagnostics:              ", do_diagnostics, "\n")
cat("  do_fit_final_models:         ", do_fit_final_models, "\n")
cat("  do_supplement:               ", do_supplement, "\n")
cat("  dpi:                         ", dpi, "\n")
cat("\n")

# =============================================================================
# Step 0: Data (build if missing)
# =============================================================================
cat("\t\tStep 0: Data (build if missing)\n")
if (!file.exists("data/all_extracted_new.rds")) {
  source("modelling/pipeline/build_all_extracted_new.R")
} else {
  cat("  Using existing data/all_extracted_new.rds\n")
}
cat("\n")

# =============================================================================
# Phase I: Warm-start (baseline pruning + tuning on deterministic folds)
#
# Optional speed/initialization phase. Robust scripts now run without this:
#   - hyperparameters fall back to model defaults when no saved configs exist
#   - covariates fall back to correlation-filtered predictors when no pruned
#     variable file exists
# =============================================================================
if (isTRUE(do_warm_start)) {
  # -------------------------------------------------------------------------
  # Step 1: Baseline covariate pruning (cv_type = pixel_grouped)
  # -------------------------------------------------------------------------
  cat("\t\tStep 1: Baseline covariate pruning (cv_type=pixel_grouped)\n")
  assign("cv_type", "pixel_grouped", envir = .GlobalEnv)
  source("modelling/pipeline/covariate_pruning_pipeline.R")
  cat("\n")

  # -------------------------------------------------------------------------
  # Step 2: Baseline hyperparameter tuning -> best_config_*.rds
  # -------------------------------------------------------------------------
  cat("\t\tStep 2: Baseline hyperparameter tuning (cv_type=pixel_grouped)\n")
  source("modelling/pipeline/hyperparameter_tuning_pipeline.R")
  cat("\n")
} else {
  cat("\t\tPhase I: Skipping warm-start (do_warm_start=FALSE)\n")
  cat("  Robust scripts will use default hyperparameters + correlation covariates where needed.\n\n")
}

# =============================================================================
# Phase II: Robust selection (multi-seed tuning + SHAP pruning)
# cv_type remains "pixel_grouped"; robust behavior comes from varying seed lists.
# =============================================================================
cat("  --- Phase II: seeded pixel_grouped robust selection ---\n\n")
assign("cv_type", "pixel_grouped", envir = .GlobalEnv)

# -----------------------------------------------------------------------------
# Step 3: Robust hyperparameter tuning (initial covariates)
#
# Uses robust_fold_seed_list to evaluate each candidate config across multiple
# fold instantiations. Covariate sets come from the fallback chain:
#   1. Robust SHAP-pruned vars (if available from a previous Step 4 run)
#   2. Baseline SHAP/perm vars from Step 1
# -----------------------------------------------------------------------------
cat("\t\tStep 3: Robust hyperparameter tuning (initial covariates)\n")
assign("robust_pruned_importance_type", robust_pruned_importance_type, envir = .GlobalEnv)
source("modelling/multiseed/robust_hyperparameter_tuning.R")
cat("\n")

# -----------------------------------------------------------------------------
# Step 4: Robust SHAP covariate pruning
#
# Computes SHAP importance on fold-specific training sets across multiple
# seeds, then selects top covariates per model based on averaged |SHAP|.
# Uses robust hyperparams from Step 3 (or baseline fallback).
# -----------------------------------------------------------------------------
cat("\t\tStep 4: Robust SHAP covariate pruning\n")
source("modelling/multiseed/robust_shap_covariate_pruning.R")
cat("\t\tStep 4b: Robust SHAP importance plots (mean +/- SD)\n")
source("modelling/plots/robust_shap_importance_plots.R")
cat("\n")

# -----------------------------------------------------------------------------
# Step 5: Robust hyperparameter tuning with SHAP-selected covariates
#
# Re-runs robust tuning using the SHAP-pruned covariate sets from Step 4.
# This iterative refinement (tune -> prune -> re-tune) is the core
# methodology for the small-dataset setting.
# -----------------------------------------------------------------------------
if (isTRUE(do_shap_refined_tuning)) {
  cat("\t\tStep 5: Robust hyperparameter tuning (SHAP-selected covariates)\n")
  source("modelling/multiseed/robust_hyperparameter_tuning.R")
  cat("\n")
} else {
  cat("\t\tStep 5: Skipping SHAP-refined tuning (do_shap_refined_tuning=FALSE)\n\n")
}

# =============================================================================
# Phase III: Evaluation & diagnostics
# =============================================================================

# -----------------------------------------------------------------------------
# Step 6: Robust evaluation across held-out seeds
#
# Evaluates the final (tuned covariates + tuned hyperparams) on
# eval_fold_seed_list, which is disjoint from robust_fold_seed_list.
# -----------------------------------------------------------------------------
cat("\t\tStep 6: Robust evaluation across eval seeds\n")
assign("cv_type", "pixel_grouped", envir = .GlobalEnv)
source("modelling/multiseed/robust_evaluation.R")
cat("\n")

# -----------------------------------------------------------------------------
# Step 7: R2 sensitivity analysis (performance vs seed count, fold count)
# -----------------------------------------------------------------------------
if (isTRUE(do_sensitivity)) {
  cat("\t\tStep 7: R2 sensitivity analysis\n")
  source("modelling/analysis/sensitivity_suite.R")
  source("modelling/plots/sensitivity_plots.R")
  if (isTRUE(do_tuning_seed_sweep)) {
    cat("  NOTE: tuning seed-count sweep is standalone.\n")
    cat("  Run separately (tmux-friendly): Rscript modelling/analysis/tuning_seed_sweep.R\n")
  }
  cat("\n")
} else {
  cat("\t\tStep 7: Skipping R2 sensitivity analysis (do_sensitivity=FALSE)\n\n")
}

# -----------------------------------------------------------------------------
# Step 8: Train/test fraction & target-variance diagnostic
# -----------------------------------------------------------------------------
if (isTRUE(do_diagnostics)) {
  cat("\t\tStep 8: Train/test fraction & target-variance diagnostics\n")
  perf_detailed_override <- file.path(run_output_dir, "by_seed_detailed.csv")
  assign("perf_detailed_csv_override", perf_detailed_override, envir = .GlobalEnv)
  assign("fold_seed_list", eval_fold_seed_list, envir = .GlobalEnv)
  source("modelling/multiseed/train_test_fraction_diagnostic.R")
  cat("\n")
} else {
  cat("\t\tStep 8: Skipping diagnostics (do_diagnostics=FALSE)\n\n")
}

# -----------------------------------------------------------------------------
# Step 9: Fit and save final models using robust tuning/covariates
# -----------------------------------------------------------------------------
if (isTRUE(do_fit_final_models)) {
  cat("\t\tStep 9: Fit and save final models (robust configs/vars)\n")
  seeds_sel_str <- paste(robust_fold_seed_list, collapse = "-")
  assign("use_robust_final_configs", TRUE, envir = .GlobalEnv)
  assign("robust_pruned_csv_override",
         file.path(
           "output", cv_regime_name, "covariate_selection", "robust_pixel_grouped",
           paste0("pruned_model_variables_shap_robust_pixel_grouped_seeds_", seeds_sel_str, ".csv")
         ),
         envir = .GlobalEnv)
  source("modelling/pipeline/fit_final_models.R")
  cat("\n")
} else {
  cat("\t\tStep 9: Skipping final model fit/save (do_fit_final_models=FALSE)\n\n")
}

# -----------------------------------------------------------------------------
# Step 10: Partial dependence plots
# -----------------------------------------------------------------------------
if (isTRUE(do_fit_final_models)) {
  cat("\t\tStep 10: Partial dependence plots (from final models)\n")
  source("deprecated/modelling/plots/partial_dependence_all_models.R")
  cat("\n")
} else {
  cat("\t\tStep 10: Skipping PDPs (requires do_fit_final_models=TRUE)\n\n")
}

# -----------------------------------------------------------------------------
# Step 11: Supplement figures
# -----------------------------------------------------------------------------
if (isTRUE(do_supplement)) {
  cat("\t\tStep 11: Supplement figures\n")
  source("modelling/plots/supplement.R")
  cat("\n")
} else {
  cat("\t\tStep 11: Skipping supplement figures (do_supplement=FALSE)\n\n")
}

# =============================================================================
# Done
# =============================================================================
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("ROBUST PIXEL_GROUPED PIPELINE COMPLETE\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")
cat("Seed policy:\n")
cat("  Selection/tuning seeds: ", paste(robust_fold_seed_list, collapse = ", "), "\n")
cat("  Evaluation seeds:       ", paste(eval_fold_seed_list, collapse = ", "), "\n")
cat("  Separation:             optimise on selection seeds, evaluate on eval seeds\n\n")
cat("Outputs (cv_regime = ", cv_regime_name, "):\n", sep = "")
cat("  ", cv_output_dir, "/covariate_selection/ – Covariate pruning (baseline + robust SHAP)\n", sep = "")
cat("  ", run_output_dir, " – This run's evaluation outputs (incl. sensitivity/diagnostics)\n", sep = "")
cat("  ", cv_output_dir, "/cv_pipeline/         – Robust tuning artifacts + all run folders\n", sep = "")
cat("\n")
