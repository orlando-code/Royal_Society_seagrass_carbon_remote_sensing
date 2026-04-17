# =============================================================================
# Robust (across-seed) pixel-grouped pipeline driver
#
# Seed-robust analogue of `modelling/run_paper.R`, focused on the
# pixel_grouped CV strategy with multiple seed instantiations. Evaluates model performance across
# multiple fold-seed instantiations to characterise split-variance in a
# small dataset where a single random split can be unrepresentative.
#
# Pipeline phases:
#   Phase I  – Robust selection: multi-seed tuning, SHAP pruning, re-tuning
#   Phase II – Evaluation & diagnostics
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
if (!exists("seagrass_init_repo", mode = "function", inherits = TRUE)) source("modelling/R/init_repo.R")
project_root <- seagrass_init_repo(
  packages = c("here", "dplyr", "readr", "ggplot2", "iml", "sf"),
  source_files = c(
    "modelling/R/assign_region_from_latlon.R",
    "modelling/pipeline_config.R",
    "modelling/plots/plot_config.R"
  ),
  include_helpers = TRUE,
  require_core_inputs = TRUE,
  check_renv = TRUE
)
set.seed(42)

cat("\n")
cat(paste(rep("=", 94), collapse = ""), "\n")
cat("ROBUST PIXEL_GROUPED: CONFIG -> ROBUST-TUNE -> ROBUST-SHAP-PRUNE -> ROBUST-TUNE-SHAP-PRUNED -> ROBUST-EVAL -> FINAL-MODELS\n")
cat(paste(rep("=", 94), collapse = ""), "\n\n")

# -----------------------------------------------------------------------------
# -1. Configuration
# -----------------------------------------------------------------------------
cat("\t\tStep -1: Configuring robust pipeline ...\n")

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

# Baseline covariate pruning (fallback when no pruned-variable files exist)
use_correlation_filter       <- isTRUE(cfg$use_correlation_filter)
correlation_filter_threshold <- cfg$correlation_filter_threshold
permutation_max_vars         <- as.integer(cfg$permutation_max_vars)
n_permutations               <- as.integer(cfg$n_permutations)
feature_coverage         <- as.numeric(cfg$feature_coverage)
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

use_robust_seeds_from_tuning_sweep <- isTRUE(cfg$use_robust_seeds_from_tuning_sweep)
use_paper_seed_registry <- isTRUE(cfg$use_paper_seed_registry)
chosen_seeds_rds <- file.path(
  "output", "tuning_seed_sweep_runs", "chosen_seeds_latest.rds"
)
if (isTRUE(use_robust_seeds_from_tuning_sweep)) {
  if (!file.exists(chosen_seeds_rds)) {
    stop("use_robust_seeds_from_tuning_sweep is TRUE but file not found:\n  ", chosen_seeds_rds)
  }
  ch <- readRDS(chosen_seeds_rds)
  robust_fold_seed_list <- as.integer(ch$robust_fold_seed_list)
  if (is.null(robust_fold_seed_list) || length(robust_fold_seed_list) == 0L || any(!is.finite(robust_fold_seed_list))) {
    stop("chosen_seeds_latest.rds missing valid robust_fold_seed_list")
  }
  if (!is.null(ch$eval_fold_seed_list) && length(ch$eval_fold_seed_list) > 0L) {
    eval_fold_seed_list <- as.integer(ch$eval_fold_seed_list)
  }
  cat("  Loaded seeds from tuning sweep:\n    ", chosen_seeds_rds, "\n", sep = "")
}

if (isTRUE(use_paper_seed_registry)) {
  sr <- cfg$seed_registry
  if (is.null(sr$paper_robust_fold_seed_list) || length(sr$paper_robust_fold_seed_list) == 0L) {
    stop("use_paper_seed_registry is TRUE but cfg$seed_registry$paper_robust_fold_seed_list is missing/empty.")
  }
  robust_fold_seed_list <- as.integer(sr$paper_robust_fold_seed_list)
}

seed_subset <- paste(robust_fold_seed_list, collapse = "-")
multiseed_run_output_id <- cfg$multiseed_run_output_id
run_output_dir <- file.path("output", paste0("pixel_grouped_", seed_subset))
# All pipeline scripts that read `cv_output_dir` (final models, maps, tuning
# configs, baseline covariate dirs, etc.) must resolve under this run folder,
# not the generic cfg$cv_output_dir (output/pixel_grouped).
cv_output_dir <- run_output_dir

# assign run_output_dir to global so all subsequent files can find it
assign("run_output_dir", run_output_dir, envir = .GlobalEnv)

# SHAP pruning controls
shap_n_points       <- as.integer(cfg$shap_n_points)
shap_folds_per_seed <- as.integer(cfg$shap_folds_per_seed)
shap_max_gpr_train  <- as.integer(cfg$shap_max_gpr_train)
shap_selection_coverage_grid <- as.numeric(cfg$shap_selection_coverage_grid)
shap_selection_max_vars_grid <- as.integer(cfg$shap_selection_max_vars_grid)

# Robust pruning importance type
robust_pruned_importance_type <- cfg$robust_pruned_importance_type

# Pipeline toggles (control runtime for development vs publication)
do_shap_refined_tuning <- isTRUE(cfg$do_shap_refined_tuning)
do_sensitivity         <- isTRUE(cfg$do_sensitivity)
do_diagnostics         <- isTRUE(cfg$do_diagnostics)
do_fit_final_models    <- isTRUE(cfg$do_fit_final_models)
do_supplement          <- isTRUE(cfg$do_supplement)
xgb_max_grid_candidates <- as.integer(cfg$xgb_max_grid_candidates)
n_lons <- as.integer(cfg$n_lons)
n_lats <- as.integer(cfg$n_lats)


config_vars <- c(
    "dpi", "show_titles",
    "target_var", "log_transform_target", "exclude_regions",
    "model_list",
    "include_seagrass_species",
    "use_correlation_filter", "correlation_filter_threshold",
    "permutation_max_vars", "n_permutations", "feature_coverage",
    "use_shap_per_model", "shap_selection_policy",
    "n_folds", "cv_type", "cv_blocksize",
    "cv_regime_name", "cv_output_dir", "cv_type_label",
    "run_output_dir",
    "robust_fold_seed_list", "eval_fold_seed_list",
    "use_paper_seed_registry",
    "multiseed_run_output_id", "use_robust_seeds_from_tuning_sweep",
    "shap_n_points", "shap_folds_per_seed", "shap_max_gpr_train",
    "shap_selection_coverage_grid", "shap_selection_max_vars_grid",
    "robust_pruned_importance_type",
    "xgb_max_grid_candidates",
    "n_lons", "n_lats"
  )
apply_pipeline_defaults(
  cfg,
  config_vars,
  envir = .GlobalEnv
)

# Create output directories
for (d in c("output", "output/cache",
            file.path(cv_output_dir, c("covariate_selection", "cv_pipeline")))) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# Persist effective run configuration for reproducibility/auditing.
run_metadata_dir <- file.path(run_output_dir, "_run_metadata")
dir.create(run_metadata_dir, recursive = TRUE, showWarnings = FALSE)
run_cfg_vars <- unique(c(
  config_vars,
  "do_shap_refined_tuning", "do_sensitivity",
  "do_diagnostics", "do_fit_final_models", "do_supplement",
  "cv_output_dir"
))
effective_run_cfg <- setNames(lapply(run_cfg_vars, function(nm) get(nm, envir = .GlobalEnv)), run_cfg_vars)
seagrass_confirm_same_config(
  current_cfg = effective_run_cfg,
  search_root = "output",
  exclude_keys = c("run_output_dir", "multiseed_run_output_id"),
  prompt_prefix = "multiseed pipeline",
  path_regex = "output/pixel_grouped_[^/]+/_run_metadata/pipeline_config_effective\\.rds$"
)
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
cat("  cv_type:          ", cv_type, "\n")
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
cat("  use_robust_seeds_from_sweep: ", use_robust_seeds_from_tuning_sweep, "\n")
cat("  use_paper_seed_registry:     ", use_paper_seed_registry, "\n")
cat("  do_shap_refined_tuning:      ", do_shap_refined_tuning, "\n")
cat("  do_sensitivity:              ", do_sensitivity, "\n")
cat("  do_diagnostics:              ", do_diagnostics, "\n")
cat("  do_fit_final_models:         ", do_fit_final_models, "\n")
cat("  do_supplement:               ", do_supplement, "\n")
cat("  dpi:                         ", dpi, "\n")
cat("\n")

# =============================================================================
# Step 0: Data (build if missing)
# =============================================================================
cat("\t\tStep 0: Data (build if missing)\n")
if (!dir.exists("data/env_rasters")) {
  stop("Required input directory missing: data/env_rasters")
}
if (!file.exists("data/all_extracted_new.rds")) {
  source("modelling/R/build_all_extracted_new.R")
} else {
  cat("  Using existing data/all_extracted_new.rds\n")
}
if (!file.exists("data/all_extracted_new.rds")) {
  stop("Required input file missing after build step: data/all_extracted_new.rds")
}
cat("\n")

# =============================================================================
# Phase I: Robust selection (multi-seed tuning + SHAP pruning)
# =============================================================================
cat("  --- Phase I: seeded pixel_grouped robust selection ---\n\n")
assign("cv_type", "pixel_grouped", envir = .GlobalEnv)

#-----------------------------------------------------------------------------
cat("\t\tStep 1: Robust hyperparameter tuning (initial covariates)\n")
assign("robust_pruned_importance_type", robust_pruned_importance_type, envir = .GlobalEnv)
source("modelling/multiseed/robust_hyperparameter_tuning.R")
cat("\n")

# -----------------------------------------------------------------------------
# Step 2: Robust SHAP covariate pruning
#
# Computes SHAP importance on fold-specific training sets across multiple
# seeds, then selects top covariates per model based on averaged |SHAP|.
# -----------------------------------------------------------------------------
cat("\t\tStep 2: Robust SHAP covariate pruning\n")
source("modelling/multiseed/robust_shap_covariate_pruning.R")
cat("\t\tStep 2b: Robust SHAP importance plots\n")
source("modelling/plots/plot_robust_shap_importance.R")
cat("\n")

# cache robust hyperparameter tuning outputs
file.rename(file.path(run_output_dir, "cv_pipeline"), file.path(run_output_dir, "cv_pipeline_pre_pruning"))

# -----------------------------------------------------------------------------
# Step 3: Robust hyperparameter tuning with SHAP-selected covariates
#
# Re-runs robust tuning using the SHAP-pruned covariate sets from Step 2.
# -----------------------------------------------------------------------------
if (isTRUE(do_shap_refined_tuning)) {
  cat("\t\tStep 3: Robust hyperparameter tuning (using SHAP-pruned covariates)\n")
  source("modelling/multiseed/robust_hyperparameter_tuning.R")
  cat("\n")
} else {
  cat("\t\tStep 3: Skipping SHAP-refined tuning (do_shap_refined_tuning=FALSE)\n\n")
}
# cache robust shap pruning outputs
file.rename(file.path(run_output_dir, "covariate_selection"), file.path(run_output_dir, "covariate_selection_pre_pruning"))

# -----------------------------------------------------------------------------
# Step 3b: SHAP importance plots for SHAP-selected covariates
#
# Plots the SHAP importance of the SHAP-selected covariates.
# -----------------------------------------------------------------------------
cat("\t\tStep 3b: SHAP importance plots for SHAP-selected covariates\n")
source("modelling/multiseed/robust_shap_covariate_pruning.R")
source("modelling/plots/plot_robust_shap_importance.R")
cat("\n")
# =============================================================================
# Phase II: Evaluation & diagnostics
# =============================================================================

# -----------------------------------------------------------------------------
# Step 4: Robust evaluation across held-out seeds
#
# Evaluates the final (tuned covariates + tuned hyperparams) on
# eval_fold_seed_list, which is disjoint from robust_fold_seed_list.
# -----------------------------------------------------------------------------
cat("\t\tStep 4: Robust evaluation across eval seeds\n")
assign("cv_type", "pixel_grouped", envir = .GlobalEnv)
source("modelling/multiseed/robust_evaluation.R")
cat("\n")

# -----------------------------------------------------------------------------
# Step 5: R2 sensitivity analysis (performance vs seed count, fold count, etc.)
# -----------------------------------------------------------------------------
if (isTRUE(do_sensitivity)) {
  cat("\t\tStep 5: R2 sensitivity analysis\n")
  source("modelling/analysis/sensitivity_suite.R")
  source("modelling/plots/plot_sensitivity_suite.R")
  cat("\n")
} else {
  cat("\t\tStep 5: Skipping R2 sensitivity analysis (do_sensitivity=FALSE)\n\n")
}

# -----------------------------------------------------------------------------
# Step 6: Train/test fraction & target-variance diagnostic
# -----------------------------------------------------------------------------
if (isTRUE(do_diagnostics)) {
  cat("\t\tStep 6: Train/test fraction & target-variance diagnostics\n")
  perf_detailed_override <- file.path(run_output_dir, "by_seed_detailed.csv")
  assign("perf_detailed_csv_override", perf_detailed_override, envir = .GlobalEnv)
  assign("fold_seed_list", eval_fold_seed_list, envir = .GlobalEnv)
  source("modelling/multiseed/train_test_fraction_diagnostic.R")
  cat("\n")
} else {
  cat("\t\tStep 6: Skipping diagnostics (do_diagnostics=FALSE)\n\n")
}

# -----------------------------------------------------------------------------
# Step 7: Fit and save final models using robust tuning/covariates
# -----------------------------------------------------------------------------
cat("\t\tStep 7: Fit and save final models (robust configs/vars)\n")
assign("use_robust_final_configs", TRUE, envir = .GlobalEnv)
source("modelling/R/fit_final_models.R")
cat("\n")

# -----------------------------------------------------------------------------
# Step 8: Directly compare model outputs with baselines
# -----------------------------------------------------------------------------
cat("\t\tStep 8: Model comparison\n")
source("modelling/analysis/model_comparison.R")
cat("\n")

# -----------------------------------------------------------------------------
# Step 9: Create model prediction maps
# -----------------------------------------------------------------------------
cat("\t\tStep 9: Create spatial prediction maps\n")
source("modelling/plots/plot_spatial_prediction_maps.R")
cat("\n")

# -----------------------------------------------------------------------------
# Step 10: Partial dependence plots
# -----------------------------------------------------------------------------
cat("\t\tStep 10: Partial dependence plots (from final models)\n")
source("modelling/plots/plot_partial_dependence.R")
cat("\n")

# -----------------------------------------------------------------------------
# Step 11: Supplement figures
# -----------------------------------------------------------------------------
cat("\t\tStep 11: Supplement figures\n")
source("modelling/plots/supplement.R")
cat("\n")


# =============================================================================
# Final summary
# =============================================================================
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("ROBUST PIXEL_GROUPED PIPELINE COMPLETE\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")
cat("Seed policy:\n")
cat("  Selection/tuning seeds: ", paste(robust_fold_seed_list, collapse = ", "), "\n")
cat("  Evaluation seeds:       ", paste(eval_fold_seed_list, collapse = ", "), "\n")
cat("  Separation:             optimise on selection seeds, evaluate on eval seeds\n\n")
cat("Outputs (cv_regime = ", cv_regime_name, ", run root = ", run_output_dir, "):\n", sep = "")
cat("  covariate_selection/ – Baseline + robust SHAP pruning\n", sep = "")
cat("  cv_pipeline/         – Robust tuning + related CV artefacts\n", sep = "")
cat("  final_models/, predictions/, analysis/, sensitivity_suite/, etc. – same run root\n", sep = "")
cat("\n")
