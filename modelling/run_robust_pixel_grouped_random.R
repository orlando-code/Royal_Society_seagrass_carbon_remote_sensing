# =============================================================================
# Robust (across-seed) pixel_grouped_random investigation driver
#
# This is a seed-robust analogue of `modelling/run_paper.R`, focused on
# pixel_grouped_random. Multi-seed scripts now live under `modelling/multiseed/`.
#   1) baseline covariate pruning (for warm-start / fallbacks)
#   2) robust warm-start hyperparameter tuning across multiple fold seeds
#   3) robust SHAP-based covariate selection across fold seeds
#   4) robust hyperparameter tuning again using SHAP-selected covariates
#   5) robust evaluation across (potentially different) fold seeds
#   6) train/test fraction + target-variance diagnostic
# =============================================================================

setwd(here::here())
set.seed(42)

cat("\n")
cat(paste(rep("=", 94), collapse = ""), "\n")
cat("ROBUST PIXEL_GROUPED_RANDOM: CONFIG -> WARM-START -> SHAP-PRUNE -> ROBUST-TUNE -> ROBUST-EVAL\n")
cat(paste(rep("=", 94), collapse = ""), "\n\n")

# -----------------------------------------------------------------------------
# -1. Configuration
# -----------------------------------------------------------------------------
dpi <- 150
assign("dpi", dpi, envir = .GlobalEnv)

# Response
target_var <- "median_carbon_density_100cm"
log_transform_target <- TRUE
assign("target_var", target_var, envir = .GlobalEnv)
assign("log_transform_target", log_transform_target, envir = .GlobalEnv)

# Region exclusion
exclude_regions <- c("Black Sea")
assign("exclude_regions", exclude_regions, envir = .GlobalEnv)

# Base pruning configuration (for warm-start fallbacks)
use_correlation_filter <- TRUE
correlation_filter_threshold <- 0.8
permutation_max_vars <- 15L
n_permutations <- 1L
permutation_coverage <- 0.99

use_shap_per_model <- TRUE  # use SHAP-derived per-model vars in baseline tuning as well

assign("use_correlation_filter", use_correlation_filter, envir = .GlobalEnv)
assign("correlation_filter_threshold", correlation_filter_threshold, envir = .GlobalEnv)
assign("permutation_max_vars", permutation_max_vars, envir = .GlobalEnv)
assign("n_permutations", n_permutations, envir = .GlobalEnv)
assign("permutation_coverage", permutation_coverage, envir = .GlobalEnv)
assign("use_shap_per_model", use_shap_per_model, envir = .GlobalEnv)

# Models
model_list <- get0("model_list", envir = .GlobalEnv, ifnotfound = c("GPR", "GAM", "XGB"))
assign("model_list", model_list, envir = .GlobalEnv)

# Fold count and base CV type for warm-start tuning/pruning
n_folds <- 10L
cv_type <- "pixel_grouped"         # supported hashing type for make_cv_folds
cv_blocksize <- 1000L             # unused for pixel_grouped

assign("n_folds", n_folds, envir = .GlobalEnv)
assign("cv_type", cv_type, envir = .GlobalEnv)
assign("cv_blocksize", cv_blocksize, envir = .GlobalEnv)

cv_regime_name <- "pixel_grouped"
cv_output_dir <- file.path("output", cv_regime_name)
assign("cv_regime_name", cv_regime_name, envir = .GlobalEnv)
assign("cv_output_dir", cv_output_dir, envir = .GlobalEnv)

# -----------------------------------------------------------------------------
# Robust seeds: (selection/tuning) vs (final evaluation)
# -----------------------------------------------------------------------------
# You can set these equal if you prefer faster runs, but for a cleaner
# "optimize on A, evaluate on B" story, keep them different.
robust_fold_seed_list <- get0("robust_fold_seed_list", envir = .GlobalEnv, ifnotfound = 42L:50L)
eval_fold_seed_list <- get0("eval_fold_seed_list", envir = .GlobalEnv, ifnotfound = 51L:60L)

assign("robust_fold_seed_list", robust_fold_seed_list, envir = .GlobalEnv)
assign("eval_fold_seed_list", eval_fold_seed_list, envir = .GlobalEnv)

# -----------------------------------------------------------------------------
# SHAP pruning settings (controls runtime)
# -----------------------------------------------------------------------------
shap_n_points <- get0("shap_n_points", envir = .GlobalEnv, ifnotfound = 50L)            # increase for more stable SHAP
shap_folds_per_seed <- get0("shap_folds_per_seed", envir = .GlobalEnv, ifnotfound = 10L)      # 1 = faster, fewer fold-specific SHAP evals
shap_max_gpr_train <- get0("shap_max_gpr_train", envir = .GlobalEnv, ifnotfound = 400L)     # cap for SHAP GPR

assign("shap_n_points", shap_n_points, envir = .GlobalEnv)
assign("shap_folds_per_seed", shap_folds_per_seed, envir = .GlobalEnv)
assign("shap_max_gpr_train", shap_max_gpr_train, envir = .GlobalEnv)

# -----------------------------------------------------------------------------
# Shared regime variables for the robust scripts
# -----------------------------------------------------------------------------
cv_type_label <- "pixel_grouped_random"
assign("cv_type_label", cv_type_label, envir = .GlobalEnv)

# -----------------------------------------------------------------------------
# 0. Data
# -----------------------------------------------------------------------------
cat("\n\n\tStep 0: Data (build if missing)\n")
if (!file.exists("data/all_extracted_new.rds")) {
  source("modelling/pipeline/build_all_extracted_new.R")
}

source("modelling/R/extract_covariates_from_rasters.R")

# -----------------------------------------------------------------------------
# 1. Baseline covariate pruning (for warm-start fallbacks + best_config_*.rds)
# -----------------------------------------------------------------------------
cat("\n\n\tStep 1: Covariate pruning (warm-start)\n")
assign("cv_output_dir", cv_output_dir, envir = .GlobalEnv)
assign("cv_type", "pixel_grouped", envir = .GlobalEnv)
source("modelling/pipeline/covariate_pruning_pipeline.R")

# -----------------------------------------------------------------------------
# 2. Baseline hyperparameter tuning (creates best_config_*.rds)
# -----------------------------------------------------------------------------
cat("\n\n\tStep 2: Baseline hyperparameter tuning (cv_type=pixel_grouped)\n")
source("modelling/pipeline/hyperparameter_tuning_pipeline.R")

# -----------------------------------------------------------------------------
# 3. Robust warm-start tuning across fold seeds (SHAP variable sets)
# -----------------------------------------------------------------------------
cat("\n\n\tStep 3: Robust warm-start hyperparameter tuning\n")
assign("cv_type", "pixel_grouped_random", envir = .GlobalEnv)
assign("robust_pruned_importance_type", "shap", envir = .GlobalEnv)
source("modelling/multiseed/robust_hyperparameter_tuning.R")

# -----------------------------------------------------------------------------
# 4. Robust SHAP-based covariate selection across fold seeds
# -----------------------------------------------------------------------------
cat("\n\n\tStep 4: Robust SHAP covariate pruning\n")
source("modelling/multiseed/robust_shap_covariate_pruning.R")

# -----------------------------------------------------------------------------
# 5. Robust tuning again using SHAP-selected covariates
# -----------------------------------------------------------------------------
cat("\n\n\tStep 5: Robust hyperparameter tuning using SHAP-selected covariates\n")
assign("robust_pruned_importance_type", "shap", envir = .GlobalEnv)
source("modelling/multiseed/robust_hyperparameter_tuning.R")

# -----------------------------------------------------------------------------
# 6. Robust evaluation across fold seeds
# -----------------------------------------------------------------------------
cat("\n\n\tStep 6: Robust evaluation across seeds\n")
assign("cv_type", "pixel_grouped_random", envir = .GlobalEnv)
assign("robust_pruned_importance_type", "shap", envir = .GlobalEnv)
source("modelling/multiseed/robust_evaluation.R")

# -----------------------------------------------------------------------------
# 7. Train/test fraction diagnostic joined to robust evaluation results
# -----------------------------------------------------------------------------
cat("\n\n\tStep 7: Train/test fraction & target-variance diagnostics\n")

seeds_sel_str <- paste(robust_fold_seed_list, collapse = "-")
perf_detailed_override <- file.path(
  "output", cv_regime_name, "cv_pipeline",
  paste0("robust_pixel_grouped_random_evaluation_robustSeeds_", seeds_sel_str),
  "by_seed_detailed.csv"
)

assign("perf_detailed_csv_override", perf_detailed_override, envir = .GlobalEnv)

assign("fold_seed_list", eval_fold_seed_list, envir = .GlobalEnv)
source("modelling/multiseed/train_test_fraction_diagnostic.R")

cat("\nRobust pixel_grouped_random driver finished.\n")

