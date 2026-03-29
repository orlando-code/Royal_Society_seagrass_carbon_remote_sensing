# Unified hyperparameter tuning for all models (GPR, GAM, XGB).
#
# Tunes GPR (kernel + nugget grid), XGB (random grid search), and GAM
# (k_covariate grid) using the same CV folds derived from the global
# cv_type setting. Saves best_config_<model>.rds to the cv_pipeline dir.
#
# Requires: covariate_pruning_pipeline.R has been run (pruned variables)
# Outputs:  output/<cv_regime>/cv_pipeline/best_config_{gpr,xgb,gam}.rds
#
# Usage: source after step 1 in run_paper.R, or run standalone.

project_root <- here::here()
source(file.path(project_root, "modelling/R/helpers.R"))
load_packages(c("here", "mgcv", "dplyr", "randomForest", "GauPro", "xgboost", "sf"))

out_dir <- file.path(get0("cv_output_dir", envir = .GlobalEnv, ifnotfound = "output"), "cv_pipeline")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Data and per-model vars (same logic as fit_final_models)
# ---------------------------------------------------------------------------
dat <- readr::read_rds(file.path(project_root, "data/all_extracted_new.rds"))
exclude_regions <- get0("exclude_regions", envir = .GlobalEnv, ifnotfound = character(0))
if (length(exclude_regions) > 0L) {
  if (!"region" %in% names(dat)) dat <- assign_region_from_latlon(dat)
  dat <- dat[is.na(dat$region) | !dat$region %in% exclude_regions, ]
}
target_var <- get0("target_var", envir = .GlobalEnv, ifnotfound = "median_carbon_density_100cm")
log_transform_target <- get0("log_transform_target", envir = .GlobalEnv, ifnotfound = TRUE)
predictor_vars <- raster_covariates[raster_covariates %in% colnames(dat)]
include_seagrass_species <- isTRUE(get0("include_seagrass_species", envir = .GlobalEnv, ifnotfound = TRUE))
extra_cols <- c(if (isTRUE(include_seagrass_species)) "seagrass_species" else character(0), "region")

core_data <- dat %>%
  dplyr::select(longitude, latitude, target_var, dplyr::all_of(predictor_vars),
                dplyr::all_of(intersect(extra_cols, names(dat)))) %>%
  dplyr::filter(complete.cases(.))
core_data$median_carbon_density <- core_data[[target_var]]
predictor_vars <- predictor_vars[predictor_vars %in% colnames(core_data)]

if (isTRUE(include_seagrass_species) && "seagrass_species" %in% names(core_data)) {
  cat("Including seagrass_species as factor in hyperparameter tuning.\n")
}

n_folds       <- as.integer(get0("n_folds", envir = .GlobalEnv, ifnotfound = 5L))
cv_type       <- get0("cv_type", envir = .GlobalEnv, ifnotfound = "spatial")
cv_blocksize <- get0("cv_blocksize", envir = .GlobalEnv, ifnotfound = 5000L)
log_response  <- isTRUE(get0("log_transform_target", envir = .GlobalEnv, ifnotfound = TRUE))
model_list    <- get0("model_list", envir = .GlobalEnv, ifnotfound = c("GPR", "GAM", "XGB"))

cov_dir <- file.path(get0("cv_output_dir", envir = .GlobalEnv, ifnotfound = "output"), "covariate_selection")
use_shap_per_model <- isTRUE(get0("use_shap_per_model", envir = .GlobalEnv, ifnotfound = FALSE))
per_model_vars <- get_per_model_vars(cov_dir, colnames(core_data), use_shap_first = use_shap_per_model)


# Build folds from the global cv_type setting (run_paper.R)
cv_fold_info <- make_cv_folds(
  core_data, predictor_vars, n_folds, cv_type,
  cv_blocksize = cv_blocksize, exclude_regions = exclude_regions,
  cache_tag = "hp_tuning"
)
block_size   <- NULL
fold_indices <- cv_fold_info$fold_indices


# ---------------------------------------------------------------------------
# GPR
# ---------------------------------------------------------------------------
if ("GPR" %in% model_list) {
  cat("\n=== GPR hyperparameter tuning (progress below) ===\n")
  pvars <- load_model_vars("GPR", per_model_vars, use_shap_first = use_shap_per_model)
  cat("  Predictor(s):", length(pvars), "\n")
  best_gpr <- tune_gpr_cv(
    train_data = core_data,
    predictor_vars = pvars,
    n_folds = n_folds,
    verbose = TRUE,
    value_var = target_var,
    fold_indices = fold_indices,
    block_size = block_size,
    cache_tag = "tune_gpr_cv",
    exclude_regions = exclude_regions
  )
  gpr_config <- list(
    kernel = best_gpr$kernel,
    nug.min = best_gpr$nug.min,
    nug.max = best_gpr$nug.max,
    nug.est = TRUE,
    cv_metrics = best_gpr$cv_metrics
  )
  saveRDS(gpr_config, file.path(out_dir, "best_config_gpr.rds"))
  cat("  Saved best_config_gpr.rds: kernel =", gpr_config$kernel,
      ", nug.min =", format(gpr_config$nug.min, scientific = TRUE),
      ", nug.max =", gpr_config$nug.max, "\n")
}

# ---------------------------------------------------------------------------
# XGBoost (grid with log_response and back-transform for metrics, same as fit_final_models)
# ---------------------------------------------------------------------------
if ("XGB" %in% model_list) {
  cat("\n=== XGBoost hyperparameter tuning ===\n")
  pvars <- load_model_vars("XGB", per_model_vars, use_shap_first = use_shap_per_model)
  cat("  Predictor(s):", length(pvars), "\n")
  best_xgb <- tune_xgb_cv(
    core_data = core_data,
    predictor_vars = pvars,
    fold_indices = fold_indices,
    log_response = log_response,
    verbose = TRUE
  )
  xgb_config <- list(
    nrounds = best_xgb$nrounds,
    max_depth = best_xgb$max_depth,
    learning_rate = best_xgb$learning_rate,
    subsample = best_xgb$subsample,
    colsample_bytree = best_xgb$colsample_bytree,
    min_child_weight = best_xgb$min_child_weight,
    min_split_loss = best_xgb$min_split_loss,
    reg_reg_lambda = best_xgb$reg_reg_lambda,
    cv_metrics = best_xgb$cv_metrics
  )
  saveRDS(xgb_config, file.path(out_dir, "best_config_xgb.rds"))
  cat("  Saved best_config_xgb.rds: nrounds =", xgb_config$nrounds,
      ", max_depth =", xgb_config$max_depth,
      ", lr =", xgb_config$learning_rate,
      ", mcw =", xgb_config$min_child_weight,
      ", min_split_loss =", xgb_config$min_split_loss,
      ", reg_lambda =", xgb_config$reg_reg_lambda, "\n")
}

# ---------------------------------------------------------------------------
# GAM (k_covariate grid — basis dimension for smooth covariate terms)
# ---------------------------------------------------------------------------
if ("GAM" %in% model_list) {
  cat("\n=== GAM hyperparameter tuning (k_covariate, no spatial smooth) ===\n")
  pvars <- load_model_vars("GAM", per_model_vars, use_shap_first = use_shap_per_model)
  cat("  Predictor(s):", length(pvars), "\n")
  best_gam <- tune_gam_cv(
    core_data = core_data,
    predictor_vars = pvars,
    fold_indices = fold_indices,
    log_response = log_response,
    verbose = TRUE
  )
  gam_config <- list(k_covariate = best_gam$k_covariate, cv_metrics = best_gam$cv_metrics)
  saveRDS(gam_config, file.path(out_dir, "best_config_gam.rds"))
  cat("  Saved best_config_gam.rds: k_covariate =", best_gam$k_covariate, "\n")
}

cat("\nHyperparameter tuning complete. Best configs in", out_dir, "\n")

