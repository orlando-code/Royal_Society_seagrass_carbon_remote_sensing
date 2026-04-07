# Fit and save final models (GPR, GAM, XGB, LR) on all training data.
#
# For each model:
#   1. Tune key hyperparameters via 5-fold random CV on the full training set.
#   2. Fit the best configuration on all data.
#   3. Save a shareable RDS to output/<cv_regime>/final_models/<model>_final.rds
#
# Each saved object is a named list with:
#   $model            – the fitted (tuned) model object
#   $predictor_vars   – character vector of predictors used (best covariate set)
#   $best_covariate_set – same as predictor_vars (explicit best feature set)
#   $hyperparams      – list of tuned hyperparameter values (from tuning pipeline)
#   $cv_metrics       – data.frame of tuning fold metrics when loaded from tuning
#   $train_metrics    – r2 / rmse on the full training set (in-sample)
#   $scale_params     – list(means, sds) for z-score re-prediction
#   $encoding         – categorical encoding (levels) for factor predictors used
#   $encoded_names    – same as predictor_vars (no one-hot expansion)
#
# Usage: sourced from run_paper.R (step 5), or run standalone.

# rm(list = ls())
if (!exists("seagrass_init_repo", mode = "function", inherits = TRUE)) source("modelling/R/init_repo.R")
project_root <- seagrass_init_repo(
  packages = c("here", "mgcv", "tidyverse", "randomForest", "GauPro", "xgboost", "sf"),
  source_files = c(
    "modelling/R/extract_covariates_from_rasters.R",
    "modelling/R/assign_region_from_latlon.R",
    "modelling/pipeline_config.R"
  ),
  include_helpers = TRUE,
  require_core_inputs = TRUE,
  check_renv = TRUE
)
set.seed(42)

cfg <- get_pipeline_config()
apply_pipeline_defaults(
  cfg,
  c(
    "cv_output_dir", "target_var", "use_shap_per_model", "log_transform_target",
    "n_folds", "cv_type", "cv_blocksize", "exclude_regions",
    "robust_pruned_csv_override", "use_robust_final_configs",
    "include_seagrass_species"
  ),
  envir = .GlobalEnv
)

cv_out  <- get("cv_output_dir", envir = .GlobalEnv)
out_dir <- file.path(cv_out, "final_models")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
use_robust_final_configs <- isTRUE(get("use_robust_final_configs", envir = .GlobalEnv))
robust_tune_dir <- file.path(cv_out, "cv_pipeline", "robust_pixel_grouped_tuning")

# ---------------------------------------------------------------------------
# Data (same pipeline as cv_pipeline / covariate_pruning_pipeline)
# ---------------------------------------------------------------------------
dat <- readr::read_rds("data/all_extracted_new.rds")
if (exists("exclude_regions", envir = .GlobalEnv)) {
  excl <- get("exclude_regions", envir = .GlobalEnv)
  if (length(excl) > 0L) {
    if (!"region" %in% names(dat)) dat <- assign_region_from_latlon(dat)
    dat <- dat[is.na(dat$region) | !dat$region %in% excl, ]
    cat("Excluded region(s):", paste(excl, collapse = ", "), "\n")
  }
}

target_var <- get("target_var", envir = .GlobalEnv)
predictor_vars <- raster_covariates[raster_covariates %in% colnames(dat)]
include_seagrass_species <- isTRUE(get("include_seagrass_species", envir = .GlobalEnv))
extra_cols <- c(if (include_seagrass_species) "seagrass_species" else character(0), "region")

core_data <- dat %>%
  dplyr::mutate(median_carbon_density = .data[[target_var]]) %>%
  dplyr::select(longitude, latitude, median_carbon_density,
                dplyr::all_of(predictor_vars),
                dplyr::all_of(intersect(extra_cols, names(dat))))
core_data <- as.data.frame(
  core_data %>%
    dplyr::select(longitude, latitude, median_carbon_density,
                  dplyr::all_of(predictor_vars),
                  dplyr::all_of(intersect(extra_cols, names(core_data)))) %>%
    dplyr::filter(complete.cases(.))
)
predictor_vars <- predictor_vars[predictor_vars %in% colnames(core_data)]

cov_dir <- file.path(cv_out, "covariate_selection")

# Load shared and per-model pruned covariates (helpers: get_per_model_vars, load_model_vars)
shared_file <- file.path(cov_dir, "pruned_variables_to_include.csv")
use_shap_per_model <- isTRUE(get("use_shap_per_model", envir = .GlobalEnv))

# Shared covariates (fallback)
shared_vars <- if (file.exists(shared_file)) {
  pv <- intersect(read.csv(shared_file, stringsAsFactors = FALSE)$variable, colnames(core_data))
  if (length(pv) >= 2) pv else predictor_vars
} else predictor_vars

per_model_vars <- get_per_model_vars(cov_dir, colnames(core_data), use_shap_first = use_shap_per_model)
robust_pruned_csv_override <- get("robust_pruned_csv_override", envir = .GlobalEnv)
robust_vars_by_model <- NULL
if (!is.na(robust_pruned_csv_override) && nzchar(robust_pruned_csv_override) && file.exists(robust_pruned_csv_override)) {
  robust_df <- read.csv(robust_pruned_csv_override, stringsAsFactors = FALSE)
  if (all(c("model", "variable") %in% names(robust_df))) {
    robust_vars_by_model <- lapply(split(robust_df$variable, robust_df$model), function(v) intersect(v, colnames(core_data)))
    cat("Loaded robust pruned covariates from:", robust_pruned_csv_override, "\n")
  }
}

load_model_vars_final <- function(model_name) {
  if (is.list(robust_vars_by_model) && model_name %in% names(robust_vars_by_model)) {
    vv <- robust_vars_by_model[[model_name]]
    if (!is.null(vv) && length(vv) >= 2L) return(vv)
  }
  load_model_vars(model_name, per_model_vars, use_shap_first = use_shap_per_model)
}

# Model/general settings
log_response   <- isTRUE(get("log_transform_target", envir = .GlobalEnv))
n_tune_folds   <- as.integer(get("n_folds", envir = .GlobalEnv))
cv_type        <- get("cv_type", envir = .GlobalEnv)
cv_blocksize  <- get("cv_blocksize", envir = .GlobalEnv)
exclude_regions <- get("exclude_regions", envir = .GlobalEnv)

cat("Training rows:", nrow(core_data),
    " | Shared predictors:", length(shared_vars),
    " | log_response:", log_response, "\n")

# Folds for fallback tuning: same cv_type as the rest of the pipeline
cv_fold_info <- make_cv_folds(
  core_data, predictor_vars, n_tune_folds, cv_type,
  cv_blocksize = cv_blocksize, exclude_regions = exclude_regions,
  cache_tag = "fit_final_fallback"
)
tune_folds <- cv_fold_info$fold_indices

# ---------------------------------------------------------------------------
# XGBoost: load best config from hyperparameter_tuning_pipeline.R
# ---------------------------------------------------------------------------
cat("=== XGBoost: hyperparameter tuning ===\n")
predictor_vars <- load_model_vars_final("XGB")
xgb_pvars <- predictor_vars
cat("  Predictors (", length(predictor_vars), "): ", paste(predictor_vars, collapse = ", "), "\n", sep = "")
config_dir <- file.path(cv_out, "cv_pipeline")
xgb_cfg <- load_best_model_config(
  model_name = "xgb",
  config_dir = config_dir,
  robust_config_dir = robust_tune_dir,
  prefer_robust = isTRUE(use_robust_final_configs),
  include_baseline = TRUE,
  include_legacy = TRUE
)
best_xgb_hp <- model_hyperparams_from_config("XGB", xgb_cfg)
xgb_cv_metrics <- NULL
if (!is.null(xgb_cfg)) {
  xgb_cv_metrics <- xgb_cfg$cv_metrics
  cat("  Loaded best config: nrounds =", best_xgb_hp$nrounds,
      ", max_depth =", best_xgb_hp$max_depth,
      ", lr =", best_xgb_hp$learning_rate,
      ", mcw =", best_xgb_hp$min_child_weight,
      ", min_split_loss =", best_xgb_hp$min_split_loss,
      ", reg_lambda =", best_xgb_hp$reg_reg_lambda, "\n\n")
} else {
  cat("  No saved XGB config found; using default XGB hyperparameters.\n")
}


xgb_data <- if (log_response) transform_response(core_data, "median_carbon_density", log = TRUE) else core_data
xgb_prep <- prepare_data_for_model("XGB", xgb_data, xgb_data, predictor_vars)
xgb_final <- fit_xgboost(xgb_prep$train, xgb_prep$test, xgb_prep$predictor_vars, best_xgb_hp)
xgb_pred <- xgb_final$predictions
if (log_response) xgb_pred <- inverse_response_transform(xgb_pred, log = TRUE)
xgb_train_metrics <- calculate_metrics(core_data$median_carbon_density, xgb_pred)
saveRDS(list(
  model            = xgb_final$model,
  predictor_vars   = predictor_vars,
  best_covariate_set = predictor_vars,
  encoded_names    = xgb_prep$predictor_vars,
  hyperparams      = best_xgb_hp,
  cv_metrics       = xgb_cv_metrics,
  train_metrics    = xgb_train_metrics,
  scale_params     = xgb_prep$scale_params,
  encoding         = xgb_prep$encoding,
  log_response     = log_response
), file.path(out_dir, "XGB_final.rds"))
cat("Saved XGB_final.rds  (train R2=", round(xgb_train_metrics$r2, 3),
    " RMSE=", round(xgb_train_metrics$rmse, 4), ")\n\n")

# ---------------------------------------------------------------------------
# GAM: load best config from hyperparameter_tuning_pipeline.R
# ---------------------------------------------------------------------------
cat("=== GAM: hyperparameter tuning ===\n")
predictor_vars <- load_model_vars_final("GAM")
gam_pvars <- predictor_vars
cat("  Predictors (", length(predictor_vars), "): ", paste(predictor_vars, collapse = ", "), "\n", sep = "")
gam_cfg <- load_best_model_config(
  model_name = "gam",
  config_dir = config_dir,
  robust_config_dir = robust_tune_dir,
  prefer_robust = isTRUE(use_robust_final_configs),
  include_baseline = TRUE,
  include_legacy = TRUE
)
best_gam_hp <- model_hyperparams_from_config("GAM", gam_cfg)
gam_cv_metrics <- NULL
if (!is.null(gam_cfg)) {
  gam_cv_metrics <- gam_cfg$cv_metrics
  cat("  Loaded best config: k_covariate =", best_gam_hp$k_covariate, "\n\n")
} else {
  cat("  No saved GAM config found; using default GAM hyperparameters.\n")
}


gam_data <- if (log_response) transform_response(core_data, "median_carbon_density", log = TRUE) else core_data
gam_prep <- prepare_data_for_model("GAM", gam_data, gam_data, predictor_vars)
gam_final <- fit_gam(gam_prep$train, gam_prep$test, predictor_vars,
                     include_spatial = FALSE, k_covariate = best_gam_hp$k_covariate)
gam_pred <- gam_final$predictions
if (log_response) gam_pred <- inverse_response_transform(gam_pred, log = TRUE)
gam_train_metrics <- calculate_metrics(core_data$median_carbon_density, gam_pred)
saveRDS(list(
  model              = gam_final$model,
  predictor_vars     = predictor_vars,
  best_covariate_set = predictor_vars,
  hyperparams        = best_gam_hp,
  cv_metrics         = gam_cv_metrics,
  train_metrics      = gam_train_metrics,
  scale_params       = gam_prep$scale_params,
  log_response       = log_response
), file.path(out_dir, "GAM_final.rds"))
cat("Saved GAM_final.rds  (train R2=", round(gam_train_metrics$r2, 3),
    " RMSE=", round(gam_train_metrics$rmse, 4), ")\n\n")

# ---------------------------------------------------------------------------
# GPR: use best config from hyperparameter_tuning_pipeline.R
# ---------------------------------------------------------------------------
cat("=== GPR: loading best config and fitting final model ===\n")
predictor_vars <- load_model_vars_final("GPR")
gpr_pvars <- predictor_vars
cat("  Predictors (", length(predictor_vars), "): ", paste(predictor_vars, collapse = ", "), "\n", sep = "")
gpr_cfg <- load_best_model_config(
  model_name = "gpr",
  config_dir = config_dir,
  robust_config_dir = robust_tune_dir,
  prefer_robust = isTRUE(use_robust_final_configs),
  include_baseline = TRUE,
  include_legacy = TRUE
)
gpr_hp <- model_hyperparams_from_config("GPR", gpr_cfg)
gpr_cv_metrics <- NULL
if (!is.null(gpr_cfg)) {
  gpr_cv_metrics <- gpr_cfg$cv_metrics
  cat("  Loaded best GPR config: kernel=", gpr_hp$kernel,
      " nug_max=", gpr_hp$nug.max, "\n")
} else {
  cat("  No saved GPR config found; using default GPR hyperparameters.\n")
}
gpr_data <- if (log_response) transform_response(core_data, "median_carbon_density", log = TRUE) else core_data
gpr_final_fit <- fit_gpr(
  train_data = gpr_data,
  test_data = gpr_data,
  predictor_vars = predictor_vars,
  hyperparams = gpr_hp
)
gpr_pred <- gpr_final_fit$predictions
if (log_response) gpr_pred <- inverse_response_transform(gpr_pred, log = TRUE)
gpr_train_metrics <- calculate_metrics(core_data$median_carbon_density, gpr_pred)
saveRDS(list(
  model              = gpr_final_fit$model,
  predictor_vars     = predictor_vars,
  best_covariate_set = predictor_vars,
  hyperparams        = gpr_hp,
  scale_params       = gpr_final_fit$scale_params,
  encoding           = gpr_final_fit$encoding,
  encoded_names      = gpr_final_fit$encoded_names,
  cv_metrics        = gpr_cv_metrics,
  train_metrics      = gpr_train_metrics,
  log_response       = log_response
), file.path(out_dir, "GPR_final.rds"))
cat("Saved GPR_final.rds  (train R2=", round(gpr_train_metrics$r2, 3),
    " RMSE=", round(gpr_train_metrics$rmse, 4), ")\n\n")

# ---------------------------------------------------------------------------
# LR baseline: no hyperparameter tuning (uses LR-specific covariates)
# ---------------------------------------------------------------------------
cat("=== LR: fitting linear baseline model ===\n")
predictor_vars <- load_model_vars_final("LR")
lm_pvars <- predictor_vars
cat("  Predictors (", length(predictor_vars), "): ", paste(predictor_vars, collapse = ", "), "\n", sep = "")

lm_cfg <- load_best_model_config(
  model_name = "lr",
  config_dir = config_dir,
  robust_config_dir = robust_tune_dir,
  prefer_robust = isTRUE(use_robust_final_configs),
  include_baseline = TRUE,
  include_legacy = TRUE
)
lm_cv_metrics <- if (!is.null(lm_cfg)) lm_cfg$cv_metrics else NULL

lm_data <- if (log_response) transform_response(core_data, "median_carbon_density", log = TRUE) else core_data
lm_prep <- prepare_data_for_model("LR", lm_data, lm_data, predictor_vars)
lm_final <- fit_lm(lm_prep$train, lm_prep$test, lm_prep$predictor_vars)
lm_pred <- lm_final$predictions
if (log_response) lm_pred <- inverse_response_transform(lm_pred, log = TRUE)
lm_train_metrics <- calculate_metrics(core_data$median_carbon_density, lm_pred)
saveRDS(list(
  model              = lm_final$model,
  predictor_vars     = predictor_vars,
  best_covariate_set = predictor_vars,
  hyperparams        = list(),
  cv_metrics         = lm_cv_metrics,
  train_metrics      = lm_train_metrics,
  scale_params       = lm_prep$scale_params,
  log_response       = log_response
), file.path(out_dir, "LR_final.rds"))
cat("Saved LR_final.rds  (train R2=", round(lm_train_metrics$r2, 3),
    " RMSE=", round(lm_train_metrics$rmse, 4), ")\n\n")

# ---------------------------------------------------------------------------
# Summary table (with variables included per model)
# ---------------------------------------------------------------------------
summary_tbl <- data.frame(
  model        = c("XGB", "GAM", "GPR", "LR"),
  train_r2     = round(c(xgb_train_metrics$r2,  gam_train_metrics$r2,  gpr_train_metrics$r2, lm_train_metrics$r2),  3),
  train_rmse   = round(c(xgb_train_metrics$rmse, gam_train_metrics$rmse, gpr_train_metrics$rmse, lm_train_metrics$rmse), 4),
  best_hp      = c(
    paste0("nrounds=", best_xgb_hp$nrounds, " max_depth=", best_xgb_hp$max_depth),
    paste0("k_covariate=", best_gam_hp$k_covariate),
    paste0("kernel=", gpr_hp$kernel, " nug_max=", gpr_hp$nug.max),
    "none (linear baseline)"
  ),
  variables    = c(
    paste(xgb_pvars, collapse = ", "),
    paste(gam_pvars, collapse = ", "),
    paste(gpr_pvars, collapse = ", "),
    paste(lm_pvars, collapse = ", ")
  ),
  stringsAsFactors = FALSE
)
write.csv(summary_tbl, file.path(out_dir, "final_models_summary.csv"), row.names = FALSE)
cat("Summary saved to", file.path(out_dir, "final_models_summary.csv"), "\n")
print(summary_tbl)
cat("\nAll final models saved to", out_dir, "\n")
cat("To re-use: obj <- readRDS('", file.path(out_dir, "XGB_final.rds"), "')\n", sep = "")
cat("           predict(obj$model, xgb.DMatrix(as.matrix(new_data[, obj$predictor_vars])))\n\n")
