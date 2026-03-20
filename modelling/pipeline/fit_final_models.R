# Fit and save final models (GPR, GAM, XGB) on all training data.
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
#   $encoding         – categorical encoding (levels) for seagrass_species, region
#   $encoded_names    – same as predictor_vars (no one-hot expansion)
#
# Usage: sourced from run_paper.R (step 5), or run standalone.

# rm(list = ls())
setwd(here::here())
set.seed(42)

source("modelling/R/helpers.R")
source("modelling/R/extract_covariates_from_rasters.R")
source("modelling/R/assign_region_from_latlon.R")
load_packages(c("here", "mgcv", "tidyverse", "randomForest", "GauPro", "xgboost", "sf"))

cv_out  <- get0("cv_output_dir", envir = .GlobalEnv, ifnotfound = "output")
out_dir <- file.path(cv_out, "final_models")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

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

target_var <- get0("target_var", envir = .GlobalEnv, ifnotfound = "median_carbon_density_100cm")
predictor_vars <- raster_covariates[raster_covariates %in% colnames(dat)]

core_data <- dat %>%
  dplyr::mutate(median_carbon_density = .data[[target_var]]) %>%
  dplyr::select(longitude, latitude, median_carbon_density,
                dplyr::all_of(predictor_vars),
                dplyr::all_of(intersect(c("seagrass_species", "region"), names(dat))))
core_data <- as.data.frame(
  core_data %>%
    dplyr::select(longitude, latitude, median_carbon_density,
                  dplyr::all_of(predictor_vars),
                  dplyr::all_of(intersect(c("seagrass_species", "region"), names(core_data)))) %>%
    dplyr::filter(complete.cases(.))
)
predictor_vars <- predictor_vars[predictor_vars %in% colnames(core_data)]

cov_dir <- file.path(cv_out, "covariate_selection")

# Load shared and per-model pruned covariates (helpers: get_per_model_vars, load_model_vars)
shared_file <- file.path(cov_dir, "pruned_variables_to_include.csv")
use_shap_per_model <- isTRUE(get0("use_shap_per_model", envir = .GlobalEnv, ifnotfound = FALSE))

# Shared covariates (fallback)
shared_vars <- if (file.exists(shared_file)) {
  pv <- intersect(read.csv(shared_file, stringsAsFactors = FALSE)$variable, colnames(core_data))
  if (length(pv) >= 2) pv else predictor_vars
} else predictor_vars

per_model_vars <- get_per_model_vars(cov_dir, colnames(core_data), use_shap_first = use_shap_per_model)

# Model/general settings
log_response   <- isTRUE(get0("log_transform_target", envir = .GlobalEnv, ifnotfound = TRUE))
n_tune_folds   <- as.integer(get0("n_folds", envir = .GlobalEnv, ifnotfound = 5L))
cv_type        <- get0("cv_type", envir = .GlobalEnv, ifnotfound = "spatial")
cv_blocksize  <- get0("cv_blocksize", envir = .GlobalEnv, ifnotfound = 5000L)
exclude_regions <- get0("exclude_regions", envir = .GlobalEnv, ifnotfound = character(0))

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
predictor_vars <- load_model_vars("XGB", per_model_vars, use_shap_first = use_shap_per_model)
xgb_pvars <- predictor_vars
cat("  Predictors (", length(predictor_vars), "): ", paste(predictor_vars, collapse = ", "), "\n", sep = "")
config_dir <- file.path(cv_out, "cv_pipeline")
xgb_config_path <- file.path(config_dir, "best_config_xgb.rds")
if (!file.exists(xgb_config_path)) xgb_config_path <- file.path(config_dir, "xgb_best_config.rds")
if (file.exists(xgb_config_path)) {
  xgb_cfg <- readRDS(xgb_config_path)
  best_xgb_hp <- list(nrounds = xgb_cfg$nrounds, max_depth = xgb_cfg$max_depth,
    learning_rate = xgb_cfg$learning_rate %||% 0.1, subsample = xgb_cfg$subsample %||% 0.8,
    colsample_bytree = xgb_cfg$colsample_bytree %||% 0.8,
    min_child_weight = xgb_cfg$min_child_weight %||% 1L,
    min_split_loss = xgb_cfg$min_split_loss %||% 0,
    reg_reg_lambda = xgb_cfg$reg_reg_lambda %||% 1)
  xgb_cv_metrics <- xgb_cfg$cv_metrics
  cat("  Loaded best config: nrounds =", best_xgb_hp$nrounds,
      ", max_depth =", best_xgb_hp$max_depth,
      ", lr =", best_xgb_hp$learning_rate,
      ", mcw =", best_xgb_hp$min_child_weight,
      ", min_split_loss =", best_xgb_hp$min_split_loss,
      ", reg_lambda =", best_xgb_hp$reg_reg_lambda, "\n\n")
} else {
  cat("  No best config found; make sure to run hyperparameter_tuning_pipeline.R first).\n")}


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
predictor_vars <- load_model_vars("GAM", per_model_vars, use_shap_first = use_shap_per_model)
gam_pvars <- predictor_vars
cat("  Predictors (", length(predictor_vars), "): ", paste(predictor_vars, collapse = ", "), "\n", sep = "")
gam_config_path <- file.path(config_dir, "best_config_gam.rds")
if (!file.exists(gam_config_path)) gam_config_path <- file.path(config_dir, "gam_best_config.rds")
if (file.exists(gam_config_path)) {
  gam_cfg <- readRDS(gam_config_path)
  best_gam_hp <- list(k_covariate = gam_cfg$k_covariate)
  gam_cv_metrics <- gam_cfg$cv_metrics
  cat("  Loaded best config: k_covariate =", best_gam_hp$k_covariate, "\n\n")
} else {
  cat("  No best config found; make sure to run hyperparameter_tuning_pipeline.R first).\n")}


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
predictor_vars <- load_model_vars("GPR", per_model_vars, use_shap_first = use_shap_per_model)
gpr_pvars <- predictor_vars
cat("  Predictors (", length(predictor_vars), "): ", paste(predictor_vars, collapse = ", "), "\n", sep = "")
best_config_path <- file.path(config_dir, "best_config_gpr.rds")
if (!file.exists(best_config_path)) best_config_path <- file.path(config_dir, "gpr_best_config.rds")
if (file.exists(best_config_path)) {
  best_config <- readRDS(best_config_path)
  gpr_kernel  <- if (!is.null(best_config$kernel))  best_config$kernel  else "matern52"
  gpr_nug_max <- if (!is.null(best_config$nug_max)) best_config$nug_max else 100
  gpr_nug_min <- if (!is.null(best_config$nug_min)) best_config$nug_min else 1e-8
  gpr_cv_metrics <- best_config$cv_metrics
  cat("  Loaded best GPR config: kernel=", gpr_kernel,
      " nug_max=", gpr_nug_max, "\n")
} else {
  cat("  No best config found; make sure to run hyperparameter_tuning_pipeline.R first).\n")}

gpr_hp <- list(kernel = gpr_kernel, nug.min = gpr_nug_min, nug.max = gpr_nug_max)
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
# Summary table (with variables included per model)
# ---------------------------------------------------------------------------
summary_tbl <- data.frame(
  model        = c("XGB", "GAM", "GPR"),
  train_r2     = round(c(xgb_train_metrics$r2,  gam_train_metrics$r2,  gpr_train_metrics$r2),  3),
  train_rmse   = round(c(xgb_train_metrics$rmse, gam_train_metrics$rmse, gpr_train_metrics$rmse), 4),
  best_hp      = c(
    paste0("nrounds=", best_xgb_hp$nrounds, " max_depth=", best_xgb_hp$max_depth),
    paste0("k_covariate=", best_gam_hp$k_covariate),
    paste0("kernel=", gpr_kernel, " nug_max=", gpr_nug_max)
  ),
  variables    = c(
    paste(xgb_pvars, collapse = ", "),
    paste(gam_pvars, collapse = ", "),
    paste(gpr_pvars, collapse = ", ")
  ),
  stringsAsFactors = FALSE
)
write.csv(summary_tbl, file.path(out_dir, "final_models_summary.csv"), row.names = FALSE)
cat("Summary saved to", file.path(out_dir, "final_models_summary.csv"), "\n")
print(summary_tbl)
cat("\nAll final models saved to", out_dir, "\n")
cat("To re-use: obj <- readRDS('", file.path(out_dir, "XGB_final.rds"), "')\n", sep = "")
cat("           predict(obj$model, xgb.DMatrix(as.matrix(new_data[, obj$predictor_vars])))\n\n")
