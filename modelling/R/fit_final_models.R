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
    "use_robust_final_configs", "include_seagrass_species",
    "robust_fold_seed_list", "robust_pruned_importance_type", "model_list"
  ),
  envir = .GlobalEnv
)

out_dir <- file.path(cv_output_dir, "final_models")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
# use_robust_final_configs <- isTRUE(get("use_robust_final_configs", envir = .GlobalEnv))
# run_output_dir <- get0("run_output_dir", envir = .GlobalEnv, ifnotfound = NA_character_)
if (is.na(run_output_dir) || !nzchar(as.character(run_output_dir)) || is.null(run_output_dir)) {
  stop("run_output_dir must be set in .GlobalEnv for fit_final_models.R")
}
# run_output_dir <- as.character(run_output_dir)
robust_tune_dir <- file.path(run_output_dir, "cv_pipeline", "robust_pixel_grouped_tuning")
if (!dir.exists(robust_tune_dir)) {
  stop("Required robust tuning directory not found: ", robust_tune_dir)
}

# ---------------------------------------------------------------------------
# Data (same pipeline as baseline CV/covariate pruning)
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
# count number of unique rows in dat for raster_covariates columns
nrow(unique(dat[raster_covariates]))
nrow(unique(core_data[raster_covariates]))


predictor_vars <- raster_covariates[raster_covariates %in% colnames(dat)]
extra_cols <- c(if (include_seagrass_species) "seagrass_species" else character(0), "region")

selected_cols <- c(
  "longitude",
  "latitude",
  "median_carbon_density",
  predictor_vars,
  intersect(extra_cols, names(dat))
)
selected_cols <- unique(selected_cols)

core_data <- dat %>%
  dplyr::mutate(median_carbon_density = .data[[target_var]]) %>%
  dplyr::select(dplyr::all_of(selected_cols)) %>%
  dplyr::filter(complete.cases(.)) %>%
  as.data.frame()

predictor_vars <- intersect(predictor_vars, names(core_data))

cov_dir <- file.path(cv_output_dir, "covariate_selection")

# Load shared and per-model pruned covariates (helpers: get_per_model_vars, load_model_vars)
shared_file <- file.path(cov_dir, "pruned_variables_to_include.csv")

# Shared covariates (fallback)
shared_vars <- if (file.exists(shared_file)) {
  pv <- intersect(read.csv(shared_file, stringsAsFactors = FALSE)$variable, colnames(core_data))
  if (length(pv) >= 2) pv else predictor_vars
} else predictor_vars

per_model_vars <- get_per_model_vars(file.path(cov_dir, "robust_pixel_grouped"), colnames(core_data), use_shap_first = use_shap_per_model)
robust_vars_by_model <- NULL
seeds_str <- paste(as.integer(robust_fold_seed_list), collapse = "-")
importance_type <- match.arg(robust_pruned_importance_type, choices = c("perm", "shap"))
robust_cov_dir <- file.path(run_output_dir, "covariate_selection", "robust_pixel_grouped")
robust_pruned_csv <- file.path(
  robust_cov_dir,
  if (identical(importance_type, "shap")) {
    paste0("pruned_model_variables_shap_robust_pixel_grouped_seeds_", seeds_str, ".csv")
  } else {
    paste0("pruned_model_variables_perm_robust_pixel_grouped_seeds_", seeds_str, ".csv")
  }
)
if (!file.exists(robust_pruned_csv)) {
  stop("Required robust pruned covariate file not found: ", robust_pruned_csv)
}
robust_df <- read.csv(robust_pruned_csv, stringsAsFactors = FALSE)
if (all(c("model", "variable") %in% names(robust_df))) {
  robust_vars_by_model <- lapply(split(robust_df$variable, robust_df$model), function(v) intersect(v, colnames(core_data)))
  cat("Loaded robust pruned covariates from:", robust_pruned_csv, "\n")
}

load_model_vars_final <- function(model_name) {
  if (is.list(robust_vars_by_model) && model_name %in% names(robust_vars_by_model)) {
    vv <- robust_vars_by_model[[model_name]]
    if (!is.null(vv) && length(vv) >= 2L) return(vv)
  }
  load_model_vars(model_name, per_model_vars, use_shap_first = use_shap_per_model)
}

# Model/general settings
cv_types_allowed <- c("random", "location_grouped", "pixel_grouped", "spatial")
# cv_type_raw <- cfg$cv_type
if (length(cv_type) != 1L || !cv_type %in% cv_types_allowed) {
  if (isTRUE(use_robust_final_configs)) {
    message("Ignoring invalid global cv_type (", encodeString(as.character(cv_type)[1], quote = "\""),
            "); using \"pixel_grouped\" for fallback CV in robust final-model fitting.")
    cv_type <- "pixel_grouped"
    # assign("cv_type", cv_type, envir = .GlobalEnv)
  } else {
    stop(
      "Invalid cv_type in .GlobalEnv: ", paste(cv_type, collapse = ", "),
      ". Expected one of: ", paste(cv_types_allowed, collapse = ", "), ".",
      call. = FALSE
    )
  }
}

cat("Training rows:", nrow(core_data),
    " | Shared predictors:", length(shared_vars),
    " | log_transform_target:", log_transform_target, "\n")

# Folds for fallback tuning: same cv_type as the rest of the pipeline
cv_fold_info <- make_cv_folds(
  core_data, predictor_vars, n_folds, cv_type,
  cv_blocksize = cv_blocksize, exclude_regions = exclude_regions,
  cache_tag = "fit_final_fallback"
)
tune_folds <- cv_fold_info$fold_indices

# Log-transform median carbon density if required
if (log_transform_target) {
  final_data <- transform_response(core_data, "median_carbon_density", log = TRUE)
}
# ---------------------------------------------------------------------------
# XGBoost: load best config from tuning outputs
# ---------------------------------------------------------------------------
if ("XGB" %in% model_list) {
  cat("=== XGBoost: hyperparameter tuning ===\n")
  xgb_pvars <- load_model_vars_final("XGB")
  cat("  Predictors (", length(xgb_pvars), "): ", paste(xgb_pvars, collapse = ", "), "\n", sep = "")
  config_dir <- file.path(cv_output_dir, "cv_pipeline")
  xgb_cfg <- load_best_model_config(
    model_name = "xgb",
    config_dir = config_dir,
    robust_config_dir = robust_tune_dir,
    prefer_robust = isTRUE(use_robust_final_configs),
    include_baseline = TRUE,
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

  xgb_prep <- prepare_data_for_model("XGB", train=final_data, test=final_data, predictor_vars=xgb_pvars)
  xgb_final <- fit_xgboost(xgb_prep$train, xgb_prep$test, xgb_prep$predictor_vars, best_xgb_hp)
  xgb_pred <- xgb_final$predictions
  if (log_transform_target) xgb_pred <- inverse_response_transform(xgb_pred, log = TRUE)
  xgb_train_metrics <- calculate_metrics(core_data$median_carbon_density, xgb_pred)
  saveRDS(list(
    model            = xgb_final$model,
    predictor_vars   = xgb_pvars,
    best_covariate_set = xgb_pvars,
    encoded_names    = xgb_prep$xgb_pvars,
    hyperparams      = best_xgb_hp,
    cv_metrics       = xgb_cv_metrics,
    train_metrics    = xgb_train_metrics,
    scale_params     = xgb_prep$scale_params,
    encoding         = xgb_prep$encoding,
    log_response     = log_transform_target
  ), file.path(out_dir, "XGB_final.rds"))
  cat("Saved XGB_final.rds  (train R2=", round(xgb_train_metrics$r2, 3),
      " RMSE=", round(xgb_train_metrics$rmse, 4), ")\n\n")
}

# ---------------------------------------------------------------------------
# GAM: load best config from tuning outputs
# ---------------------------------------------------------------------------
if ("GAM" %in% model_list) {
  cat("=== GAM: hyperparameter tuning ===\n")
  gam_pvars <- load_model_vars_final("GAM")
  cat("  Predictors (", length(gam_pvars), "): ", paste(gam_pvars, collapse = ", "), "\n", sep = "")
  gam_cfg <- load_best_model_config(
    model_name = "gam",
    config_dir = config_dir,
    robust_config_dir = robust_tune_dir,
    prefer_robust = isTRUE(use_robust_final_configs),
    include_baseline = TRUE,
  )
  best_gam_hp <- model_hyperparams_from_config("GAM", gam_cfg)
  gam_cv_metrics <- NULL
  if (!is.null(gam_cfg)) {
    gam_cv_metrics <- gam_cfg$cv_metrics
    cat("  Loaded best config: k_covariate =", best_gam_hp$k_covariate, "\n\n")
  } else {
    cat("  No saved GAM config found; using default GAM hyperparameters.\n")
  }

  gam_prep <- prepare_data_for_model("GAM", train=final_data, test=final_data, predictor_vars=gam_pvars)
  gam_final <- fit_gam(gam_prep$train, gam_prep$test, gam_prep$predictor_vars,
                      include_spatial = FALSE, k_covariate = best_gam_hp$k_covariate)
  gam_pred <- gam_final$predictions
  if (log_transform_target) gam_pred <- inverse_response_transform(gam_pred, log = TRUE)
  gam_train_metrics <- calculate_metrics(core_data$median_carbon_density, gam_pred)
  saveRDS(list(
    model              = gam_final$model,
    predictor_vars     = gam_pvars,
    best_covariate_set = gam_pvars,
    hyperparams        = best_gam_hp,
    cv_metrics         = gam_cv_metrics,
    train_metrics      = gam_train_metrics,
    scale_params       = gam_prep$scale_params,
    log_response       = log_transform_target
  ), file.path(out_dir, "GAM_final.rds"))
  cat("Saved GAM_final.rds  (train R2=", round(gam_train_metrics$r2, 3),
      " RMSE=", round(gam_train_metrics$rmse, 4), ")\n\n")
}
# ---------------------------------------------------------------------------
# GPR: use best config from tuning outputs
# ---------------------------------------------------------------------------
if ("GPR" %in% model_list) {
  cat("=== GPR: loading best config and fitting final model ===\n")
  gpr_pvars <- load_model_vars_final("GPR")
  cat("  Predictors (", length(gpr_pvars), "): ", paste(gpr_pvars, collapse = ", "), "\n", sep = "")
  gpr_cfg <- load_best_model_config(
    model_name = "gpr",
    config_dir = config_dir,
    robust_config_dir = robust_tune_dir,
    prefer_robust = isTRUE(use_robust_final_configs),
    include_baseline = TRUE,
  )
  best_gpr_hp <- model_hyperparams_from_config("GPR", gpr_cfg)
  gpr_cv_metrics <- NULL
  if (!is.null(gpr_cfg)) {
    gpr_cv_metrics <- gpr_cfg$cv_metrics
    cat("  Loaded best GPR config: kernel=", best_gpr_hp$kernel,
        " nug_max=", best_gpr_hp$nug.max, "\n")
  } else {
    cat("  No saved GPR config found; using default GPR hyperparameters.\n")
  }
  gpr_final_fit <- fit_gpr(
    train_data = final_data,
    test_data = final_data,
    predictor_vars = gpr_pvars,
    hyperparams = best_gpr_hp
  )
  gpr_pred <- gpr_final_fit$predictions
  if (log_transform_target) gpr_pred <- inverse_response_transform(gpr_pred, log = TRUE)
  gpr_train_metrics <- calculate_metrics(core_data$median_carbon_density, gpr_pred)
  saveRDS(list(
    model              = gpr_final_fit$model,
    predictor_vars     = gpr_pvars,
    best_covariate_set = gpr_pvars,
    hyperparams        = best_gpr_hp,
    scale_params       = gpr_final_fit$scale_params,
    encoding           = gpr_final_fit$encoding,
    encoded_names      = gpr_final_fit$encoded_names,
    cv_metrics        = gpr_cv_metrics,
    train_metrics      = gpr_train_metrics,
    log_response       = log_transform_target
  ), file.path(out_dir, "GPR_final.rds"))
  cat("Saved GPR_final.rds  (train R2=", round(gpr_train_metrics$r2, 3),
      " RMSE=", round(gpr_train_metrics$rmse, 4), ")\n\n")
}
# ---------------------------------------------------------------------------
# LR baseline: no hyperparameter tuning (uses LR-specific covariates)
# ---------------------------------------------------------------------------
if ("LR" %in% model_list) {
  cat("=== LR: fitting linear baseline model ===\n")
  lm_pvars <- load_model_vars_final("LR")
  cat("  Predictors (", length(lm_pvars), "): ", paste(lm_pvars, collapse = ", "), "\n", sep = "")

  lm_cfg <- load_best_model_config(
    model_name = "lr",
    config_dir = config_dir,
    robust_config_dir = robust_tune_dir,
    prefer_robust = isTRUE(use_robust_final_configs),
    include_baseline = TRUE,
  )
  lm_cv_metrics <- if (!is.null(lm_cfg)) lm_cfg$cv_metrics else NULL

  lm_prep <- prepare_data_for_model("LR", train=final_data, test=final_data, predictor_vars=lm_pvars)
  lm_final <- fit_lm(lm_prep$train, lm_prep$test, lm_prep$predictor_vars)
  lm_pred <- lm_final$predictions
  if (log_transform_target) lm_pred <- inverse_response_transform(lm_pred, log = TRUE)
  lm_train_metrics <- calculate_metrics(core_data$median_carbon_density, lm_pred)
  saveRDS(list(
    model              = lm_final$model,
    predictor_vars     = lm_pvars,
    best_covariate_set = lm_pvars,
    hyperparams        = list(),
    cv_metrics         = lm_cv_metrics,
    train_metrics      = lm_train_metrics,
    scale_params       = lm_prep$scale_params,
    log_response       = log_transform_target
  ), file.path(out_dir, "LR_final.rds"))
  cat("Saved LR_final.rds  (train R2=", round(lm_train_metrics$r2, 3),
      " RMSE=", round(lm_train_metrics$rmse, 4), ")\n\n")
}
# ---------------------------------------------------------------------------
# Summary table (with variables included per model)
# ---------------------------------------------------------------------------

summary_rows <- list()
add_summary_row <- function(model, metrics, best_hp, vars) {
  summary_rows[[length(summary_rows) + 1L]] <<- data.frame(
    model = model,
    train_r2 = round(metrics$r2, 4),
    train_rmse = round(metrics$rmse, 4),
    best_hp = best_hp,
    variables = paste(vars, collapse = ", "),
    stringsAsFactors = FALSE
  )
}

if ("XGB" %in% model_list) {
  add_summary_row(
    "XGB",
    xgb_train_metrics,
    paste0(
      "nrounds=", best_xgb_hp$nrounds,
      " max_depth=", best_xgb_hp$max_depth,
      " learning_rate=", best_xgb_hp$learning_rate,
      " subsample=", best_xgb_hp$subsample,
      " colsample_bytree=", best_xgb_hp$colsample_bytree,
      " min_child_weight=", best_xgb_hp$min_child_weight,
      " min_split_loss=", best_xgb_hp$min_split_loss,
      " reg_reg_lambda=", best_xgb_hp$reg_reg_lambda
    ),
    xgb_pvars
  )
}

if ("GAM" %in% model_list) {
  add_summary_row("GAM", gam_train_metrics, paste0("k_covariate=", best_gam_hp$k_covariate), gam_pvars)
}

if ("GPR" %in% model_list) {
  add_summary_row(
    "GPR",
    gpr_train_metrics,
    paste0(
      "kernel=", best_gpr_hp$kernel,
      " nug_min=", best_gpr_hp$nug.min,
      " nug_max=", best_gpr_hp$nug.max,
      " nug_est=", best_gpr_hp$nug.est
    ),
    gpr_pvars
  )
}

if ("LR" %in% model_list) {
  add_summary_row("LR", lm_train_metrics, "none (linear baseline)", lm_pvars)
}

summary_tbl <- dplyr::bind_rows(summary_rows)

write.csv(summary_tbl, file.path(out_dir, "final_models_summary.csv"), row.names = FALSE)
cat("Summary saved to", file.path(out_dir, "final_models_summary.csv"), "\n")
print(summary_tbl)
cat("\nAll final models saved to", out_dir, "\n")
cat("To re-use: obj <- readRDS('", file.path(out_dir, "XGB_final.rds"), "')\n", sep = "")
cat("           predict(obj$model, xgb.DMatrix(as.matrix(new_data[, obj$predictor_vars])))\n\n")
