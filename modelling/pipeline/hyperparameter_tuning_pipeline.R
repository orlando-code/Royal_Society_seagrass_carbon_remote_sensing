# Unified hyperparameter tuning for all models (GPR, GAM, XGB).
#
# Uses the best hyperparameters and best covariate set per model (from
# hyperparameter_tuning_pipeline and pruned_model_variables_shap/perm).
# Fits each model with best config on full data and computes mean |phi| per
# variable via iml::Shapley, then saves CSV and bar plot per model.
#
# Requires: hyperparameter_tuning_pipeline.R has been run (best_config_*.rds),
#           pruned_model_variables_shap.csv or pruned_model_variables_perm.csv
# Outputs:  output/cv_pipeline/importance_shap_<model>.csv and .png
#
# Usage: source after step 3 in run_paper.R, or run standalone.

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

core_data <- dat %>%
  dplyr::select(longitude, latitude, target_var, dplyr::all_of(predictor_vars),
                dplyr::all_of(intersect(c("seagrass_species", "region"), names(dat)))) %>%
  dplyr::filter(complete.cases(.))
core_data$median_carbon_density <- core_data[[target_var]]
predictor_vars <- predictor_vars[predictor_vars %in% colnames(core_data)]

if ("seagrass_species" %in% names(core_data)) {
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
  best_gpr <- tune_gpr(
    train_data = core_data,
    predictor_vars = pvars,
    n_folds = n_folds,
    verbose = TRUE,
    value_var = target_var,
    fold_indices = fold_indices,
    block_size = block_size,
    cache_tag = "tune_gpr",
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
  xgb_grid <- expand.grid(
    nrounds = c(100L, 300L, 500L),
    max_depth = c(3L, 4L, 6L),
    learning_rate = c(0.01, 0.05, 0.1),
    subsample = c(0.7, 0.8),
    colsample_bytree = c(0.6, 0.8),
    min_child_weight = c(1L, 5L, 10L),
    stringsAsFactors = FALSE
  )
  cv_tune_xgb <- function(core_data, predictor_vars, folds, hp, log_response = TRUE) {
    n_folds <- max(folds)
    metrics <- vector("list", n_folds)
    for (k in seq_len(n_folds)) {
      train <- core_data[folds != k, , drop = FALSE]
      test  <- core_data[folds == k, , drop = FALSE]
      observed_orig <- test$median_carbon_density
      if (log_response) {
        train <- transform_response(train, "median_carbon_density", log = TRUE)
        test  <- transform_response(test, "median_carbon_density", log = TRUE)
      }
      tryCatch({
        prep <- prepare_data_for_model("XGB", train, test, predictor_vars)
        res  <- fit_xgboost(prep$train, prep$test, prep$predictor_vars, hyperparams = hp)
        pred <- res$predictions
        if (log_response) pred <- inverse_response_transform(pred, log = TRUE)
        metrics[[k]] <- calculate_metrics(observed_orig, pred)
      }, error = function(e) NULL)
    }
    metrics <- dplyr::bind_rows(metrics)
    data.frame(r2 = mean(metrics$r2, na.rm = TRUE), rmse = mean(metrics$rmse, na.rm = TRUE))
  }
  xgb_results <- vector("list", nrow(xgb_grid))
  for (i in seq_len(nrow(xgb_grid))) {
    hp <- list(
      nrounds = xgb_grid$nrounds[i], max_depth = xgb_grid$max_depth[i],
      learning_rate = xgb_grid$learning_rate[i], subsample = xgb_grid$subsample[i],
      colsample_bytree = xgb_grid$colsample_bytree[i],
      min_child_weight = xgb_grid$min_child_weight[i]
    )
    m <- cv_tune_xgb(core_data, pvars, fold_indices, hp, log_response = log_response)
    xgb_results[[i]] <- cbind(xgb_grid[i, ], m)
    cat("  nrounds=", hp$nrounds, " depth=", hp$max_depth, " lr=", hp$learning_rate,
        " mcw=", hp$min_child_weight, " -> R2=", round(m$r2, 4), " RMSE=", round(m$rmse, 4), "\n")
  }
  xgb_cv_metrics <- dplyr::bind_rows(xgb_results)
  best_idx <- which.max(xgb_cv_metrics$r2)
  xgb_config <- list(
    nrounds = xgb_cv_metrics$nrounds[best_idx],
    max_depth = xgb_cv_metrics$max_depth[best_idx],
    learning_rate = xgb_cv_metrics$learning_rate[best_idx],
    subsample = xgb_cv_metrics$subsample[best_idx],
    colsample_bytree = xgb_cv_metrics$colsample_bytree[best_idx],
    min_child_weight = xgb_cv_metrics$min_child_weight[best_idx],
    cv_metrics = xgb_cv_metrics
  )
  saveRDS(xgb_config, file.path(out_dir, "best_config_xgb.rds"))
  cat("  Saved best_config_xgb.rds: nrounds =", xgb_config$nrounds,
      ", max_depth =", xgb_config$max_depth,
      ", lr =", xgb_config$learning_rate,
      ", mcw =", xgb_config$min_child_weight, "\n")
}

# ---------------------------------------------------------------------------
# GAM (k_covariate grid — basis dimension for smooth covariate terms)
# ---------------------------------------------------------------------------
if ("GAM" %in% model_list) {
  cat("\n=== GAM hyperparameter tuning (k_covariate, no spatial smooth) ===\n")
  pvars <- load_model_vars("GAM", per_model_vars, use_shap_first = use_shap_per_model)
  cat("  Predictor(s):", length(pvars), "\n")
  gam_k_grid <- c(2L, 4L, 6L, 8L, 10L, 12L, 15L, 20L)

  cv_tune_gam <- function(core_data, predictor_vars, folds, k_covariate, log_response = TRUE) {
    metrics <- vector("list", max(folds))
    for (k in seq_len(max(folds))) {
      train <- core_data[folds != k, , drop = FALSE]
      test  <- core_data[folds == k, , drop = FALSE]
      observed_orig <- test$median_carbon_density
      if (log_response) {
        train <- transform_response(train, "median_carbon_density", log = TRUE)
        test  <- transform_response(test, "median_carbon_density", log = TRUE)
      }
      metric_k <- tryCatch({
        prep <- prepare_data_for_model("GAM", train, test, predictor_vars)
        res  <- fit_gam(prep$train, prep$test, prep$predictor_vars,
                        include_spatial = FALSE, k_covariate = k_covariate)
        pred <- res$predictions
        if (log_response) pred <- inverse_response_transform(pred, log = TRUE)
        keep <- test_rows_with_factors_in_train(train, test, intersect(predictor_vars, c("seagrass_species", "region")))
        obs_sub <- observed_orig[keep]
        pred_sub <- pred[keep]
        if (sum(keep) < 2L) data.frame(r2 = NA_real_, rmse = NA_real_, mae = NA_real_, bias = NA_real_) else calculate_metrics(obs_sub, pred_sub)
      }, error = function(e) NULL)
      metrics[[k]] <- metric_k
    }
    metrics <- dplyr::bind_rows(metrics)
    data.frame(r2 = mean(metrics$r2, na.rm = TRUE), rmse = mean(metrics$rmse, na.rm = TRUE))
  }

  gam_results <- vector("list", length(gam_k_grid))
  for (i in seq_along(gam_k_grid)) {
    k_cov <- gam_k_grid[i]
    m     <- cv_tune_gam(core_data, pvars, fold_indices, k_cov, log_response = log_response)
    r2_val <- if (is.nan(m$r2) || is.infinite(m$r2)) NA_real_ else m$r2
    rmse_val <- if (is.nan(m$rmse) || is.infinite(m$rmse)) NA_real_ else m$rmse
    gam_results[[i]] <- data.frame(k_covariate = k_cov, r2 = r2_val, rmse = rmse_val)
    cat("  k_covariate =", k_cov,
        " -> R2 =", if (is.na(r2_val)) "NA" else round(r2_val, 4),
        ", RMSE =", if (is.na(rmse_val)) "NA" else round(rmse_val, 4), "\n")
  }
  gam_cv_metrics <- dplyr::bind_rows(gam_results)
  if (all(is.na(gam_cv_metrics$r2))) {
    best_idx <- which.min(gam_cv_metrics$rmse)
    warning("All R2 values are NA for GAM k_covariate grid search, selecting best k by minimum RMSE.")
  } else {
    best_idx <- which.max(gam_cv_metrics$r2)
  }
  best_gam_k     <- gam_cv_metrics$k_covariate[best_idx]
  gam_config     <- list(k_covariate = best_gam_k, cv_metrics = gam_cv_metrics)
  saveRDS(gam_config, file.path(out_dir, "best_config_gam.rds"))
  cat("  Saved best_config_gam.rds: k_covariate =", best_gam_k, "\n")
}

cat("\nHyperparameter tuning complete. Best configs in", out_dir, "\n")

