# =============================================================================
# Lightweight showcase: optimal species encoding per model (XGB, GPR, GAM, LR)
#
# Species handling (recommended / "Option A" where applicable):
#   GAM  - unordered factor(seagrass_species) + scaled env covariates
#   LR   - unordered factor(seagrass_species) + scaled env covariates
#   XGB  - one-hot species dummies + z-scaled numeric env (dummies not scaled)
#   GPR  - species fixed effect on training fold, GPR on env residuals only
#
# Usage (from repo root):
#   Rscript modelling/analysis/showcase_species_encoding.R
# =============================================================================

if (!exists("seagrass_init_repo", mode = "function", inherits = TRUE)) {
  init_path <- file.path("modelling", "R", "init_repo.R")
  if (!file.exists(init_path)) {
    ff <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
    if (!length(ff)) stop("Run from repo root.", call. = FALSE)
    script_path <- normalizePath(sub("^--file=", "", ff[[1]]), winslash = "/", mustWork = FALSE)
    init_path <- normalizePath(file.path(dirname(script_path), "..", "R", "init_repo.R"),
                               winslash = "/", mustWork = FALSE)
  }
  sys.source(init_path, envir = .GlobalEnv)
}

project_root <- seagrass_init_repo(
  packages = c("dplyr", "readr", "ggplot2"),
  source_files = "modelling/pipeline_config.R",
  include_helpers = TRUE,
  require_core_inputs = TRUE,
  check_renv = FALSE
)

cfg <- get_pipeline_config()
SPECIES_VAR <- "seagrass_species"
TARGET_VAR  <- cfg$target_var
CV_SEED     <- 42L
N_FOLDS     <- 5L
LOG_RESPONSE <- TRUE

# Small env set (stable core from sweep_001); keeps runtime low.
ENV_PREDICTORS <- c(
  "fe_mean_monthly_1.5m_mmol_m3",
  "vo_mean_1.5m_m_s",
  "surf_spco2_p95_uatm",
  "nppv_mean_monthly_0.5m_mg_m3_day",
  "sst_daily_p95_k",
  "bottomt_p95_daily_c"
)
ALL_PREDICTORS <- c(ENV_PREDICTORS, SPECIES_VAR)

# ---------------------------------------------------------------------------
# Encoding helpers (showcase-only; documents recommended handling)
# ---------------------------------------------------------------------------

describe_species_encoding <- function() {
  data.frame(
    model = c("GAM", "LR", "XGB", "GPR"),
    species_handling = c(
      "factor(seagrass_species) in formula",
      "factor(seagrass_species) in formula",
      "one-hot dummy columns (not z-scaled)",
      "training-fold species means subtracted before GPR; added back at predict"
    ),
    env_handling = c(
      "z-scaled numerics + optional smooths",
      "z-scaled numerics (linear terms)",
      "z-scaled numerics",
      "z-scaled numerics in GP kernel only"
    ),
    stringsAsFactors = FALSE
  )
}

compute_train_species_means <- function(train, response_var = TARGET_VAR,
                                        species_var = SPECIES_VAR) {
  means <- tapply(train[[response_var]], train[[species_var]], mean, na.rm = TRUE)
  structure(as.numeric(means), names = names(means))
}

lookup_species_means <- function(species_means, species_vec, fallback = 0) {
  out <- as.numeric(species_means[as.character(species_vec)])
  nas <- is.na(out)
  if (any(nas)) out[nas] <- fallback
  out
}

prepare_xgb_onehot_train <- function(data, predictor_vars,
                                     categorical_vars = c(SPECIES_VAR, "region", "Region")) {
  X <- as.data.frame(data[, predictor_vars, drop = FALSE])
  cat_vars <- intersect(predictor_vars, categorical_vars)
  cat_vars <- cat_vars[vapply(X[cat_vars], function(col) is.factor(col) || is.character(col), logical(1))]
  num_vars <- setdiff(predictor_vars, cat_vars)

  for (v in cat_vars) {
    if (is.character(X[[v]])) X[[v]] <- factor(X[[v]])
    X[[v]] <- factor(as.character(X[[v]]), levels = levels(factor(X[[v]])))
  }

  dummy_df <- if (length(cat_vars) > 0L) {
    mm <- model.matrix(as.formula(paste("~", paste(cat_vars, collapse = " + "), "- 1")), data = X)
    colnames(mm) <- make.names(colnames(mm), unique = TRUE)
    as.data.frame(mm)
  } else {
    data.frame()
  }

  num_df <- if (length(num_vars) > 0L) X[, num_vars, drop = FALSE] else data.frame()
  if (ncol(num_df) > 0L) {
    sp <- compute_scale_params(num_df, names(num_df))
    num_df <- apply_scaling(num_df, sp, names(num_df))
  } else {
    sp <- list(means = numeric(), sds = numeric())
  }

  out <- cbind(num_df, dummy_df)
  list(
    data = out,
    scale_params = sp,
    encoding = list(
      type = "xgb_one_hot",
      levels = stats::setNames(lapply(cat_vars, function(v) levels(X[[v]])), cat_vars),
      numeric_vars = num_vars,
      dummy_vars = names(dummy_df),
      predictor_vars = predictor_vars
    ),
    encoded_names = names(out)
  )
}

prepare_xgb_onehot_new <- function(data, encoding, encoded_names) {
  pvars <- encoding$predictor_vars
  X <- as.data.frame(data[, pvars, drop = FALSE])
  cat_vars <- intersect(names(encoding$levels), pvars)

  for (v in cat_vars) {
    X[[v]] <- factor(as.character(X[[v]]), levels = encoding$levels[[v]])
  }

  dummy_df <- if (length(cat_vars) > 0L) {
    mm <- model.matrix(as.formula(paste("~", paste(cat_vars, collapse = " + "), "- 1")), data = X)
    colnames(mm) <- make.names(colnames(mm), unique = TRUE)
    as.data.frame(mm)
  } else {
    data.frame()
  }

  num_vars <- encoding$numeric_vars
  num_df <- if (length(num_vars) > 0L) X[, num_vars, drop = FALSE] else data.frame()
  if (ncol(num_df) > 0L && length(encoding$scale_params$means) > 0L) {
    num_df <- apply_scaling(num_df, encoding$scale_params, num_vars)
  }

  out <- cbind(num_df, dummy_df)
  miss <- setdiff(encoded_names, names(out))
  for (m in miss) out[[m]] <- 0
  out[, encoded_names, drop = FALSE]
}

ensure_model_response <- function(df, target_var = TARGET_VAR) {
  df <- as.data.frame(df)
  if (!"median_carbon_density" %in% names(df) && target_var %in% names(df)) {
    df$median_carbon_density <- df[[target_var]]
  }
  df
}

ensure_species_factor <- function(df, species_var = SPECIES_VAR) {
  df <- as.data.frame(df)
  if (!species_var %in% names(df)) return(df)
  df[[species_var]] <- factor(as.character(df[[species_var]]))
  df
}

align_factor_levels <- function(train, test, factor_cols) {
  train <- as.data.frame(train)
  test  <- as.data.frame(test)
  for (col in factor_cols) {
    if (!col %in% names(train) || !is.factor(train[[col]])) next
    tr_levels <- levels(train[[col]])
    test[[col]] <- factor(as.character(test[[col]]), levels = tr_levels)
    nas <- is.na(test[[col]])
    if (any(nas)) test[[col]][nas] <- tr_levels[1L]
  }
  list(train = train, test = test)
}

rmse_vec <- function(obs, pred) {
  ok <- is.finite(obs) & is.finite(pred)
  if (!any(ok)) return(NA_real_)
  sqrt(mean((obs[ok] - pred[ok])^2))
}

r2_vec <- function(obs, pred) {
  ok <- is.finite(obs) & is.finite(pred)
  if (sum(ok) < 2L) return(NA_real_)
  ss_res <- sum((obs[ok] - pred[ok])^2)
  ss_tot <- sum((obs[ok] - mean(obs[ok]))^2)
  if (ss_tot <= 0) return(NA_real_)
  1 - ss_res / ss_tot
}

# GPR Option A: species fixed effect + env-only Gaussian process on residuals.
fit_gpr_species_adjusted <- function(train_raw, test_raw, env_predictors,
                                     hyperparams = NULL, log_response = LOG_RESPONSE) {
  train <- ensure_species_factor(train_raw)
  test  <- ensure_species_factor(test_raw)

  species_means <- compute_train_species_means(train, "median_carbon_density", SPECIES_VAR)
  global_mean   <- mean(train$median_carbon_density, na.rm = TRUE)

  train_adj <- train
  test_adj  <- test
  train_adj$median_carbon_density <- train$median_carbon_density -
    lookup_species_means(species_means, train[[SPECIES_VAR]], fallback = global_mean)

  fit <- fit_gpr(
    train_data = train_adj[, c(env_predictors, "median_carbon_density"), drop = FALSE],
    test_data = NULL,
    predictor_vars = env_predictors,
    value_var = "median_carbon_density",
    hyperparams = hyperparams %||% list(kernel = "matern52", nug.min = 1e-8, nug.max = 100, nug.est = TRUE)
  )

  if (is.null(fit$model)) {
    return(list(predictions = rep(NA_real_, nrow(test)), species_means = species_means, fit = fit))
  }

  test_X <- test[, env_predictors, drop = FALSE]
  ok <- complete.cases(test_X)
  preds <- rep(NA_real_, nrow(test))
  if (any(ok)) {
    test_sc <- prepare_predictors_new_onehot(test_X[ok, , drop = FALSE], fit$encoding, fit$encoded_names)
    test_sc <- as.data.frame(apply_scaling(test_sc, fit$scale_params, fit$encoded_names))
    X_mat <- as.matrix(test_sc[, fit$encoded_names, drop = FALSE])
    storage.mode(X_mat) <- "double"
    resid_pred <- as.numeric(fit$model$pred(X_mat, se.fit = FALSE))
    preds[ok] <- resid_pred +
      lookup_species_means(species_means, test[[SPECIES_VAR]][ok], fallback = global_mean)
  }

  if (log_response) preds <- inverse_response_transform(preds, log = TRUE)

  list(
    predictions = preds,
    species_means = species_means,
    fit = fit,
    encoding_note = "GPR env-only kernel + additive species means"
  )
}

fit_fold_model <- function(model, train_raw, test_raw, env_predictors, all_predictors,
                           log_response = LOG_RESPONSE) {
  keep <- test_rows_with_factors_in_train(train_raw, test_raw, SPECIES_VAR)
  if (!any(keep)) {
    return(list(predictions = numeric(), observed = numeric(), n_eval = 0L))
  }
  train_raw <- ensure_model_response(train_raw)
  test_raw  <- ensure_model_response(test_raw)
  train_raw <- train_raw[complete.cases(train_raw[, c(all_predictors, TARGET_VAR), drop = FALSE]), , drop = FALSE]
  test_raw  <- test_raw[keep, , drop = FALSE]
  observed_orig <- test_raw[[TARGET_VAR]]

  train <- train_raw
  test  <- test_raw
  if (log_response) {
    train <- transform_response(train, "median_carbon_density", log = TRUE)
    test  <- transform_response(test, "median_carbon_density", log = TRUE)
  }

  if (model == "GAM") {
    prep <- prepare_data_for_model("GAM", ensure_species_factor(train), ensure_species_factor(test), all_predictors)
    fit  <- fit_gam(prep$train, prep$test, all_predictors, k_covariate = 3L)
    pred_log <- fit$predictions
    preds <- if (log_response) inverse_response_transform(pred_log, log = TRUE) else pred_log
  } else if (model == "LR") {
    prep <- prepare_data_for_model("LR", ensure_species_factor(train), ensure_species_factor(test), all_predictors)
    fit  <- fit_lm(prep$train, prep$test, all_predictors)
    pred_log <- fit$predictions
    preds <- if (log_response) inverse_response_transform(pred_log, log = TRUE) else pred_log
  } else if (model == "XGB") {
    xgb_train <- ensure_species_factor(train)
    xgb_test  <- ensure_species_factor(test)
    prep <- prepare_xgb_onehot_train(xgb_train, all_predictors)
    test_mat <- prepare_xgb_onehot_new(xgb_test, c(prep$encoding, list(scale_params = prep$scale_params)), prep$encoded_names)
    train_mat <- prep$data
    train_mat$median_carbon_density <- xgb_train$median_carbon_density
    test_mat$median_carbon_density  <- xgb_test$median_carbon_density
    fit <- fit_xgboost(train_mat, test_mat, prep$encoded_names,
                       hyperparams = list(nrounds = 80L, max_depth = 4L, learning_rate = 0.1))
    pred_log <- fit$predictions
    preds <- if (log_response) inverse_response_transform(pred_log, log = TRUE) else pred_log
  } else if (model == "GPR") {
    fit <- fit_gpr_species_adjusted(train, test, env_predictors, log_response = log_response)
    preds <- fit$predictions
  } else {
    stop("Unknown model: ", model)
  }

  list(predictions = preds, observed = observed_orig, n_eval = length(observed_orig))
}

run_species_encoding_showcase <- function() {
  cat("\n=== Species encoding showcase ===\n\n")
  print(describe_species_encoding(), row.names = FALSE)

  dat <- readRDS(file.path(project_root, "data", "all_extracted_new.rds"))
  dat <- ensure_model_response(dat)
  dat <- dat[complete.cases(dat[, c(TARGET_VAR, ALL_PREDICTORS), drop = FALSE]), , drop = FALSE]
  dat[[SPECIES_VAR]] <- factor(as.character(dat[[SPECIES_VAR]]))

  # One-fold encoding preview (first training split, lightweight)
  pf_preview <- make_pixel_grouped_folds(dat, ALL_PREDICTORS, n_folds = N_FOLDS, seed = CV_SEED)
  tr_preview <- ensure_species_factor(ensure_model_response(dat[pf_preview$fold_indices != 1, , drop = FALSE]))
  cat("\n--- Encoding preview (fold 1 train) ---\n")
  cat("GAM/LR species column class:", class(tr_preview[[SPECIES_VAR]]), "\n")
  cat("GAM/LR species levels:", paste(levels(tr_preview[[SPECIES_VAR]]), collapse = ", "), "\n")
  xgb_prep <- prepare_xgb_onehot_train(tr_preview, ALL_PREDICTORS)
  cat("XGB encoded columns (", length(xgb_prep$encoded_names), "): ",
      paste(xgb_prep$encoded_names, collapse = ", "), "\n", sep = "")
  cat("GPR env-only kernel columns (species excluded): ",
      paste(ENV_PREDICTORS, collapse = ", "), "\n", sep = "")
  cat("GPR species handled via training-fold means (n = ",
      length(compute_train_species_means(tr_preview)), " levels)\n\n", sep = "")

  pf <- make_pixel_grouped_folds(dat, ALL_PREDICTORS, n_folds = N_FOLDS, seed = CV_SEED)
  fold_ids <- pf$fold_indices

  models <- c("GAM", "LR", "XGB", "GPR")
  rows <- vector("list", length(models) * N_FOLDS)
  pooled_store <- stats::setNames(vector("list", length(models)), models)
  for (m in models) pooled_store[[m]] <- list(obs = numeric(), pred = numeric())

  for (k in seq_len(N_FOLDS)) {
    train_raw <- dat[fold_ids != k, , drop = FALSE]
    test_raw  <- dat[fold_ids == k, , drop = FALSE]
    for (m in models) {
      res <- fit_fold_model(m, train_raw, test_raw, ENV_PREDICTORS, ALL_PREDICTORS)
      rows[[length(rows) + 1L]] <- data.frame(
        model = m,
        fold = k,
        rmse = rmse_vec(res$observed, res$predictions),
        r2 = r2_vec(res$observed, res$predictions),
        n_eval = res$n_eval,
        stringsAsFactors = FALSE
      )
      pooled_store[[m]]$obs <- c(pooled_store[[m]]$obs, res$observed)
      pooled_store[[m]]$pred <- c(pooled_store[[m]]$pred, res$predictions)
    }
  }

  by_fold <- dplyr::bind_rows(rows)
  pooled_df <- dplyr::bind_rows(lapply(models, function(m) {
    data.frame(
      model = m,
      pooled_rmse = rmse_vec(pooled_store[[m]]$obs, pooled_store[[m]]$pred),
      pooled_r2 = r2_vec(pooled_store[[m]]$obs, pooled_store[[m]]$pred),
      stringsAsFactors = FALSE
    )
  }))

  out_dir <- file.path(project_root, "output", "analysis", "species_encoding_showcase")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  encoding_path <- file.path(out_dir, "species_encoding_strategies.csv")
  fold_path     <- file.path(out_dir, "cv_by_fold.csv")
  summary_path  <- file.path(out_dir, "cv_summary.csv")

  write.csv(describe_species_encoding(), encoding_path, row.names = FALSE)
  write.csv(by_fold, fold_path, row.names = FALSE)
  write.csv(dplyr::left_join(
    by_fold |>
      dplyr::group_by(model) |>
      dplyr::summarise(mean_rmse = mean(rmse, na.rm = TRUE),
                       mean_r2 = mean(r2, na.rm = TRUE), .groups = "drop"),
    pooled_df,
    by = "model"
  ), summary_path, row.names = FALSE)

  cat("\nPredictors:", paste(ALL_PREDICTORS, collapse = ", "), "\n")
  cat("Folds:", N_FOLDS, " pixel-grouped; n =", nrow(dat), "\n\n")
  cat("Fold-wise mean metrics:\n")
  print(dplyr::left_join(
    by_fold |>
      dplyr::group_by(model) |>
      dplyr::summarise(mean_rmse = mean(rmse, na.rm = TRUE),
                       mean_r2 = mean(r2, na.rm = TRUE), .groups = "drop"),
    pooled_df,
    by = "model"
  ), row.names = FALSE)

  cat("\nSaved:\n  ", encoding_path, "\n  ", fold_path, "\n  ", summary_path, "\n", sep = "")

  invisible(list(by_fold = by_fold, pooled = pooled_df, encoding = describe_species_encoding()))
}

run_species_encoding_showcase()
