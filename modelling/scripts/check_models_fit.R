#!/usr/bin/env Rscript
# =============================================================================
# Lightweight check that each model (GPR, GAM, XGB, RF) fits correctly
# with and without categorical variables. Run from project root.
#
# Usage: Rscript modelling/scripts/check_models_fit.R
#        Or: setwd(here::here()); source("modelling/scripts/check_models_fit.R")
# =============================================================================

setwd(here::here())
source("modelling/R/helpers.R")
source("modelling/R/assign_region_from_latlon.R")
# Ensure model packages available for fit checks
if (!requireNamespace("xgboost", quietly = TRUE)) install.packages("xgboost", quiet = TRUE)
if (!requireNamespace("randomForest", quietly = TRUE)) install.packages("randomForest", quiet = TRUE)
suppressPackageStartupMessages({
  library(xgboost)
  library(randomForest)
})

# -----------------------------------------------------------------------------
# Data: use real data if available, else minimal synthetic
# -----------------------------------------------------------------------------
target_var <- "median_carbon_density_100cm"
if (file.exists("data/all_extracted_new.rds")) {
  dat <- readr::read_rds("data/all_extracted_new.rds")
  dat <- process_rs_covariates(dat)
  if (!"region" %in% names(dat)) dat <- assign_region_from_latlon(dat)
  # Subset for speed and complete cases
  env_candidates <- setdiff(
    colnames(dat),
    c("latitude", "longitude", "number_id_final_version", "seagrass_species", "region", target_var)
  )
  env_candidates <- env_candidates[vapply(dat[env_candidates], is.numeric, logical(1))]
  env_vars <- head(env_candidates, 5L)  # use 5 env vars for speed
  need <- c("longitude", "latitude", target_var, env_vars)
  if ("seagrass_species" %in% names(dat)) need <- c(need, "seagrass_species")
  if ("region" %in% names(dat)) need <- c(need, "region")
  need <- intersect(need, colnames(dat))
  dat <- dat[complete.cases(dat[, need, drop = FALSE]), need, drop = FALSE]
  if (nrow(dat) > 400L) dat <- dat[sample(nrow(dat), 400L), , drop = FALSE]
  dat$median_carbon_density <- dat[[target_var]]
  env_vars <- intersect(env_vars, colnames(dat))
  has_species <- "seagrass_species" %in% colnames(dat)
  has_region  <- "region" %in% colnames(dat)
  cat("Using real data: n =", nrow(dat), ", env vars =", length(env_vars),
      ", species =", has_species, ", region =", has_region, "\n")
} else {
  # Minimal synthetic data
  set.seed(42)
  n <- 120
  dat <- data.frame(
    longitude = runif(n, -5, 10),
    latitude  = runif(n, 38, 55),
    median_carbon_density = exp(rnorm(n, 0, 0.5)),
    v1 = rnorm(n), v2 = rnorm(n), v3 = rnorm(n)
  )
  dat$median_carbon_density_100cm <- dat$median_carbon_density
  env_vars <- c("v1", "v2", "v3")
  has_species <- has_region <- FALSE
  cat("Using synthetic data: n =", nrow(dat), "\n")
}

stopifnot(length(env_vars) >= 2L, nrow(dat) >= 20L)

# -----------------------------------------------------------------------------
# Train/test split (single split for fit checks)
# -----------------------------------------------------------------------------
set.seed(1)
idx <- sample(nrow(dat), max(50L, floor(0.7 * nrow(dat))))
train <- dat[idx, , drop = FALSE]
test  <- dat[-idx, , drop = FALSE]
train$median_carbon_density <- train[[target_var]]
test$median_carbon_density  <- test[[target_var]]
train <- transform_response(train, "median_carbon_density", log = TRUE)
test  <- transform_response(test, "median_carbon_density", log = TRUE)

# -----------------------------------------------------------------------------
# Helpers: run one fit and return status
# -----------------------------------------------------------------------------
check_fit <- function(model_name, predictor_vars, label) {
  preds_len <- nrow(test)
  out <- list(ok = FALSE, msg = "", preds = NULL, obj = NULL)
  tryCatch({
    prep <- prepare_data_for_model(model_name, train, test, predictor_vars)
    fit <- switch(model_name,
      GPR = fit_gpr(prep$train, prep$predictor_vars, test_data = prep$test),
      XGB = fit_xgboost(prep$train, prep$test, prep$predictor_vars),
      GAM = fit_gam(prep$train, prep$test, prep$predictor_vars, k_spatial = 20L),
      RF  = fit_rf(prep$train, prep$test, prep$predictor_vars),
      stop("unknown model ", model_name)
    )
    preds <- fit$predictions
    if (is.null(fit$model)) {
      out$msg <- "model is NULL"
      return(out)
    }
    if (length(preds) != preds_len) {
      out$msg <- sprintf("predictions length %d != test rows %d", length(preds), preds_len)
      return(out)
    }
    n_na <- sum(is.na(preds))
    n_finite <- sum(is.finite(preds))
    if (n_finite == 0L) {
      out$msg <- "no finite predictions"
      return(out)
    }
    # Build object for predict_model round-trip (same shape as saved RDS).
    # GPR does its own prep inside fit_gpr, so use fit$scale_params and fit$encoding for GPR.
    scale_params <- if (model_name == "GPR" && !is.null(fit$scale_params)) fit$scale_params else prep$scale_params
    encoding     <- if (model_name == "GPR" && !is.null(fit$encoding)) fit$encoding else prep$encoding
    encoded_names <- if (model_name == "GPR" && !is.null(fit$encoded_names)) fit$encoded_names else prep$predictor_vars
    obj <- list(
      model = fit$model,
      predictor_vars = predictor_vars,
      scale_params = scale_params,
      encoding = encoding,
      encoded_names = encoded_names,
      model_type = model_name,
      log_response = TRUE
    )
    out$ok <- TRUE
    out$msg <- sprintf("n_finite=%d, n_na=%d", n_finite, n_na)
    out$preds <- preds
    out$obj <- obj
  }, error = function(e) {
    out$msg <- conditionMessage(e)
  })
  out
}

check_predict_model <- function(obj, newdata) {
  if (is.null(obj)) return(list(ok = FALSE, msg = "no obj"))
  tryCatch({
    res <- predict_model(obj, newdata = newdata, se = FALSE)
    preds <- res$mean
    if (length(preds) != nrow(newdata)) return(list(ok = FALSE, msg = "length mismatch"))
    if (all(!is.finite(preds))) return(list(ok = FALSE, msg = "no finite predictions"))
    list(ok = TRUE, msg = sprintf("n_finite=%d", sum(is.finite(preds))))
  }, error = function(e) list(ok = FALSE, msg = conditionMessage(e)))
}

# -----------------------------------------------------------------------------
# Run checks per model
# -----------------------------------------------------------------------------
models <- c("GPR", "GAM", "XGB", "RF")
results <- matrix("", nrow = length(models), ncol = 4L)
rownames(results) <- models
colnames(results) <- c("Env-only fit", "With categoricals fit", "predict_model round-trip", "Notes")

pv_env <- env_vars
pv_cat <- if (has_species || has_region) c(env_vars, "seagrass_species", "region") else env_vars
pv_cat <- intersect(pv_cat, colnames(dat))

for (m in models) {
  # Env-only
  r1 <- check_fit(m, pv_env, "env_only")
  results[m, 1L] <- if (r1$ok) paste0("OK (", r1$msg, ")") else paste0("FAIL: ", r1$msg)

  # With categoricals (if present)
  if (length(pv_cat) > length(pv_env)) {
    r2 <- check_fit(m, pv_cat, "with_cat")
    results[m, 2L] <- if (r2$ok) paste0("OK (", r2$msg, ")") else paste0("FAIL: ", r2$msg)
    obj_for_roundtrip <- r2$obj
  } else {
    results[m, 2L] <- "skip (no species/region in data)"
    obj_for_roundtrip <- r1$obj
  }

  # predict_model round-trip
  r3 <- check_predict_model(obj_for_roundtrip, test)
  results[m, 3L] <- if (r3$ok) paste0("OK (", r3$msg, ")") else paste0("FAIL: ", r3$msg)

  # Notes
  if (m == "GAM" && (!"longitude" %in% pv_env || !"latitude" %in% pv_env))
    results[m, 4L] <- "GAM uses s(lon,lat) from data; lon/lat added in prepare_data_for_model"
}

# -----------------------------------------------------------------------------
# Print summary
# -----------------------------------------------------------------------------
cat("\n")
cat(paste(rep("=", 72), collapse = ""), "\n")
cat("Model fit check summary\n")
cat(paste(rep("=", 72), collapse = ""), "\n\n")
print(as.data.frame(results))
cat("\n")

all_ok <- all(grepl("^OK", results[, 1L])) &&
  all(grepl("^OK|^skip", results[, 2L])) &&
  all(grepl("^OK", results[, 3L]))
if (all_ok) {
  cat("All checks passed.\n")
} else {
  cat("Some checks failed (see table above).\n")
}

# -----------------------------------------------------------------------------
# Optional: quick 2-fold CV for full pipeline (prepare_data_for_model + fit)
# -----------------------------------------------------------------------------
cat("\n--- Quick 2-fold CV (env-only) ---\n")
folds <- rep(1:2, length.out = nrow(dat))
for (m in models) {
  err <- tryCatch({
    rmse <- numeric(2)
    for (k in 1:2) {
      tr <- dat[folds != k, , drop = FALSE]
      te <- dat[folds == k, , drop = FALSE]
      tr$median_carbon_density <- tr[[target_var]]
      te$median_carbon_density <- te[[target_var]]
      tr <- transform_response(tr, "median_carbon_density", log = TRUE)
      te <- transform_response(te, "median_carbon_density", log = TRUE)
      prep <- prepare_data_for_model(m, tr, te, pv_env)
      fit <- switch(m,
        GPR = fit_gpr(prep$train, prep$predictor_vars, test_data = prep$test),
        XGB = fit_xgboost(prep$train, prep$test, prep$predictor_vars),
        GAM = fit_gam(prep$train, prep$test, prep$predictor_vars, k_spatial = 20L),
        RF  = fit_rf(prep$train, prep$test, prep$predictor_vars)
      )
      pred <- fit$predictions
      pred <- inverse_response_transform(pred, log = TRUE)
      rmse[k] <- sqrt(mean((te[[target_var]] - pred)^2, na.rm = TRUE))
    }
    cat("  ", m, ": mean RMSE =", round(mean(rmse), 5), "\n")
    NA_character_
  }, error = function(e) conditionMessage(e))
  if (!is.na(err)) cat("  ", m, ": ERROR -", err, "\n")
}

cat("\nDone.\n")
