# =============================================================================
# ML models and data preparation for GPR, RF, XGB, GAM.
# Shared scaling/encoding; fit_*; predict_model (generic); fit_gpr (unified).
# =============================================================================
load_packages(c("GauPro", "xgboost", "mgcv", "randomForest"))

#' Null-coalescing operator: returns x if not NULL, otherwise y.
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

# ================================ SCALING & ENCODING ================================
# XGB/GPR: categoricals -> integer codes, then z-score scale (prepare_predictors_train / prepare_predictors_new).
# GAM/RF: scale numerics only, categoricals stay as factors (prepare_predictors_train_numeric_only).

#' Compute z-score scaling parameters from training data (numeric columns only).
#' @param data Training data frame.
#' @param predictor_vars Character vector of predictor column names to scale.
#' @return Named list with \code{means} and \code{sds} (named numeric vectors).
compute_scale_params <- function(data, predictor_vars) {
  pv_num <- predictor_vars[vapply(data[predictor_vars], is.numeric, logical(1))]
  if (length(pv_num) == 0L) return(list(means = numeric(), sds = numeric()))
  X <- data[, pv_num, drop = FALSE]
  means <- colMeans(X, na.rm = TRUE)
  sds   <- apply(X, 2, sd, na.rm = TRUE)
  sds[!is.finite(sds) | sds <= 0] <- 1
  list(means = means, sds = sds)
}

#' Apply z-score scaling using pre-computed parameters.
apply_scaling <- function(data, scale_params, predictor_vars = NULL) {
  data <- as.data.frame(data)
  vars_to_scale <- names(scale_params$means)
  if (length(vars_to_scale) == 0L) return(data)
  for (v in vars_to_scale) {
    if (!v %in% names(data)) next
    data[[v]] <- (data[[v]] - scale_params$means[v]) / scale_params$sds[v]
  }
  data
}

#' Log-transform the response variable in place (optional, default TRUE).
#' Use before fitting; pair with \code{inverse_response_transform} when predicting.
#' @param data Data frame containing the response column.
#' @param response_var Name of the response column (default \code{median_carbon_density}).
#' @param log If TRUE (default), apply log-transform; if FALSE, return data unchanged.
#' @param min_val Minimum value before log (avoids log(0)); default \code{.Machine$double.eps}.
#' @return Data frame with \code{data[[response_var]]} replaced by \code{log(data[[response_var]])} when \code{log = TRUE}.
transform_response <- function(data, response_var = "median_carbon_density", log = TRUE,
                               min_val = .Machine$double.eps) {
  if (!log || !response_var %in% names(data)) return(data)
  y <- data[[response_var]]
  data[[response_var]] <- log(pmax(as.numeric(y), min_val))
  data
}

#' Inverse of response log-transform for predictions (optional, default TRUE).
#' @param pred Numeric vector of predictions on log scale.
#' @param log If TRUE (default), apply exp; if FALSE, return pred unchanged.
#' @return Predictions on original scale when \code{log = TRUE}.
inverse_response_transform <- function(pred, log = TRUE) {
  if (!log) return(pred)
  exp(pred)
}

#' Convert categorical predictors to integer codes (one column each). No one-hot expansion.
#' @param data Data frame with predictor columns.
#' @param predictor_vars Character vector of predictor names.
#' @param categorical_vars Categorical columns to convert (default \code{seagrass_species}, \code{region}).
#' @return List with \code{data} (modified in place; categoricals as integer) and \code{encoding} (\code{levels} list for each categorical).
convert_categoricals_to_codes <- function(data, predictor_vars, categorical_vars = c("seagrass_species", "region")) {
  X <- as.data.frame(data[, predictor_vars, drop = FALSE])
  to_convert <- intersect(predictor_vars, categorical_vars)
  to_convert <- to_convert[vapply(X[to_convert], function(col) is.factor(col) || is.character(col), logical(1))]
  encoding <- list(levels = list())
  for (v in to_convert) {
    if (is.character(X[[v]])) X[[v]] <- factor(X[[v]])
    lv <- levels(factor(X[[v]]))
    encoding$levels[[v]] <- lv
    X[[v]] <- as.integer(factor(as.character(X[[v]]), levels = lv))
    X[[v]][is.na(X[[v]])] <- 1L
  }
  list(data = X, encoding = encoding)
}

#' Apply categorical encoding (factor to integer) to new data using stored levels.
apply_categorical_encoding <- function(data, encoding, predictor_vars) {
  if (is.null(encoding) || is.null(encoding$levels) || length(encoding$levels) == 0L)
    return(data[, predictor_vars, drop = FALSE])
  X <- as.data.frame(data[, predictor_vars, drop = FALSE])
  for (v in names(encoding$levels)) {
    if (!v %in% names(X)) next
    lev <- encoding$levels[[v]]
    X[[v]] <- as.integer(factor(as.character(X[[v]]), levels = lev))
    X[[v]][is.na(X[[v]])] <- 1L
  }
  X
}

#' Prepare predictor matrix from training data: categoricals to codes + scale. Single source of truth for all models.
#' @param data Training data (predictor columns only or full frame with predictor_vars).
#' @param predictor_vars Character vector of predictor names.
#' @param categorical_vars Categorical columns to convert to integer codes (default includes \code{Region} for compatibility).
#' @return List \code{data} (prepared df), \code{scale_params}, \code{encoding}.
prepare_predictors_train <- function(data, predictor_vars, categorical_vars = c("seagrass_species", "region", "Region")) {
  data <- as.data.frame(data[, predictor_vars, drop = FALSE])
  conv <- convert_categoricals_to_codes(data, predictor_vars, categorical_vars)
  sp   <- compute_scale_params(conv$data, predictor_vars)
  list(
    data          = apply_scaling(conv$data, sp, predictor_vars),
    scale_params  = sp,
    encoding      = conv$encoding
  )
}

#' Apply train-derived encoding and scaling to new data. Use with \code{prepare_predictors_train} for consistent prep.
prepare_predictors_new <- function(data, predictor_vars, scale_params, encoding) {
  X <- apply_categorical_encoding(data, encoding, predictor_vars)
  apply_scaling(X, scale_params, predictor_vars)
}

#' Prepare predictors for GAM/RF: scale numerics only, leave categoricals as factors.
prepare_predictors_train_numeric_only <- function(data, predictor_vars) {
  data <- as.data.frame(data[, predictor_vars, drop = FALSE])
  sp   <- compute_scale_params(data, predictor_vars)
  list(
    data         = apply_scaling(data, sp, predictor_vars),
    scale_params = sp,
    encoding     = NULL
  )
}

#' Prepare train/test data for a given model.
#' XGB/GPR: categoricals -> integer codes + scale. GAM/RF: scale numerics only, categoricals stay as factors.
#' GPR returns raw so \code{fit_gpr} can do prep internally (formula, prediction_grid, etc.).
prepare_data_for_model <- function(model_name, train, test, predictor_vars) {
  if (model_name == "GPR")
    return(list(train = train, test = test, predictor_vars = predictor_vars, scale_params = NULL, encoding = NULL))

  train <- as.data.frame(train)
  test  <- as.data.frame(test)
  if (model_name == "GAM" || model_name == "RF") {
    prep     <- prepare_predictors_train_numeric_only(train, predictor_vars)
    test_prep <- as.data.frame(apply_scaling(test[, predictor_vars, drop = FALSE], prep$scale_params, predictor_vars))
  } else {
    prep     <- prepare_predictors_train(train, predictor_vars)
    test_prep <- prepare_predictors_new(test, predictor_vars, prep$scale_params, prep$encoding)
  }
  # Attach response so fit_* have median_carbon_density; GAM also needs longitude, latitude for s(longitude, latitude)
  train_out <- prep$data
  if ("median_carbon_density" %in% names(train)) {
    train_out$median_carbon_density <- train$median_carbon_density
    test_prep$median_carbon_density <- test$median_carbon_density
  }
  if (model_name == "GAM" && all(c("longitude", "latitude") %in% names(train))) {
    train_out$longitude <- train$longitude
    train_out$latitude  <- train$latitude
    test_prep$longitude <- test$longitude
    test_prep$latitude  <- test$latitude
  }
  # Align factor levels in test to train so predict() does not fail on new levels (e.g. spatial CV)
  if (model_name == "GAM" || model_name == "RF") {
    for (col in predictor_vars) {
      if (col %in% names(train_out) && is.factor(train_out[[col]])) {
        tr_levels <- levels(train_out[[col]])
        test_prep[[col]] <- factor(as.character(test_prep[[col]]), levels = tr_levels)
        nas <- is.na(test_prep[[col]])
        if (any(nas)) test_prep[[col]][nas] <- tr_levels[1L]
      }
    }
  }
  list(
    train          = train_out,
    test           = test_prep,
    predictor_vars = predictor_vars,
    scale_params   = prep$scale_params,
    encoding       = prep$encoding
  )
}

# ================================ GENERIC PREDICTOR ================================

#' Normalize scale_params to shared format (means, sds). Handles legacy x_center/x_scale.
normalize_scale_params <- function(sp) {
  if (is.null(sp)) return(NULL)
  if (!is.null(sp$means)) return(sp)
  if (!is.null(sp$x_center)) return(list(means = sp$x_center, sds = sp$x_scale))
  NULL
}

#' Infer model type from saved object.
infer_model_type <- function(obj) {
  if (!is.null(obj$model_type)) return(obj$model_type)
  m <- obj$model
  if (inherits(m, "GauPro")) return("GPR")
  if (inherits(m, "xgb.Booster")) return("XGB")
  if (inherits(m, "gam")) return("GAM")
  if (inherits(m, "randomForest")) return("RF")
  stop("Could not infer model type from object.")
}

#' Generic predictor for saved models (GPR, XGB, GAM, RF).
#' @param obj Saved model list (model, predictor_vars, scale_params, encoding, encoded_names, log_response, ...)
#' @param newdata Data frame with predictor columns
#' @param se If TRUE and model supports it (GPR), return standard errors
#' @return List with \code{mean} (predictions on original scale) and \code{se} (NULL or on original scale if GPR)
predict_model <- function(obj, newdata, se = TRUE) {
  type <- infer_model_type(obj)
  pvars <- obj$predictor_vars
  X <- newdata[, pvars, drop = FALSE]
  # GAM formula uses s(longitude, latitude); ensure they are in X if present in newdata
  if (type == "GAM" && all(c("longitude", "latitude") %in% names(newdata))) {
    if (!"longitude" %in% names(X)) X$longitude <- newdata$longitude
    if (!"latitude" %in% names(X)) X$latitude <- newdata$latitude
  }
  if (!is.null(obj$encoding)) {
    X <- apply_categorical_encoding(X, obj$encoding, pvars)
  }
  model_cols <- obj$encoded_names %||% pvars

  # Scale if params present
  sp <- normalize_scale_params(obj$scale_params)
  if (!is.null(sp) && length(sp$means) > 0L) {
    X <- apply_scaling(X, sp, model_cols)
  }

  if (type == "GPR") {
    X_mat <- as.matrix(X[, model_cols, drop = FALSE])
    storage.mode(X_mat) <- "double"
    pred_out <- obj$model$pred(X_mat, se.fit = se, return_df = TRUE)
    if (is.list(pred_out) && !is.null(pred_out$mean)) {
      out <- list(mean = as.numeric(pred_out$mean), se = if (se) as.numeric(pred_out$se) else NULL)
    } else {
      out <- list(mean = as.numeric(pred_out), se = NULL)
    }
  } else if (type == "XGB") {
    X_mat <- as.matrix(X[, model_cols, drop = FALSE])
    storage.mode(X_mat) <- "double"
    out <- list(mean = as.numeric(predict(obj$model, newdata = X_mat)), se = NULL)
  } else if (type == "GAM") {
    out <- list(mean = as.numeric(mgcv::predict.gam(obj$model, newdata = X, type = "response")), se = NULL)
  } else if (type == "RF") {
    out <- list(mean = as.numeric(predict(obj$model, newdata = X)), se = NULL)
  } else {
    stop("Unsupported model type: ", type)
  }
  if (isTRUE(obj$log_response)) {
    out$mean <- inverse_response_transform(out$mean, log = TRUE)
    if (!is.null(out$se)) out$se <- out$se * out$mean
  }
  out
}

#' @rdname predict_model
#' @description \code{predict_gpr} is a backward-compatible alias for GPR objects.
predict_gpr <- function(gpr, newdata, se = TRUE) {
  predict_model(gpr, newdata, se = se)
}

# ================================ GPR ================================

#' Fit GPR on training data. Optionally predict on test_data (CV) or prediction_grid (spatial).
#'
#' @param train_data Training data.
#' @param test_data Optional test data for CV-style evaluation.
#' @param predictor_vars Predictor variable names.
#' @param value_var Response column (default: median_carbon_density).
#' @param prediction_grid Optional; if provided, predict on grid and return prediction_grid with gpr_mean, gpr_se.
#' @param coords Optional; for prediction_grid only; required if include_spatial=TRUE.
#' @param include_spatial If TRUE, add coords to predictor_vars (for spatial smoothing).
#' @param formula Optional formula; if NULL, value_var ~ all encoded predictors.
#' @param hyperparams kernel, nug.min, nug.max, nug.est.
#' @return Model list with model, scale_params (means/sds), encoding, encoded_names, predictor_vars, model_type.
#'   If test_data: predictions, r2, rmse. If prediction_grid: predictions, se, prediction_grid, n_train, n_pred.
fit_gpr <- function(train_data,
                    test_data = NULL,
                    predictor_vars,
                    value_var = "median_carbon_density",
                    prediction_grid = NULL,
                    coords = NULL, include_spatial = FALSE, formula = NULL,
                    hyperparams = NULL) {
  if (!requireNamespace("GauPro", quietly = TRUE))
    stop("GauPro package is required.")

  hp <- hyperparams %||% list()
  kernel  <- hp$kernel  %||% "matern52"
  nug_min <- hp$nug.min %||% 1e-8
  nug_max <- hp$nug.max %||% 100
  nug_est <- hp$nug.est %||% TRUE

  if (value_var != "median_carbon_density") {
    train_data$median_carbon_density <- train_data[[value_var]]
    if (!is.null(test_data)) test_data$median_carbon_density <- test_data[[value_var]]
  }

  if (include_spatial && length(coords) > 0L) {
    coords_to_add <- setdiff(coords, predictor_vars)
    if (length(coords_to_add) > 0) predictor_vars <- c(coords_to_add, predictor_vars)
  }

  train_data <- as.data.frame(train_data)
  train_X <- train_data[, predictor_vars, drop = FALSE]
  train_y <- train_data$median_carbon_density
  ok_train <- complete.cases(train_X) & !is.na(train_y)
  train_X <- train_X[ok_train, , drop = FALSE]
  train_y <- train_y[ok_train]
  if (nrow(train_X) < 2L) {
    out <- list(model = NULL, scale_params = NULL, encoding = NULL, encoded_names = predictor_vars,
                predictor_vars = predictor_vars, model_type = "GPR")
    if (!is.null(test_data)) out <- c(out, list(predictions = rep(NA_real_, nrow(test_data)), r2 = NA_real_, rmse = NA_real_))
    if (!is.null(prediction_grid)) out <- c(out, list(predictions = numeric(), se = numeric(), prediction_grid = prediction_grid, n_train = 0L, n_pred = 0L))
    return(out)
  }

  prep <- prepare_predictors_train(train_X, predictor_vars)
  train_sc <- prep$data
  sp <- prep$scale_params

  model_pvars <- if (!is.null(formula)) {
    intersect(setdiff(all.vars(formula), value_var), predictor_vars)
  } else predictor_vars
  if (length(model_pvars) == 0L) model_pvars <- predictor_vars

  form <- as.formula(paste(value_var, "~", paste(model_pvars, collapse = " + ")))
  train_df <- cbind(setNames(data.frame(train_y), value_var), train_sc)
  train_df <- train_df[, c(value_var, model_pvars), drop = FALSE]

  mdl <- tryCatch(
    GauPro::gpkm(form, data = train_df, kernel = kernel,
                 nug.min = nug_min, nug.max = nug_max, nug.est = nug_est),
    error = function(e) NULL
  )
  if (is.null(mdl)) {
    out <- list(model = NULL, scale_params = sp, encoding = prep$encoding, encoded_names = predictor_vars,
                predictor_vars = predictor_vars, model_type = "GPR")
    if (!is.null(test_data)) out <- c(out, list(predictions = rep(NA_real_, nrow(test_data)), r2 = NA_real_, rmse = NA_real_))
    if (!is.null(prediction_grid)) out <- c(out, list(predictions = numeric(), se = numeric(), prediction_grid = prediction_grid, n_train = nrow(train_data), n_pred = 0L))
    return(out)
  }

  base <- list(model = mdl, scale_params = sp, encoding = prep$encoding, encoded_names = predictor_vars,
               predictor_vars = predictor_vars, model_type = "GPR")

  if (!is.null(test_data)) {
    test_data <- as.data.frame(test_data)
    test_X <- test_data[, predictor_vars, drop = FALSE]
    ok_test <- complete.cases(test_X)
    na_pred <- rep(NA_real_, nrow(test_data))
    if (sum(ok_test) > 0L) {
      test_sc <- as.data.frame(prepare_predictors_new(test_X, predictor_vars, sp, prep$encoding))
      X_mat <- as.matrix(test_sc[ok_test, model_pvars, drop = FALSE])
      storage.mode(X_mat) <- "double"
      preds <- na_pred
      preds[ok_test] <- as.numeric(mdl$pred(X_mat, se.fit = FALSE))
    } else preds <- na_pred
    obs <- test_data$median_carbon_density
    valid <- !is.na(preds) & !is.na(obs) & is.finite(preds) & is.finite(obs)
    r2 <- if (sum(valid) >= 2) cor(preds[valid], obs[valid])^2 else NA_real_
    rmse <- if (sum(valid) >= 1) sqrt(mean((obs[valid] - preds[valid])^2)) else NA_real_
    return(c(base, list(predictions = preds, r2 = r2, rmse = rmse)))
  }

  if (!is.null(prediction_grid)) {
    missing_g <- setdiff(predictor_vars, names(prediction_grid))
    if (length(missing_g) > 0) stop("Missing from prediction grid: ", paste(missing_g, collapse = ", "))
    pred_ok <- complete.cases(prediction_grid[, predictor_vars, drop = FALSE])
    mu <- se <- rep(NA_real_, nrow(prediction_grid))
    if (sum(pred_ok) > 0L) {
      pred_sc <- as.data.frame(prepare_predictors_new(
        prediction_grid[pred_ok, predictor_vars, drop = FALSE], predictor_vars, sp, prep$encoding))
      XX <- as.matrix(pred_sc[, model_pvars, drop = FALSE])
      storage.mode(XX) <- "double"
      pr <- mdl$pred(XX, se.fit = TRUE, return_df = TRUE)
      mu[pred_ok] <- as.numeric(pr$mean)
      se[pred_ok] <- as.numeric(pr$se)
    }
    prediction_grid$gpr_mean <- mu
    prediction_grid$gpr_se   <- se
    return(c(base, list(predictions = mu[pred_ok], se = se[pred_ok], variance = se[pred_ok]^2,
                        prediction_grid = prediction_grid, n_train = nrow(train_data), n_pred = sum(pred_ok))))
  }

  base
}

#' @rdname fit_gpr
#' @description \code{fit_gaussian_process_regression} is a backward-compatible alias for the spatial (prediction_grid) case.
fit_gaussian_process_regression <- function(dat, value_var, coords, predictor_vars,
                                            formula = NULL, prediction_grid,
                                            include_spatial = FALSE,
                                            kernel = "matern52", nug_min = 1e-8, nug_max = 100, nug_est = TRUE) {
  fit_gpr(
    train_data = dat,
    test_data = NULL,
    predictor_vars = predictor_vars,
    value_var = value_var,
    prediction_grid = prediction_grid,
    coords = coords,
    include_spatial = include_spatial,
    formula = formula,
    hyperparams = list(kernel = kernel, nug.min = nug_min, nug.max = nug_max, nug.est = nug_est)
  )
}

# ================================ RF, XGB, GAM ================================

fit_rf <- function(train_data, test_data, predictor_vars, hyperparams = NULL) {
  formula_obj <- as.formula(paste("median_carbon_density ~", paste(predictor_vars, collapse = " + ")))
  ntree    <- hyperparams$ntree    %||% 500L
  mtry     <- hyperparams$mtry     %||% floor(sqrt(length(predictor_vars)))
  nodesize <- hyperparams$nodesize %||% NULL
  model <- randomForest(formula_obj, data = train_data,
    ntree = ntree, mtry = mtry, nodesize = nodesize, importance = TRUE)
  list(model = model, predictions = as.numeric(predict(model, newdata = test_data)))
}

fit_xgboost <- function(train_data, test_data, predictor_vars, hyperparams = NULL) {
  X_train <- as.matrix(train_data[, predictor_vars, drop = FALSE])
  X_test  <- as.matrix(test_data[, predictor_vars, drop = FALSE])
  storage.mode(X_train) <- "double"
  storage.mode(X_test)  <- "double"
  y_train <- as.numeric(train_data$median_carbon_density)

  nrounds          <- hyperparams$nrounds          %||% 100L
  max_depth        <- hyperparams$max_depth        %||% 6L
  learning_rate    <- hyperparams$learning_rate    %||% 0.3
  subsample        <- hyperparams$subsample        %||% 0.8
  colsample_bytree <- hyperparams$colsample_bytree %||% 0.8
  nthreads         <- hyperparams$nthreads         %||% 1L

  model <- xgboost(x = X_train, y = y_train, nrounds = nrounds,
    max_depth = max_depth, learning_rate = learning_rate,
    subsample = subsample, colsample_bytree = colsample_bytree,
    nthread = as.integer(nthreads), objective = "reg:squarederror")

  list(model = model, predictions = as.numeric(predict(model, newdata = X_test)))
}

fit_gam <- function(train_data, test_data, predictor_vars, k_spatial = 80,
                    include_spatial = TRUE) {
  if (!requireNamespace("mgcv", quietly = TRUE))
    return(list(model = NULL, predictions = rep(NA_real_, nrow(test_data))))

  y_train <- train_data$median_carbon_density
  family_used <- if (any(y_train <= 0, na.rm = TRUE)) gaussian() else Gamma(link = "log")
  covar_terms  <- paste(predictor_vars, collapse = " + ")

  if (include_spatial && all(c("longitude", "latitude") %in% names(train_data))) {
    spatial_term <- paste0("s(longitude, latitude, k = ", k_spatial, ")")
    form <- as.formula(paste("median_carbon_density ~", spatial_term, "+", covar_terms))
  } else {
    form <- as.formula(paste("median_carbon_density ~", covar_terms))
  }

  fit <- try(mgcv::gam(form, data = train_data, family = family_used, method = "REML"), silent = TRUE)
  if (inherits(fit, "try-error"))
    return(list(model = NULL, predictions = rep(NA_real_, nrow(test_data))))

  pred <- tryCatch(
    as.numeric(mgcv::predict.gam(fit, newdata = test_data, type = "response")),
    error = function(e) rep(mean(y_train, na.rm = TRUE), nrow(test_data))
  )
  list(model = fit, predictions = pred)
}
