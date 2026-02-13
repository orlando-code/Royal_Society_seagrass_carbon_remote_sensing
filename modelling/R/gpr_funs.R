# =============================================================================
# Shared GPR utilities for fitting, partial dependence, and importance.
# No side effects or rm(list = ls()). Callers load data and options.
# Dependencies: GauPro; for plotting, ggplot2; for bind_rows, dplyr (or use base).
# Predictors are scaled (center + scale) before fitting for stable length-scales.
# Categorical predictors are one-hot encoded so GauPro receives numeric X only.
# =============================================================================

# -----------------------------------------------------------------------------
# One-hot encoding (GauPro requires numeric X)
# -----------------------------------------------------------------------------

#' One-hot encode only designated categorical predictors (e.g. seagrass_species, Region).
#' Other predictors are left as-is (must be numeric for GauPro).
#' @param X Data frame with predictor_vars (numeric and/or factor/character)
#' @param predictor_vars Character vector of column names
#' @param categorical_vars Character vector of names to one-hot encode (default: seagrass_species, Region)
#' @return List with X_numeric (data frame, all numeric), encoded_names (col names), encoding (for gpr_apply_onehot)
gpr_onehot_encode <- function(X, predictor_vars,
                              categorical_vars = c("seagrass_species", "Region")) {
  X <- X[, predictor_vars, drop = FALSE]
  # Only encode these specific variables when they are factor/character
  to_encode <- intersect(predictor_vars, categorical_vars)
  to_encode <- to_encode[vapply(X[to_encode, drop = FALSE], function(col) is.factor(col) || is.character(col), logical(1))]
  if (length(to_encode) == 0) {
    return(list(
      X_numeric = as.data.frame(X),
      encoded_names = predictor_vars,
      encoding = NULL
    ))
  }
  for (v in to_encode) {
    if (is.character(X[[v]])) X[[v]] <- factor(X[[v]])
  }
  mm <- model.matrix(~ . - 1, data = X)
  raw_names <- colnames(mm)
  # Sanitize names for formula (spaces/special chars break parsing)
  encoded_names <- gsub("[^A-Za-z0-9._]", "_", raw_names)
  encoded_names <- make.names(encoded_names, unique = TRUE)
  colnames(mm) <- encoded_names
  levels_list <- lapply(X[to_encode], function(x) levels(factor(x)))
  list(
    X_numeric = as.data.frame(mm),
    encoded_names = encoded_names,
    encoding = list(factor_vars = to_encode, levels = levels_list, encoded_names = encoded_names, raw_names = raw_names)
  )
}

#' Apply the same one-hot encoding to new data (test/grid).
#' @param X Data frame with same predictor_vars (numeric and factor/character)
#' @param encoding From gpr_onehot_encode; NULL = no categoricals, return X as-is
#' @param predictor_vars Original predictor var names (used when encoding is NULL)
#' @return Data frame with numeric columns only, same column names as encoding$encoded_names
gpr_apply_onehot <- function(X, encoding, predictor_vars) {
  if (is.null(encoding)) {
    return(X[, predictor_vars, drop = FALSE])
  }
  X <- X[, predictor_vars, drop = FALSE]
  for (v in encoding$factor_vars) {
    if (!v %in% names(X)) next
    lev <- encoding$levels[[v]]
    if (is.null(lev)) next
    x <- X[[v]]
    if (is.character(x)) x <- factor(x, levels = lev)
    else x <- factor(x, levels = lev)
    x[is.na(x)] <- lev[1]
    X[[v]] <- x
  }
  mm <- model.matrix(~ . - 1, data = X)
  raw_names <- if (!is.null(encoding$raw_names)) encoding$raw_names else encoding$encoded_names
  out <- matrix(0, nrow = nrow(X), ncol = length(encoding$encoded_names))
  colnames(out) <- encoding$encoded_names
  for (i in seq_along(encoding$encoded_names)) {
    nm_enc <- encoding$encoded_names[i]
    nm_raw <- raw_names[i]
    if (nm_raw %in% colnames(mm)) out[, nm_enc] <- mm[, nm_raw]
  }
  as.data.frame(out)
}

# -----------------------------------------------------------------------------
# Scaling helpers (train stats only; apply same transform to test/grid)
# -----------------------------------------------------------------------------

#' Compute scale parameters from training predictor matrix (numeric columns only).
#' @param X Data frame or matrix of predictors (training data only)
#' @param predictor_vars Character vector of column names to scale
#' @return List with x_center (named vector), x_scale (named vector; 1 where sd=0)
gpr_compute_scale_params <- function(X, predictor_vars) {
  x_center <- numeric(0)
  x_scale  <- numeric(0)
  for (v in predictor_vars) {
    if (!v %in% names(X)) next
    col <- X[[v]]
    if (!is.numeric(col)) next
    m <- mean(col, na.rm = TRUE)
    s <- sd(col, na.rm = TRUE)
    if (!is.finite(s) || s <= 0) s <- 1
    x_center[v] <- m
    x_scale[v]  <- s
  }
  list(x_center = x_center, x_scale = x_scale)
}

#' Apply scale parameters to a matrix or data frame (same columns as in params).
#' Only numeric columns in scale_params are scaled; factors/character are left unchanged.
#' @param X Matrix or data frame with predictor columns
#' @param scale_params List from gpr_compute_scale_params (x_center, x_scale)
#' @return Same type as X with numeric columns scaled; others unchanged
gpr_apply_scale <- function(X, scale_params) {
  if (is.null(scale_params) || is.null(scale_params$x_center)) return(X)
  is_df <- is.data.frame(X)
  if (is_df) {
    out <- as.data.frame(X)
    for (v in names(scale_params$x_center)) {
      if (!v %in% names(out)) next
      col <- out[[v]]
      if (!is.numeric(col)) next
      s <- scale_params$x_scale[v]
      if (!is.finite(s) || s <= 0) s <- 1
      out[[v]] <- (col - scale_params$x_center[v]) / s
    }
    return(out)
  }
  X <- as.matrix(X)
  out <- X
  for (v in names(scale_params$x_center)) {
    if (!v %in% colnames(out)) next
    col <- X[, v, drop = TRUE]
    if (!is.numeric(col)) next
    s <- scale_params$x_scale[v]
    if (!is.finite(s) || s <= 0) s <- 1
    out[, v] <- (col - scale_params$x_center[v]) / s
  }
  out
}

# -----------------------------------------------------------------------------
# Fit GPR for one CV fold (train/test)
# -----------------------------------------------------------------------------

#' Fit GPR on training data and predict on test data (for CV).
#'
#' @param train_data Data frame with response and predictors
#' @param test_data Data frame with same predictor columns
#' @param value_var Response variable name
#' @param predictor_vars Predictor variable names
#' @param kernel Kernel type ("matern52", "matern32", "gaussian")
#' @param nug.min Minimum nugget (passed to GauPro::gpkm; only if nug.est = TRUE)
#' @param nug.max Maximum nugget (passed to GauPro::gpkm)
#' @param nug.est Whether to estimate nugget (TRUE) or use fixed nug (FALSE with nug)
#' @return List with predictions, r2, rmse, and optional error message
fit_gpr_cv <- function(train_data, test_data, value_var, predictor_vars,
                       kernel = "matern52", nug.min = 1e-8, nug.max = 100, nug.est = TRUE) {
  if (!requireNamespace("GauPro", quietly = TRUE)) {
    stop("GauPro package is required. Install with: install.packages('GauPro')")
  }
  train_data <- as.data.frame(train_data)
  test_data <- as.data.frame(test_data)

  train_X <- train_data[, predictor_vars, drop = FALSE]
  train_y <- train_data[[value_var]]
  complete_train <- complete.cases(train_X) & !is.na(train_y)
  train_X <- train_X[complete_train, , drop = FALSE]
  train_y <- train_y[complete_train]

  if (nrow(train_X) < 2) {
    return(list(
      predictions = rep(NA, nrow(test_data)),
      r2 = NA, rmse = NA, error = "Insufficient training data"
    ))
  }

  test_X <- test_data[, predictor_vars, drop = FALSE]
  test_complete <- complete.cases(test_X)
  if (sum(test_complete) == 0) {
    return(list(
      predictions = rep(NA, nrow(test_data)),
      r2 = NA, rmse = NA, error = "No complete test cases"
    ))
  }

  # One-hot encode categoricals so GauPro receives numeric X only
  enc_train <- gpr_onehot_encode(train_X, predictor_vars)
  train_X_num <- enc_train$X_numeric
  encoded_names <- enc_train$encoded_names
  test_X_num <- gpr_apply_onehot(test_X, enc_train$encoding, predictor_vars)

  # Scale predictors from training stats (improves kernel length-scale optimization)
  scale_params <- gpr_compute_scale_params(train_X_num, encoded_names)
  train_X_scaled <- as.data.frame(gpr_apply_scale(train_X_num, scale_params))
  test_X_scaled <- as.data.frame(gpr_apply_scale(test_X_num, scale_params))

  formula_str <- paste(value_var, "~", paste(encoded_names, collapse = " + "))
  formula_obj <- as.formula(formula_str)
  train_combined <- cbind(train_y, train_X_scaled)
  names(train_combined)[1] <- value_var

  result <- tryCatch({
    gpr_model <- GauPro::gpkm(formula_obj, data = train_combined, kernel = kernel,
                              nug.min = nug.min, nug.max = nug.max, nug.est = nug.est)
    test_X_complete <- test_X_scaled[test_complete, , drop = FALSE]
    test_X_mat <- as.matrix(test_X_complete)
    storage.mode(test_X_mat) <- "double"
    test_predictions <- gpr_model$pred(test_X_mat, se.fit = FALSE)
    predictions <- rep(NA, nrow(test_data))
    predictions[test_complete] <- test_predictions

    test_y <- test_data[[value_var]][test_complete]
    valid <- !is.na(predictions[test_complete]) & !is.na(test_y) &
             is.finite(predictions[test_complete]) & is.finite(test_y)
    if (sum(valid) < 2) {
      return(list(
        predictions = predictions, r2 = NA, rmse = NA,
        error = "Insufficient valid predictions"
      ))
    }
    pred_valid <- predictions[test_complete][valid]
    obs_valid <- test_y[valid]
    r2 <- cor(pred_valid, obs_valid)^2
    rmse <- sqrt(mean((obs_valid - pred_valid)^2))
    list(predictions = predictions, r2 = r2, rmse = rmse, error = NULL)
  }, error = function(e) {
    list(
      predictions = rep(NA, nrow(test_data)),
      r2 = NA, rmse = NA, error = e$message
    )
  })
  result
}

# -----------------------------------------------------------------------------
# Fit full GPR for prediction on a grid (with optional kernel)
# -----------------------------------------------------------------------------

#' Fit GPR on full data and predict on a prediction grid.
#'
#' @param dat Training data (data frame)
#' @param value_var Response variable name
#' @param coords Coordinate column names (e.g. c("longitude", "latitude"))
#' @param predictor_vars Predictor variable names
#' @param formula Optional formula; if NULL, value_var ~ predictor_vars
#' @param prediction_grid Data frame with same predictor columns (and optionally coords)
#' @param include_spatial If TRUE, add coords to predictor_vars
#' @param kernel Kernel type ("matern52", "matern32", "gaussian"); default matern52
#' @return List with model, predictions, se, variance, prediction_grid, n_train, n_pred
fit_gaussian_process_regression <- function(dat, value_var, coords, predictor_vars,
                                           formula = NULL, prediction_grid,
                                           include_spatial = FALSE,
                                           kernel = "matern52") {
  if (!requireNamespace("GauPro", quietly = TRUE)) {
    stop("GauPro package is required. Install with: install.packages('GauPro')")
  }

  if (include_spatial) {
    coords_to_add <- setdiff(coords, predictor_vars)
    if (length(coords_to_add) > 0) {
      predictor_vars <- c(coords_to_add, predictor_vars)
    }
  }

  need_vars <- c(coords, value_var, predictor_vars)
  missing_in_dat <- setdiff(need_vars, names(dat))
  if (length(missing_in_dat) > 0) {
    stop("Missing from training data: ", paste(missing_in_dat, collapse = ", "))
  }

  train_data <- dat[, need_vars, drop = FALSE]
  train_data <- train_data[complete.cases(train_data), ]
  if (nrow(train_data) < 2) {
    stop("Insufficient complete cases for GPR fitting (need at least 2)")
  }

  # Restrict predictors to those actually present in the training data
  available_predictors <- intersect(predictor_vars, names(train_data))
  predictor_vars <- available_predictors
  if (length(predictor_vars) == 0) stop("No valid predictor variables remaining.")

  # Ensure the prediction grid has the same predictor columns
  missing_in_grid <- setdiff(predictor_vars, names(prediction_grid))
  if (length(missing_in_grid) > 0) {
    stop("Missing from prediction grid: ", paste(missing_in_grid, collapse = ", "))
  }

  # One-hot encode categorical predictors (e.g. seagrass_species) so GauPro gets numeric X only
  train_X <- train_data[, predictor_vars, drop = FALSE]
  enc <- gpr_onehot_encode(train_X, predictor_vars)
  train_X_num <- enc$X_numeric
  encoded_names <- enc$encoded_names
  encoding <- enc$encoding

  # Scale predictors from training data (improves kernel optimization and stability)
  scale_params <- gpr_compute_scale_params(train_X_num, encoded_names)
  train_X_scaled <- as.data.frame(gpr_apply_scale(train_X_num, scale_params))
  train_data_subset <- data.frame(
    y = train_data[[value_var]],
    train_X_scaled,
    stringsAsFactors = FALSE
  )
  names(train_data_subset)[1] <- value_var

  if (is.null(formula)) {
    formula_str <- paste(value_var, "~", paste(encoded_names, collapse = " + "))
  } else {
    formula_vars <- setdiff(all.vars(formula), value_var)
    formula_vars <- intersect(formula_vars, encoded_names)
    if (length(formula_vars) == 0) stop("No valid predictor variables in formula.")
    formula_str <- paste(value_var, "~", paste(formula_vars, collapse = " + "))
  }
  formula_obj <- as.formula(formula_str)
  formula_vars_needed <- all.vars(formula_obj)
  train_data_subset <- train_data_subset[, formula_vars_needed, drop = FALSE]

  # Keep only rows in the prediction grid with complete predictor data
  pred_complete <- complete.cases(prediction_grid[, predictor_vars, drop = FALSE])
  n_pred <- sum(pred_complete)

  # If no prediction rows have complete covariates, return NA predictions
  if (n_pred == 0L) {
    warning("No complete rows in prediction grid for the chosen predictor set; returning NA predictions.")
    gpr_model <- GauPro::gpkm(formula_obj, data = train_data_subset, kernel = kernel)
    mu <- rep(NA_real_, nrow(prediction_grid))
    se <- rep(NA_real_, nrow(prediction_grid))

    prediction_grid$gpr_mean <- mu
    prediction_grid$gpr_se   <- se

    return(list(
      model = gpr_model,
      mu = mu,
      se = se,
      var = se^2,
      prediction_grid = prediction_grid,
      n_train = nrow(train_data),
      n_pred = length(mu),
      scale_params = scale_params,
      encoding = encoding,
      encoded_names = encoded_names,
      predictor_vars = predictor_vars
    ))
  }

  # One-hot encode prediction grid with same encoding as training
  pred_X <- prediction_grid[pred_complete, predictor_vars, drop = FALSE]
  pred_X_num <- gpr_apply_onehot(pred_X, encoding, predictor_vars)
  XX_scaled <- gpr_apply_scale(pred_X_num, scale_params)

  gpr_model <- GauPro::gpkm(formula_obj, data = train_data_subset, kernel = kernel)
  # GauPro::pred expects a numeric matrix (same column order as model)
  XX_mat <- as.matrix(XX_scaled[, encoded_names, drop = FALSE])
  storage.mode(XX_mat) <- "double"
  pred_result <- gpr_model$pred(XX_mat, se.fit = TRUE, return_df = TRUE)
  mu <- as.numeric(pred_result$mean)
  se <- as.numeric(pred_result$se)

  prediction_grid$gpr_mean <- NA_real_
  prediction_grid$gpr_se   <- NA_real_
  prediction_grid$gpr_mean[pred_complete] <- mu
  prediction_grid$gpr_se[pred_complete]   <- se

  list(
    model           = gpr_model,
    predictions     = mu,
    se              = se,
    variance        = se^2,
    prediction_grid = prediction_grid,
    n_train         = nrow(train_data),
    n_pred          = sum(pred_complete),
    scale_params    = scale_params,
    encoding        = encoding,
    encoded_names   = encoded_names,
    predictor_vars  = predictor_vars
  )
}

# -----------------------------------------------------------------------------
# Partial dependence (training-data ranges; optional range_quantiles)
# -----------------------------------------------------------------------------

#' Create partial dependence data for one variable over training data range.
#'
#' All evaluation is over the training data variable range (optionally restricted
#' by range_quantiles). Use create_partial_dependence + plot_partial_dependence
#' for paper-ready PDPs.
#'
#' @param gpr_model Fitted GauPro model
#' @param dat Training data used to fit the model
#' @param var_name Name of the variable for partial dependence
#' @param predictor_vars Vector of all predictor variable names
#' @param n_points Number of points to evaluate (default 50)
#' @param reference_method "median", "mean", "sample", or "conditional"
#' @param n_samples If reference_method = "sample", number of samples
#' @param use_conditional If TRUE, use conditional partial dependence (actual data bins)
#' @param range_quantiles Optional c(low, high) e.g. c(0.01, 0.99) to restrict to percentile range
#' @param scale_params Optional list from fit (x_center, x_scale); if provided, data is scaled before prediction
#' @return Data frame with var_value, prediction_mean, prediction_se (and optionally n_points)
create_partial_dependence <- function(gpr_model, dat, var_name, predictor_vars,
                                      n_points = 50, reference_method = "median",
                                      n_samples = 100, use_conditional = FALSE,
                                      range_quantiles = NULL, scale_params = NULL) {
  if (!var_name %in% predictor_vars) {
    stop("Variable '", var_name, "' not found in predictor_vars")
  }
  missing_in_dat <- setdiff(predictor_vars, names(dat))
  if (length(missing_in_dat) > 0) {
    stop("Missing predictor variables in data: ", paste(missing_in_dat, collapse = ", "))
  }

  var_values <- dat[[var_name]]
  if (all(is.na(var_values))) {
    stop("Variable '", var_name, "' contains only NA values.")
  }
  var_range <- range(var_values, na.rm = TRUE)
  if (!is.null(range_quantiles) && length(range_quantiles) >= 2) {
    var_range <- quantile(var_values, range_quantiles, na.rm = TRUE)
  }
  if (diff(var_range) == 0) {
    warning("Variable '", var_name, "' has no variation in range.")
  }

  if (use_conditional || reference_method == "conditional") {
    # For conditional PDP, restrict the sequence of evaluation points to the
    # central quantile range if requested. This trims off extreme values
    # before constructing the conditional bins.
    if (!is.null(range_quantiles) && length(range_quantiles) >= 2) {
      q <- quantile(var_values, range_quantiles, na.rm = TRUE)
      var_values_use <- var_values[var_values >= q[1] & var_values <= q[2]]
      if (length(var_values_use) == 0L) {
        stop("No values of '", var_name, "' remain after applying range_quantiles.")
      }
    } else {
      var_values_use <- var_values
    }

    quantile_probs <- seq(0, 1, length.out = n_points + 1)
    quantile_probs <- quantile_probs[-length(quantile_probs)]
    var_seq <- quantile(var_values_use, quantile_probs, na.rm = TRUE)
    var_seq <- unique(var_seq)
    n_points <- length(var_seq)
  } else {
    var_seq <- seq(var_range[1], var_range[2], length.out = n_points)
  }

  if (use_conditional || reference_method == "conditional") {
    pred_list <- lapply(seq_along(var_seq), function(i) {
      v_val <- var_seq[i]
      if (i < length(var_seq)) bin_width <- var_seq[i + 1] - v_val
      else bin_width <- v_val - var_seq[i - 1]
      tolerance <- bin_width * 0.5
      in_bin <- abs(dat[[var_name]] - v_val) <= tolerance
      if (sum(in_bin, na.rm = TRUE) == 0) {
        distances <- abs(dat[[var_name]] - v_val)
        in_bin <- distances <= quantile(distances, 0.1, na.rm = TRUE)
      }
      bin_data <- dat[in_bin & !is.na(in_bin), predictor_vars, drop = FALSE]
      if (nrow(bin_data) == 0) return(NULL)
      pred_input <- as.matrix(bin_data)
      if (!is.null(scale_params)) pred_input <- gpr_apply_scale(pred_input, scale_params)
      pred_result <- tryCatch({
        gpr_model$pred(pred_input, se.fit = TRUE, return_df = TRUE)
      }, error = function(e) NULL)
      if (!is.null(pred_result) && nrow(pred_result) > 0) {
        data.frame(
          var_value = v_val,
          prediction_mean = mean(pred_result$mean, na.rm = TRUE),
          prediction_se = sqrt(mean(pred_result$se^2, na.rm = TRUE)),
          n_points = nrow(bin_data)
        )
      } else NULL
    })
  } else {
    if (reference_method == "sample") {
      sample_rows <- sample(nrow(dat), min(n_samples, nrow(dat)))
      reference_data <- dat[sample_rows, predictor_vars, drop = FALSE]
    } else if (!reference_method %in% c("median", "mean")) {
      stop("reference_method must be 'median', 'mean', 'sample', or 'conditional'")
    }

    if (reference_method == "sample") {
      pred_list <- lapply(var_seq, function(v_val) {
        ref_data_mod <- reference_data[, predictor_vars, drop = FALSE]
        ref_data_mod[[var_name]] <- v_val
        pred_input <- as.matrix(ref_data_mod)
        if (!is.null(scale_params)) pred_input <- gpr_apply_scale(pred_input, scale_params)
        pred_result <- tryCatch({
          gpr_model$pred(pred_input, se.fit = TRUE, return_df = TRUE)
        }, error = function(e) NULL)
        if (!is.null(pred_result) && nrow(pred_result) > 0) {
          data.frame(var_value = v_val,
                    prediction_mean = mean(pred_result$mean, na.rm = TRUE),
                    prediction_se = sqrt(mean(pred_result$se^2, na.rm = TRUE)))
        } else NULL
      })
    } else {
      n_ref_points <- min(50, nrow(dat))
      sample_indices <- sample(nrow(dat), n_ref_points)
      reference_base <- dat[sample_indices, predictor_vars, drop = FALSE]
      pred_list <- lapply(var_seq, function(v_val) {
        pred_data <- reference_base
        pred_data[[var_name]] <- v_val
        pred_input <- as.matrix(pred_data)
        if (!is.null(scale_params)) pred_input <- gpr_apply_scale(pred_input, scale_params)
        pred_result <- tryCatch({
          gpr_model$pred(pred_input, se.fit = TRUE, return_df = TRUE)
        }, error = function(e) NULL)
        if (!is.null(pred_result) && nrow(pred_result) > 0) {
          data.frame(var_value = v_val,
                    prediction_mean = mean(pred_result$mean, na.rm = TRUE),
                    prediction_se = sqrt(mean(pred_result$se^2, na.rm = TRUE)))
        } else NULL
      })
    }
  }

  pred_list <- pred_list[!sapply(pred_list, is.null)]
  if (length(pred_list) == 0) {
    stop("No successful predictions for partial dependence.")
  }
  result_df <- do.call(rbind, pred_list)
  if (diff(range(result_df$prediction_mean, na.rm = TRUE)) < 1e-6) {
    warning("Partial dependence is nearly constant; variable may have little effect.")
  }
  result_df
}

#' Plot partial dependence for one or more variables (training data ranges).
#'
#' @param gpr_model Fitted GauPro model
#' @param dat Training data
#' @param var_names Variable name(s) to plot
#' @param predictor_vars All predictor variable names
#' @param n_points Evaluation points per variable
#' @param reference_method Reference method for other variables
#' @param n_samples If reference_method = "sample", number of samples
#' @param use_conditional If TRUE, conditional partial dependence
#' @param target_var_name Label for y-axis
#' @param range_quantiles Optional c(low, high) to restrict to percentile range
#' @param show_training_range If TRUE, add rug of training values (requires ggplot2)
#' @param scale_params Optional scale params from fit (for scaled models)
#' @return ggplot object (requires ggplot2)
plot_partial_dependence <- function(gpr_model, dat, var_names, predictor_vars,
                                     n_points = 50, reference_method = "median",
                                     n_samples = 100, use_conditional = FALSE,
                                     target_var_name = "Response",
                                     range_quantiles = NULL,
                                     show_training_range = FALSE,
                                     scale_params = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plot_partial_dependence.")
  }

  if (length(var_names) == 1) {
    pd_data <- create_partial_dependence(
      gpr_model, dat, var_names[1], predictor_vars,
      n_points = n_points, reference_method = reference_method,
      n_samples = n_samples, use_conditional = use_conditional,
      range_quantiles = range_quantiles, scale_params = scale_params
    )
    # Y-range of PDP ribbon
    y_min <- min(pd_data$prediction_mean - 1.96 * pd_data$prediction_se, na.rm = TRUE)
    y_max <- max(pd_data$prediction_mean + 1.96 * pd_data$prediction_se, na.rm = TRUE)

    # Background histogram of training values, scaled into lower band of y-range
    hist_layer <- NULL
    if (var_names[1] %in% names(dat)) {
      vals <- dat[[var_names[1]]]
      vals <- vals[is.finite(vals)]
      if (length(vals) > 0) {
        h <- graphics::hist(vals, plot = FALSE, breaks = "FD")
        if (length(h$mids) > 0 && max(h$counts) > 0) {
          h_df <- data.frame(
            x = h$mids,
            count = h$counts
          )
          h_df$density <- h_df$count / max(h_df$count)
          band_frac <- 0.25
          h_df$y <- y_min + h_df$density * (y_max - y_min) * band_frac
          hist_layer <- ggplot2::geom_col(
            data = h_df,
            ggplot2::aes(x = .data$x, y = .data$y),
            inherit.aes = FALSE,
            fill = "grey80",
            alpha = 0.6,
            width = diff(h$breaks)[1]
          )
        }
      }
    }

    p <- ggplot2::ggplot(pd_data, ggplot2::aes(x = .data$var_value, y = .data$prediction_mean)) +
      hist_layer +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data$prediction_mean - 1.96 * .data$prediction_se,
                    ymax = .data$prediction_mean + 1.96 * .data$prediction_se),
        alpha = 0.2, fill = "steelblue"
      ) +
      ggplot2::geom_line(color = "steelblue", linewidth = 1.2) +
      # Heavily smoothed, *linear* loess to avoid ridiculous wiggles
      ggplot2::geom_smooth(
        se = FALSE,
        color = "black",
        linewidth = 0.9,
        method = "loess",
        span = 1.5,
        formula = y ~ x,
        degree = 1
      ) +
      ggplot2::labs(
        x = var_names[1],
        y = paste("Partial effect on", target_var_name),
        title = paste("Partial Dependence:", var_names[1])
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 12, face = "bold"))
    if (show_training_range && var_names[1] %in% names(dat)) {
      rug_dat <- data.frame(x = dat[[var_names[1]]])
      p <- p + ggplot2::geom_rug(
        data = rug_dat, ggplot2::aes(x = .data$x),
        y = NULL, inherit.aes = FALSE, sides = "b", alpha = 0.3
      )
    }
  } else {
    pd_list <- lapply(var_names, function(v) {
      pd <- create_partial_dependence(
        gpr_model, dat, v, predictor_vars,
        n_points = n_points, reference_method = reference_method,
        n_samples = n_samples, use_conditional = use_conditional,
        range_quantiles = range_quantiles, scale_params = scale_params
      )
      pd$variable <- v
      pd
    })
    pd_data <- do.call(rbind, pd_list)

    # Precompute per-row lower/upper PDP bounds for later scaling of histograms
    pd_data$y_lower <- pd_data$prediction_mean - 1.96 * pd_data$prediction_se
    pd_data$y_upper <- pd_data$prediction_mean + 1.96 * pd_data$prediction_se

    # Build per-variable histograms and scale into lower y-band
    hist_dfs <- lapply(var_names, function(v) {
      if (!v %in% names(dat)) return(NULL)
      vals <- dat[[v]]
      vals <- vals[is.finite(vals)]
      if (length(vals) == 0) return(NULL)
      h <- graphics::hist(vals, plot = FALSE, breaks = "FD")
      if (length(h$mids) == 0 || max(h$counts) == 0) return(NULL)
      data.frame(
        variable = v,
        x = h$mids,
        count = h$counts,
        width = diff(h$breaks)[1]
      )
    })
    hist_dfs <- hist_dfs[!vapply(hist_dfs, is.null, logical(1))]
    hist_df <- if (length(hist_dfs) > 0) do.call(rbind, hist_dfs) else NULL

    if (!is.null(hist_df)) {
      # Per-variable min/max of PDP ribbon (mean ± 1.96 * se)
      yr_min <- aggregate(y_lower ~ variable, data = pd_data, FUN = min, na.rm = TRUE)
      names(yr_min)[2] <- "y_min"
      yr_max <- aggregate(y_upper ~ variable, data = pd_data, FUN = max, na.rm = TRUE)
      names(yr_max)[2] <- "y_max"
      yr <- merge(yr_min, yr_max, by = "variable", all.x = TRUE)

      hist_df <- merge(hist_df, yr, by = "variable", all.x = TRUE)
      hist_df$density <- with(hist_df, ifelse(count > 0, count / ave(count, variable, FUN = max), 0))
      band_frac <- 0.25
      hist_df$y <- with(hist_df, y_min + density * (y_max - y_min) * band_frac)
    }

    p <- ggplot2::ggplot(pd_data, ggplot2::aes(x = .data$var_value, y = .data$prediction_mean)) +
      (if (!is.null(hist_df)) ggplot2::geom_col(
        data = hist_df,
        ggplot2::aes(x = .data$x, y = .data$y),
        inherit.aes = FALSE,
        fill = "grey80",
        alpha = 0.6,
        width = hist_df$width
      ) else NULL) +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data$prediction_mean - 1.96 * .data$prediction_se,
                    ymax = .data$prediction_mean + 1.96 * .data$prediction_se),
        alpha = 0.2, fill = "steelblue"
      ) +
      ggplot2::geom_line(color = "steelblue", linewidth = 1) +
      ggplot2::geom_smooth(se = FALSE, color = "black", linewidth = 0.8, method = "loess", span = 1.5) +
      ggplot2::facet_wrap(~ variable, scales = "free_x", ncol = 2) +
      ggplot2::labs(
        x = "Variable value",
        y = paste("Partial effect on", target_var_name),
        title = "Partial Dependence Plots (training data ranges)"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 14, face = "bold"),
        strip.text = ggplot2::element_text(size = 10, face = "bold")
      )
    if (show_training_range) {
      rug_dfs <- lapply(var_names, function(v) {
        if (v %in% names(dat)) {
          data.frame(variable = v, x = dat[[v]])
        } else NULL
      })
      rug_dfs <- rug_dfs[!sapply(rug_dfs, is.null)]
      if (length(rug_dfs) > 0) {
        rug_df <- do.call(rbind, rug_dfs)
        p <- p + ggplot2::geom_rug(
          data = rug_df, ggplot2::aes(x = .data$x),
          y = NULL, inherit.aes = FALSE, sides = "b", alpha = 0.3
        )
      }
    }
  }
  p
}

# -----------------------------------------------------------------------------
# Permutation importance
# -----------------------------------------------------------------------------

#' Compute permutation importance for a fitted GPR model.
#'
#' @param gpr_model Fitted GauPro model
#' @param data Data frame with value_var and predictor_vars (complete cases used)
#' @param value_var Response variable name
#' @param predictor_vars Predictor variable names
#' @param n_permutations Number of permutations per variable (default 5)
#' @param n_val Max number of rows to use for validation (default 500)
#' @param scale_params Optional scale params from fit (for scaled models)
#' @return Data frame with variable, rmse_increase, mae_increase, rmse_pct_increase
gpr_permutation_importance <- function(gpr_model, data, value_var, predictor_vars,
                                      n_permutations = 5, n_val = 500, scale_params = NULL) {
  train_data <- data[, c(value_var, predictor_vars), drop = FALSE]
  train_data <- train_data[complete.cases(train_data), ]
  n_val <- min(n_val, nrow(train_data))
  val_indices <- sample(nrow(train_data), n_val)
  val_data <- train_data[val_indices, ]

  val_X <- as.matrix(val_data[, predictor_vars, drop = FALSE])
  if (!is.null(scale_params)) val_X <- gpr_apply_scale(val_X, scale_params)
  baseline_pred <- gpr_model$pred(val_X, se.fit = FALSE)
  baseline_rmse <- sqrt(mean((val_data[[value_var]] - baseline_pred)^2, na.rm = TRUE))
  baseline_mae <- mean(abs(val_data[[value_var]] - baseline_pred), na.rm = TRUE)

  numeric_vars <- predictor_vars[!predictor_vars %in% c("seagrass_species", "Region")]
  importance_results <- data.frame(
    variable = character(),
    rmse_increase = numeric(),
    mae_increase = numeric(),
    rmse_pct_increase = numeric(),
    stringsAsFactors = FALSE
  )

  for (var in numeric_vars) {
    rmse_increases <- numeric(n_permutations)
    mae_increases <- numeric(n_permutations)
    for (perm in seq_len(n_permutations)) {
      val_perm <- val_data
      val_perm[[var]] <- sample(val_perm[[var]])
      val_perm_X <- as.matrix(val_perm[, predictor_vars, drop = FALSE])
      if (!is.null(scale_params)) val_perm_X <- gpr_apply_scale(val_perm_X, scale_params)
      pred_perm <- tryCatch({
        gpr_model$pred(val_perm_X, se.fit = FALSE)
      }, error = function(e) rep(NA, nrow(val_perm)))
      if (!all(is.na(pred_perm))) {
        rmse_increases[perm] <- sqrt(mean((val_data[[value_var]] - pred_perm)^2, na.rm = TRUE)) - baseline_rmse
        mae_increases[perm] <- mean(abs(val_data[[value_var]] - pred_perm), na.rm = TRUE) - baseline_mae
      }
    }
    avg_rmse <- mean(rmse_increases, na.rm = TRUE)
    avg_mae <- mean(mae_increases, na.rm = TRUE)
    importance_results <- rbind(importance_results, data.frame(
      variable = var,
      rmse_increase = avg_rmse,
      mae_increase = avg_mae,
      rmse_pct_increase = 100 * avg_rmse / baseline_rmse,
      stringsAsFactors = FALSE
    ))
  }
  importance_results[order(importance_results$rmse_increase, decreasing = TRUE), ]
}

# -----------------------------------------------------------------------------
# Optional: interval coverage (calibration) on held-out predictions
# -----------------------------------------------------------------------------

#' Compute coverage of prediction intervals (e.g. 90% or 95%) on test predictions.
#'
#' @param observed Numeric vector of observed values
#' @param pred_mean Numeric vector of predicted means
#' @param pred_se Numeric vector of prediction standard errors
#' @param level Coverage level (e.g. 0.9 or 0.95)
#' @return List with coverage proportion, level, and interval bounds used
gpr_interval_coverage <- function(observed, pred_mean, pred_se, level = 0.95) {
  valid <- !is.na(observed) & !is.na(pred_mean) & !is.na(pred_se) & is.finite(pred_se) & pred_se > 0
  if (sum(valid) == 0) return(list(coverage = NA, level = level, n = 0))
  z <- qnorm(0.5 + level / 2)
  lo <- pred_mean[valid] - z * pred_se[valid]
  hi <- pred_mean[valid] + z * pred_se[valid]
  covered <- observed[valid] >= lo & observed[valid] <= hi
  list(coverage = mean(covered), level = level, n = sum(valid))
}
