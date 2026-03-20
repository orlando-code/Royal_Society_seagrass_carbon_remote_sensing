library(gridExtra)
library(tidyverse)
library(grid)
library(ggplot2)
library(maps)
library(glmnet)

# Get world map data
world_map <- map_data("world")
min_lon <- min(dat$longitude, na.rm = TRUE)
max_lon <- max(dat$longitude, na.rm = TRUE)
min_lat <- min(dat$latitude, na.rm = TRUE)
max_lat <- max(dat$latitude, na.rm = TRUE)

plot_world_map <- function(dat) {
  ggplot() +
    geom_polygon(data = world_map, aes(x = long, y = lat, group = group),
                 fill = "lightgray", color = "white", alpha = 0.5) +
    coord_cartesian(xlim = c(min_lon, max_lon), ylim = c(min_lat, max_lat))
}

plot_n_splits_on_grid <- function(dat, split_dat_list, n = 6, world_map, min_lon, max_lon, min_lat, max_lat) {
  plots <- list()
  legends <- list()
  n <- min(n, length(split_dat_list))
  for (i in 1:n) {
    split_dat <- split_dat_list[[i]]
    dat$split_group <- ifelse(dat$id %in% split_dat$train$id, "Train", "Test")
    # Get train/test split ratio (use ratio element from split_dat if available)
    ratio <- if (!is.null(split_dat$ratio)) split_dat$ratio else {
      nrow(split_dat$train) / (nrow(split_dat$train) + nrow(split_dat$test))
    }
    # Construct plot with legend only for first plot, remove for others.
    p <- plot_world_map(dat) +
      geom_point(data = dat, aes(x = longitude, y = latitude, color = split_group), size = 1) +
      scale_color_manual(values = c("Train" = "blue", "Test" = "red")) +
      theme_minimal()
    # Extract or suppress legend
    if (i == 1) {
      p <- p + guides(color = guide_legend(override.aes = list(size = 2)))
      legend <- cowplot::get_legend(p)
      legends[[1]] <- legend
    } else {
      p <- p + theme(legend.position = "none")
    }
    plots[[i]] <- p
  }
  n_row <- ceiling(sqrt(n))
  n_col <- ceiling(n / n_row)
  # Arrange plots with shared legend at the bottom
  grid_plots <- arrangeGrob(grobs = plots, nrow = n_row, ncol = n_col)
  if (length(legends) > 0) {
    grid.arrange(grid_plots, legends[[1]], nrow = 2, heights = unit.c(unit(1, "npc") - unit(1, "lines"), unit(1, "lines")))
  } else {
    grid.arrange(grid_plots)
  }
}


### PREPROCESSING (copied from modelling/R/helpers.R)

#' Process remote sensing covariates by clipping negative values to zero.
#' Call this after loading extracted raster data and before any training or prediction.
#'
#' @param dat Data frame with remote sensing (and possibly other) covariates.
#' @param rs_covariates Character vector of column names to clip (default: common OLCI product names).
#' @return The same data frame with the specified columns clipped to \code{pmax(x, 0)}; columns not present are ignored.
process_rs_covariates <- function(dat, rs_covariates = c(
                                    "S3A_OLCI_ERRNT.20160425_20241231.L3m.CU.CHL.chlor_a.4km",
                                    "S3A_OLCI_ERRNT.20160425_20241231.L3m.CU.IOP.adg_443.4km",
                                    "S3A_OLCI_ERRNT.20160425_20241231.L3m.CU.IOP.bbp_443.4km",
                                    "S3A_OLCI_ERRNT.20160425_20241231.L3m.CU.KD.Kd_490.4km",
                                    "S3A_OLCI_ERRNT.20160425_20241231.L3m.CU.RRS.Rrs_443.4km",
                                    "S3A_OLCI_ERRNT.20160425_20241231.L3m.CU.RRS.Rrs_620.4km"
                                  )) {
  existing_vars <- intersect(rs_covariates, colnames(dat))
  if (length(existing_vars) == 0) {
    return(dat)
  }
  dat <- as.data.frame(dat)
  for (v in existing_vars) {
    if (is.numeric(dat[[v]])) {
      dat[[v]] <- pmax(dat[[v]], 0, na.rm = FALSE)
    }
  }
  dat
}

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

# Scale covariates using z-score parameters estimated on the training fold.
# Only numeric covariates are scaled.
scale_covariates_train_test <- function(train_dat, test_dat, covariate_vars) {
  sp <- compute_scale_params(train_dat, covariate_vars)
  list(
    train = apply_scaling(train_dat, sp),
    test  = apply_scaling(test_dat,  sp)
  )
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

# check if train and test data have enough rows to run a model
check_enough <- function(f, train_dat, test_dat, status = "too_few_rows") {
  if (nrow(train_dat) < 2L || nrow(test_dat) < 1L) {
    return(data.frame(
      fold = f,
      n_train = nrow(train_dat),
      n_test = nrow(test_dat),
      rmse = NA_real_,
      r2 = NA_real_,
      status = status
    ))
  }
  NULL
}

### ML

compute_r2 <- function(y, yhat) {
  ok <- is.finite(y) & is.finite(yhat)
  y <- y[ok]
  yhat <- yhat[ok]
  if (length(y) < 2L) return(NA_real_)
  ss_res <- sum((y - yhat)^2)
  ss_tot <- sum((y - mean(y))^2)
  if (ss_tot <= 0) return(NA_real_)
  1 - ss_res / ss_tot
}

compute_rmse <- function(y, yhat) {
  sqrt(mean((y - yhat)^2, na.rm = TRUE))
}

compute_rmse_r2 <- function(y, yhat) {
  list(rmse = compute_rmse(y, yhat), r2 = compute_r2(y, yhat))
}

maybe_inverse_for_metrics <- function(y, yhat, inverse_response_for_metrics) {
  if (!inverse_response_for_metrics) return(list(y = y, yhat = yhat))
  list(
    y = inverse_response_transform(y, log = TRUE),
    yhat = inverse_response_transform(yhat, log = TRUE)
  )
}

finalize_fold_metrics <- function(f, n_train, n_test, y, yhat, inverse_response_for_metrics,
                                    status = NULL) {
  inv <- maybe_inverse_for_metrics(y, yhat, inverse_response_for_metrics)
  m <- compute_rmse_r2(inv$y, inv$yhat)
  df <- data.frame(
    fold = f,
    n_train = n_train,
    n_test = n_test,
    rmse = m$rmse,
    r2 = m$r2
  )
  if (!is.null(status)) df$status <- status
  df
}

status_only_df <- function(f, status, n_train = NA_integer_, n_test = NA_integer_) {
  data.frame(
    fold = f,
    n_train = n_train,
    n_test = n_test,
    rmse = NA_real_,
    r2 = NA_real_,
    status = status
  )
}

### SPLITTING

cov_id_split <- function(dat, split_prop = 0.8, seed = 42) {
  set.seed(seed)
  unique_cov_ids <- sort(unique(na.omit(dat$id)))
  n_unique_cov <- length(unique_cov_ids) # should be ~155
  n_train_cov <- floor(0.8 * n_unique_cov)
  train_cov_ids <- sample(unique_cov_ids, size = n_train_cov, replace = FALSE)
  train <- dat[dat$id %in% train_cov_ids, , drop = FALSE]
  test  <- dat[!(dat$id %in% train_cov_ids), , drop = FALSE]
  ratio <- dim(train)[1] / (dim(train)[1] + dim(test)[1])
  return(list(train = train, test = test, ratio = ratio))
}

# Split a full dataset into train/test for a given fold.
split_fold <- function(dat, f, fold_col = "fold") {
  list(
    train = dat[dat[[fold_col]] != f, , drop = FALSE],
    test  = dat[dat[[fold_col]] == f, , drop = FALSE]
  )
}

# Keep only complete rows on (target + predictors/covariates).
prepare_fold_complete <- function(dat, f, fold_col, target_var, x_vars = NULL) {
  sp <- split_fold(dat, f, fold_col = fold_col)
  needed <- unique(c(target_var, x_vars))
  sp$train <- sp$train[complete.cases(sp$train[, needed, drop = FALSE]), , drop = FALSE]
  sp$test  <- sp$test[complete.cases(sp$test[, needed, drop = FALSE]), , drop = FALSE]
  sp
}

align_factor_levels <- function(train_dat, test_dat, var) {
  train_dat[[var]] <- factor(train_dat[[var]])
  test_dat[[var]] <- factor(test_dat[[var]], levels = levels(train_dat[[var]]))
  list(train = train_dat, test = test_dat)
}


### MODELS

evaluate_fold_lm <- function(f, dat, target_var, predictor_vars, lm_formula) {
  train_dat <- dat[dat$fold != f, , drop = FALSE]
  test_dat  <- dat[dat$fold == f, , drop = FALSE]

  needed_cols <- unique(c(target_var, predictor_vars))
  train_dat <- train_dat[complete.cases(train_dat[, needed_cols, drop = FALSE]), , drop = FALSE]
  test_dat  <- test_dat[complete.cases(test_dat[, needed_cols, drop = FALSE]), , drop = FALSE]

  enough <- check_enough(f, train_dat, test_dat)
  if (!is.null(enough)) return(enough)

  # Ensure species factor is consistent across train/test.
  if ("seagrass_species" %in% predictor_vars) {
    test_dat$seagrass_species <- factor(test_dat$seagrass_species,
                                         levels = unique(train_dat$seagrass_species))
  }

  model <- lm(lm_formula, data = train_dat)
  y <- test_dat[[target_var]]
  yhat <- as.numeric(predict(model, newdata = test_dat))

  m <- compute_rmse_r2(y, yhat)

  data.frame(fold = f, n_train = nrow(train_dat), n_test = nrow(test_dat),
             rmse = m$rmse, r2 = m$r2)
}

evaluate_fold_null_mean <- function(f, dat, target_var, fold_col = "fold",
                                     inverse_response_for_metrics = FALSE) {
  prep <- prepare_fold_complete(dat, f, fold_col = fold_col, target_var = target_var, x_vars = NULL)
  train_dat <- prep$train
  test_dat  <- prep$test

  enough <- check_enough(f, train_dat, test_dat)
  if (!is.null(enough)) return(enough)

  y <- test_dat[[target_var]]
  yhat <- rep(mean(train_dat[[target_var]], na.rm = TRUE), nrow(test_dat))
  finalize_fold_metrics(
    f = f,
    n_train = nrow(train_dat),
    n_test = nrow(test_dat),
    y = y,
    yhat = yhat,
    inverse_response_for_metrics = inverse_response_for_metrics
  )
}

evaluate_fold_lm_covariates <- function(f, dat, target_var, covariate_vars, fold_col = "fold",
                                         inverse_response_for_metrics = FALSE) {
  if (length(covariate_vars) < 1L) stop("evaluate_fold_lm_covariates: covariate_vars is empty")

  prep <- prepare_fold_complete(dat, f, fold_col = fold_col, target_var = target_var, x_vars = covariate_vars)
  train_dat <- prep$train
  test_dat  <- prep$test

  enough <- check_enough(f, train_dat, test_dat)
  if (!is.null(enough)) return(enough)

  st <- scale_covariates_train_test(train_dat, test_dat, covariate_vars)
  train_dat <- st$train
  test_dat  <- st$test

  fml <- as.formula(paste(target_var, "~", paste(covariate_vars, collapse = " + ")))
  model <- lm(fml, data = train_dat)

  y <- test_dat[[target_var]]
  yhat <- as.numeric(predict(model, newdata = test_dat))
  finalize_fold_metrics(
    f = f,
    n_train = nrow(train_dat),
    n_test = nrow(test_dat),
    y = y,
    yhat = yhat,
    inverse_response_for_metrics = inverse_response_for_metrics
  )
}

evaluate_fold_gam <- function(f, dat, target_var, covariate_vars, species_var = NULL, fold_col = "fold",
                               gam_k = 6L, gam_bs = "cr", gam_select = TRUE, gam_method = "REML",
                               inverse_response_for_metrics = FALSE) {
  if (!requireNamespace("mgcv", quietly = TRUE)) {
    return(status_only_df(f, status = "gam_missing"))
  }
  if (length(covariate_vars) < 1L) stop("evaluate_fold_gam: covariate_vars is empty")

  # Only numeric covariates can appear in s(); keep formula well-defined.
  numeric_covs <- covariate_vars[vapply(dat[covariate_vars], is.numeric, logical(1))]
  if (length(numeric_covs) < 1L) {
    return(status_only_df(f, status = "gam_no_numeric_covariates"))
  }

  x_vars <- numeric_covs
  if (!is.null(species_var) && species_var %in% names(dat)) x_vars <- unique(c(x_vars, species_var))

  prep <- prepare_fold_complete(dat, f, fold_col = fold_col, target_var = target_var, x_vars = x_vars)
  train_dat <- prep$train
  test_dat  <- prep$test

  enough <- check_enough(f, train_dat, test_dat, status = "too_few_rows")
  if (!is.null(enough)) return(enough)

  # Scale numeric covariates (helps smoothers be more stable across folds).
  st <- scale_covariates_train_test(train_dat, test_dat, numeric_covs)
  train_dat <- st$train
  test_dat  <- st$test

  if (!is.null(species_var) && species_var %in% names(dat)) {
    train_dat[[species_var]] <- factor(train_dat[[species_var]])
    test_dat[[species_var]]  <- factor(test_dat[[species_var]], levels = levels(train_dat[[species_var]]))
  }

  smooth_terms <- paste0("s(", numeric_covs, ", k = ", gam_k, ", bs = '", gam_bs, "')")
  fixed_terms <- paste(smooth_terms, collapse = " + ")
  if (!is.null(species_var) && species_var %in% names(dat)) {
    fml <- as.formula(paste0(target_var, " ~ ", fixed_terms, " + ", species_var))
  } else {
    fml <- as.formula(paste0(target_var, " ~ ", fixed_terms))
  }

  fit <- tryCatch(
    {
      mgcv::gam(
        formula = fml,
        data = train_dat,
        family = gaussian(),
        method = gam_method,
        select = gam_select
      )
    },
    error = function(e) NULL
  )

  if (is.null(fit)) {
    return(status_only_df(f, status = "gam_fit_failed",
                           n_train = nrow(train_dat), n_test = nrow(test_dat)))
  }

  preds <- tryCatch(
    as.numeric(predict(fit, newdata = test_dat)),
    error = function(e) rep(NA_real_, nrow(test_dat))
  )

  y_test <- test_dat[[target_var]]
  finalize_fold_metrics(
    f = f,
    n_train = nrow(train_dat),
    n_test = nrow(test_dat),
    y = y_test,
    yhat = preds,
    inverse_response_for_metrics = inverse_response_for_metrics,
    status = "ok"
  )
}

evaluate_fold_species_only_lm <- function(f, dat, target_var, species_var, fold_col = "fold",
                                           inverse_response_for_metrics = FALSE) {
  prep <- prepare_fold_complete(dat, f, fold_col = fold_col, target_var = target_var, x_vars = species_var)
  train_dat <- prep$train
  test_dat  <- prep$test

  enough <- check_enough(f, train_dat, test_dat)
  if (!is.null(enough)) return(enough)

  train_dat[[species_var]] <- factor(train_dat[[species_var]])
  test_dat[[species_var]] <- factor(test_dat[[species_var]], levels = levels(train_dat[[species_var]]))

  fml <- as.formula(paste(target_var, "~", species_var))
  model <- lm(fml, data = train_dat)

  y <- test_dat[[target_var]]
  yhat <- as.numeric(predict(model, newdata = test_dat))
  finalize_fold_metrics(
    f = f,
    n_train = nrow(train_dat),
    n_test = nrow(test_dat),
    y = y,
    yhat = yhat,
    inverse_response_for_metrics = inverse_response_for_metrics
  )
}

evaluate_fold_elastic_net <- function(f, dat, target_var, covariate_vars, fold_col = "fold",
                                       alpha = 0.5, nfolds_inner = 5L, glmnet_seed = 42,
                                       inverse_response_for_metrics = FALSE) {
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    return(status_only_df(f, status = "glmnet_missing"))
  }

  prep <- prepare_fold_complete(dat, f, fold_col = fold_col, target_var = target_var, x_vars = covariate_vars)
  train_dat <- prep$train
  test_dat  <- prep$test

  enough <- check_enough(f, train_dat, test_dat, status = "too_few_rows")
  if (!is.null(enough)) return(enough)
  if (length(covariate_vars) < 1L) stop("evaluate_fold_elastic_net: covariate_vars is empty")

  X_train <- as.matrix(train_dat[, covariate_vars, drop = FALSE])
  X_test  <- as.matrix(test_dat[, covariate_vars, drop = FALSE])
  storage.mode(X_train) <- "double"
  storage.mode(X_test)  <- "double"
  y_train <- as.numeric(train_dat[[target_var]])

  set.seed(glmnet_seed + f)
  fit <- tryCatch(
    glmnet::cv.glmnet(
      x = X_train, y = y_train,
      alpha = alpha, nfolds = as.integer(nfolds_inner),
      family = "gaussian", standardize = TRUE
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(status_only_df(
      f,
      status = "glmnet_fit_failed",
      n_train = nrow(train_dat),
      n_test = nrow(test_dat)
    ))
  }

  yhat <- as.numeric(predict(fit, newx = X_test, s = "lambda.min"))
  y <- test_dat[[target_var]]
  finalize_fold_metrics(
    f = f,
    n_train = nrow(train_dat),
    n_test = nrow(test_dat),
    y = y,
    yhat = yhat,
    inverse_response_for_metrics = inverse_response_for_metrics,
    status = "ok"
  )
}

evaluate_fold_random_forest <- function(f, dat, target_var, covariate_vars, fold_col = "fold",
                                         rf_ntree = 300L, rf_seed = 42,
                                         inverse_response_for_metrics = FALSE) {
  if (!requireNamespace("randomForest", quietly = TRUE)) {
    return(status_only_df(f, status = "rf_missing"))
  }
  prep <- prepare_fold_complete(dat, f, fold_col = fold_col, target_var = target_var, x_vars = covariate_vars)
  train_dat <- prep$train
  test_dat  <- prep$test
  enough <- check_enough(f, train_dat, test_dat, status = "too_few_rows")
  if (!is.null(enough)) return(enough)
  if (length(covariate_vars) < 1L) stop("evaluate_fold_random_forest: covariate_vars is empty")

  st <- scale_covariates_train_test(train_dat, test_dat, covariate_vars)
  train_dat <- st$train
  test_dat  <- st$test

  X_train <- train_dat[, covariate_vars, drop = FALSE]
  X_test  <- test_dat[, covariate_vars, drop = FALSE]
  y_train <- train_dat[[target_var]]
  y_test  <- test_dat[[target_var]]

  set.seed(rf_seed + f)
  p <- max(2L, floor(sqrt(length(covariate_vars))))
  fit <- tryCatch(
    randomForest::randomForest(
      x = X_train, y = y_train,
      ntree = as.integer(rf_ntree),
      mtry = p,
      importance = FALSE
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(status_only_df(f, status = "rf_fit_failed",
                           n_train = nrow(train_dat), n_test = nrow(test_dat)))
  }

  yhat <- as.numeric(predict(fit, newdata = X_test))
  finalize_fold_metrics(
    f = f,
    n_train = nrow(train_dat),
    n_test = nrow(test_dat),
    y = y_test,
    yhat = yhat,
    inverse_response_for_metrics = inverse_response_for_metrics,
    status = "ok"
  )
}

evaluate_fold_xgb <- function(f, dat, target_var, covariate_vars, fold_col = "fold",
                               nrounds = 200L, xgb_seed = 42,
                               inverse_response_for_metrics = FALSE) {
  if (!requireNamespace("xgboost", quietly = TRUE)) {
    return(status_only_df(f, status = "xgb_missing"))
  }
  prep <- prepare_fold_complete(dat, f, fold_col = fold_col, target_var = target_var, x_vars = covariate_vars)
  train_dat <- prep$train
  test_dat  <- prep$test

  enough <- check_enough(f, train_dat, test_dat, status = "too_few_rows")
  if (!is.null(enough)) return(enough)
  if (length(covariate_vars) < 1L) stop("evaluate_fold_xgb: covariate_vars is empty")

  st <- scale_covariates_train_test(train_dat, test_dat, covariate_vars)
  train_dat <- st$train
  test_dat  <- st$test

  X_train <- as.matrix(train_dat[, covariate_vars, drop = FALSE])
  X_test  <- as.matrix(test_dat[, covariate_vars, drop = FALSE])
  storage.mode(X_train) <- "double"
  storage.mode(X_test)  <- "double"
  y_train <- as.numeric(train_dat[[target_var]])
  y_test  <- test_dat[[target_var]]

  # drop any rows with non-finite predictors/targets
  ok_train <- is.finite(y_train) & apply(X_train, 1, function(r) all(is.finite(r)))
  ok_test  <- is.finite(y_test)  & apply(X_test,  1, function(r) all(is.finite(r)))
  X_train <- X_train[ok_train, , drop = FALSE]
  y_train <- y_train[ok_train]
  X_test  <- X_test[ok_test, , drop = FALSE]
  y_test  <- as.numeric(y_test)[ok_test]

  if (nrow(X_train) < 2L || nrow(X_test) < 1L) {
    return(status_only_df(
      f,
      status = "nonfinite_removed_too_few_rows",
      n_train = nrow(X_train),
      n_test = nrow(X_test)
    ))
  }

  objective <- "reg:squarederror"
  max_depth <- 3L
  learning_rate <- 0.05
  subsample <- 0.8
  colsample_bytree <- 0.8
  min_child_weight <- 1
  reg_lambda <- 1

  set.seed(xgb_seed + f)
  last_err <- NULL
  fit <- tryCatch(
    xgboost::xgboost(
      x = X_train,
      y = y_train,
      nrounds = as.integer(nrounds),
      nthread = 1L,
      objective = objective,
      max_depth = max_depth,
      learning_rate = learning_rate,
      subsample = subsample,
      colsample_bytree = colsample_bytree,
      min_child_weight = min_child_weight,
      reg_lambda = reg_lambda
    ),
    error = function(e) {
      last_err <<- conditionMessage(e)
      NULL
    }
  )
  if (is.null(fit)) {
    err_short <- if (!is.null(last_err)) substr(last_err, 1, 80) else "unknown"
    return(status_only_df(
      f,
      status = paste0("xgb_fit_failed:", err_short),
      n_train = nrow(train_dat),
      n_test = nrow(test_dat)
    ))
  }

  yhat <- as.numeric(predict(fit, newdata = X_test))
  finalize_fold_metrics(
    f = f,
    n_train = nrow(train_dat),
    n_test = nrow(test_dat),
    y = y_test,
    yhat = yhat,
    inverse_response_for_metrics = inverse_response_for_metrics,
    status = "ok"
  )
}

evaluate_fold_gpr_numeric <- function(f, dat, target_var, covariate_vars, fold_col = "fold",
                                       kernel = "matern52", nug_min = 1e-6, nug_max = 50, nug_est = TRUE,
                                       gpr_seed = 42,
                                       inverse_response_for_metrics = FALSE) {
  if (!requireNamespace("GauPro", quietly = TRUE)) {
    return(status_only_df(f, status = "gpr_missing"))
  }
  prep <- prepare_fold_complete(dat, f, fold_col = fold_col, target_var = target_var, x_vars = covariate_vars)
  train_dat <- prep$train
  test_dat  <- prep$test

  enough <- check_enough(f, train_dat, test_dat, status = "too_few_rows")
  if (!is.null(enough)) return(enough)
  if (length(covariate_vars) < 1L) stop("evaluate_fold_gpr_numeric: covariate_vars is empty")

  # Standardize numeric predictors (GPR is sensitive to scale).
  X_train_raw <- train_dat[, covariate_vars, drop = FALSE]
  X_test_raw  <- test_dat[, covariate_vars, drop = FALSE]
  means <- colMeans(X_train_raw, na.rm = TRUE)
  sds <- apply(X_train_raw, 2, sd, na.rm = TRUE)
  sds[!is.finite(sds) | sds <= 0] <- 1

  X_train <- scale(X_train_raw, center = means, scale = sds)
  X_test  <- scale(X_test_raw, center = means, scale = sds)

  train_df <- data.frame(train_dat[[target_var]])
  names(train_df)[1] <- target_var
  train_df <- cbind(train_df, as.data.frame(X_train))
  test_df <- data.frame(test_dat[[target_var]])
  names(test_df)[1] <- target_var
  test_df <- cbind(test_df, as.data.frame(X_test))

  form <- as.formula(paste(target_var, "~", paste(covariate_vars, collapse = " + ")))
  set.seed(gpr_seed + f)
  fit <- tryCatch(
    GauPro::gpkm(form, data = train_df,
                 kernel = kernel, nug.min = nug_min, nug.max = nug_max, nug.est = nug_est),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(status_only_df(f, status = "gpr_fit_failed",
                           n_train = nrow(train_dat), n_test = nrow(test_dat)))
  }

  X_test_mat <- as.matrix(test_df[, covariate_vars, drop = FALSE])
  preds <- as.numeric(fit$pred(X_test_mat, se.fit = FALSE))
  y_test <- test_dat[[target_var]]
  finalize_fold_metrics(
    f = f,
    n_train = nrow(train_dat),
    n_test = nrow(test_dat),
    y = y_test,
    yhat = preds,
    inverse_response_for_metrics = inverse_response_for_metrics,
    status = "ok"
  )
}

evaluate_fold_mixed_species_random <- function(f, dat, target_var, covariate_vars, species_var,
                                                fold_col = "fold",
                                                inverse_response_for_metrics = FALSE) {
  if (!requireNamespace("lme4", quietly = TRUE)) {
    return(status_only_df(f, status = "lme4_missing"))
  }
  prep <- prepare_fold_complete(dat, f, fold_col = fold_col, target_var = target_var,
                                  x_vars = c(covariate_vars, species_var))
  train_dat <- prep$train
  test_dat  <- prep$test

  enough <- check_enough(f, train_dat, test_dat, status = "too_few_rows")
  if (!is.null(enough)) return(enough)

  train_dat[[species_var]] <- factor(train_dat[[species_var]])
  test_dat[[species_var]] <- factor(test_dat[[species_var]], levels = levels(train_dat[[species_var]]))

  st <- scale_covariates_train_test(train_dat, test_dat, covariate_vars)
  train_dat <- st$train
  test_dat  <- st$test

  fixed_terms <- if (length(covariate_vars) > 0L) paste(covariate_vars, collapse = " + ") else "1"
  form <- as.formula(paste0(target_var, " ~ ", fixed_terms, " + (1|", species_var, ")"))

  fit <- tryCatch(
    lme4::lmer(form, data = train_dat, REML = FALSE),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(status_only_df(f, status = "lmer_fit_failed",
                           n_train = nrow(train_dat), n_test = nrow(test_dat)))
  }

  preds <- tryCatch(
    as.numeric(predict(fit, newdata = test_dat, allow.new.levels = TRUE)),
    error = function(e) rep(NA_real_, nrow(test_dat))
  )

  y_test <- test_dat[[target_var]]
  finalize_fold_metrics(
    f = f,
    n_train = nrow(train_dat),
    n_test = nrow(test_dat),
    y = y_test,
    yhat = preds,
    inverse_response_for_metrics = inverse_response_for_metrics,
    status = "ok"
  )
}


### WRAPPERS

run_cv_simple_model_suite <- function(dat, target_var, covariate_vars, species_var,
                                       fold_col = "fold",
                                       models = c("null_mean", "linear", "elastic_net", "random_forest", "xgb", "gpr", "gam",
                                                  "mixed_species_random", "species_only_lm"),
                                       glmnet_alpha = 0.5, glmnet_inner_folds = 5L,
                                       rf_ntree = 300L,
                                       xgb_nrounds = 200L,
                                       gpr_kernel = "matern52",
                                       gpr_nug_min = 1e-6, gpr_nug_max = 50,
                                       gam_k = 6L, gam_bs = "cr", gam_select = TRUE, gam_method = "REML",
                                       inverse_response_for_metrics = FALSE) {

  if (!fold_col %in% names(dat)) stop("run_cv_simple_model_suite: fold column not found in dat")
  fold_ids <- sort(unique(na.omit(dat[[fold_col]])))
  if (length(fold_ids) < 2L) stop("run_cv_simple_model_suite: need >= 2 folds")

  if (is.null(species_var) || !species_var %in% names(dat)) species_var <- NULL

  eval_fns <- list(
    null_mean = function(f) {
      evaluate_fold_null_mean(f, dat, target_var, fold_col = fold_col,
                               inverse_response_for_metrics = inverse_response_for_metrics)
    },
    linear = function(f) {
      evaluate_fold_lm_covariates(f, dat, target_var, covariate_vars, fold_col = fold_col,
                                   inverse_response_for_metrics = inverse_response_for_metrics)
    },
    elastic_net = function(f) {
      evaluate_fold_elastic_net(
        f, dat, target_var, covariate_vars, fold_col = fold_col,
        alpha = glmnet_alpha, nfolds_inner = glmnet_inner_folds, glmnet_seed = 42,
        inverse_response_for_metrics = inverse_response_for_metrics
      )
    },
    random_forest = function(f) {
      evaluate_fold_random_forest(
        f, dat, target_var, covariate_vars, fold_col = fold_col,
        rf_ntree = rf_ntree, rf_seed = 42,
        inverse_response_for_metrics = inverse_response_for_metrics
      )
    },
    xgb = function(f) {
      evaluate_fold_xgb(
        f, dat, target_var, covariate_vars, fold_col = fold_col,
        nrounds = xgb_nrounds, xgb_seed = 42,
        inverse_response_for_metrics = inverse_response_for_metrics
      )
    },
    gpr = function(f) {
      evaluate_fold_gpr_numeric(
        f, dat, target_var, covariate_vars, fold_col = fold_col,
        kernel = gpr_kernel, nug_min = gpr_nug_min, nug_max = gpr_nug_max, nug_est = TRUE,
        gpr_seed = 42,
        inverse_response_for_metrics = inverse_response_for_metrics
      )
    },
    gam = function(f) {
      evaluate_fold_gam(
        f, dat, target_var, covariate_vars, species_var = species_var, fold_col = fold_col,
        gam_k = gam_k, gam_bs = gam_bs, gam_select = gam_select, gam_method = gam_method,
        inverse_response_for_metrics = inverse_response_for_metrics
      )
    },
    mixed_species_random = function(f) {
      if (is.null(species_var)) {
        data.frame(fold = f, n_train = NA_integer_, n_test = NA_integer_,
                   rmse = NA_real_, r2 = NA_real_, status = "species_missing")
      } else {
        evaluate_fold_mixed_species_random(
          f, dat, target_var, covariate_vars, species_var, fold_col = fold_col,
          inverse_response_for_metrics = inverse_response_for_metrics
        )
      }
    },
    species_only_lm = function(f) {
      if (is.null(species_var)) {
        data.frame(fold = f, n_train = NA_integer_, n_test = NA_integer_,
                   rmse = NA_real_, r2 = NA_real_, status = "species_missing")
      } else {
        evaluate_fold_species_only_lm(
          f, dat, target_var, species_var, fold_col = fold_col,
          inverse_response_for_metrics = inverse_response_for_metrics
        )
      }
    }
  )

  eval_one <- function(model_name) {
    fn <- eval_fns[[model_name]]
    if (is.null(fn)) stop("Unknown model: ", model_name)

    out_by_fold <- lapply(fold_ids, function(f) {
      tmp <- fn(f)

      # Normalize schema so `rbind()` never fails.
      if (!"status" %in% names(tmp)) tmp$status <- "ok"
      tmp$model <- model_name
      required_cols <- c("model", "fold", "n_train", "n_test", "rmse", "r2", "status")
      for (col in required_cols) if (!col %in% names(tmp)) tmp[[col]] <- NA
      tmp <- tmp[, required_cols, drop = FALSE]
      tmp
    })
    do.call(rbind, out_by_fold)
  }

  per_fold <- do.call(rbind, lapply(models, eval_one))

  summary <- per_fold %>%
    dplyr::group_by(model) %>%
    dplyr::summarise(
      mean_rmse = mean(rmse, na.rm = TRUE),
      sd_rmse = stats::sd(rmse, na.rm = TRUE),
      mean_r2 = mean(r2, na.rm = TRUE),
      sd_r2 = stats::sd(r2, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(mean_rmse)

  list(per_fold = per_fold, summary = summary)
}

# Run the same CV suite multiple times with different fold-assignment seeds.
# Fold assignment is over unique values of `id_col` (unique covariate-vector clusters),
# ensuring that changing the seed actually changes train/test partitions.
run_cv_model_suite_over_seed_list <- function(dat, target_var, covariate_vars, species_var,
                                                id_col = "id",
                                                k_folds = 5L,
                                                seed_list = 1:20,
                                                fold_col = "fold",
                                                models = c("null_mean", "linear", "elastic_net", "random_forest", "xgb", "gpr", "gam",
                                                           "mixed_species_random", "species_only_lm"),
                                                inverse_response_for_metrics = FALSE,
                                                progress_every = 1L,
                                                ...) {
  if (!id_col %in% names(dat)) stop("run_cv_model_suite_over_seed_list: id_col not found in dat")
  fold_ids <- sort(unique(na.omit(dat[[id_col]])))
  if (length(fold_ids) < 2L) stop("run_cv_model_suite_over_seed_list: need >= 2 unique id clusters")

  seed_list <- as.integer(seed_list)
  n_seeds <- length(seed_list)
  if (n_seeds < 1L) stop("run_cv_model_suite_over_seed_list: seed_list is empty")

  cat(sprintf("\nRunning CV suite over %d seed(s) (k_folds=%d, id clusters=%d)\n",
              n_seeds, as.integer(k_folds), length(fold_ids)))

  per_seed_summary <- vector("list", n_seeds)
  for (i in seq_len(n_seeds)) {
    seed <- seed_list[i]
    if (progress_every > 0L && (i %% progress_every == 0L || i == 1L || i == n_seeds)) {
      cat(sprintf("  Seed %d/%d: %d\n", i, n_seeds, seed))
    }

    set.seed(seed)
    # ensure at least 2 distinct folds are present, else the downstream CV
    # would fail when it expects multiple fold IDs.
    fold_assignment <- NULL
    for (attempt in 1:10) {
      fold_assignment <- as.integer(sample(rep(seq_len(k_folds), length.out = length(fold_ids))))
      if (length(unique(fold_assignment)) >= 2L) break
    }
    names(fold_assignment) <- as.character(fold_ids)

    dat_seed <- dat
    dat_seed[[fold_col]] <- as.integer(fold_assignment[as.character(dat_seed[[id_col]])])

    res <- run_cv_simple_model_suite(
      dat = dat_seed,
      target_var = target_var,
      covariate_vars = covariate_vars,
      species_var = species_var,
      fold_col = fold_col,
      models = models,
      inverse_response_for_metrics = inverse_response_for_metrics,
      ...
    )

    s <- res$summary
    s$seed <- seed
    per_seed_summary[[i]] <- s
  }

  per_seed_summary <- dplyr::bind_rows(per_seed_summary)

  summary_over_seeds <- per_seed_summary %>%
    dplyr::group_by(model) %>%
    dplyr::summarise(
      mean_mean_rmse = mean(mean_rmse, na.rm = TRUE),
      sd_mean_rmse = stats::sd(mean_rmse, na.rm = TRUE),
      mean_mean_r2 = mean(mean_r2, na.rm = TRUE),
      sd_mean_r2 = stats::sd(mean_r2, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(desc(mean_mean_r2))

  list(
    per_seed_summary = per_seed_summary,
    summary_over_seeds = summary_over_seeds
  )
}
