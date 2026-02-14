#' Prune predictors by pairwise correlation threshold
#'
#' Removes one variable from each pair with |correlation| > threshold, keeping the one with higher absolute correlation to the target.
#' @param data data.frame with predictors and target
#' @param predictor_vars character vector of predictor names
#' @param target_var name of target variable
#' @param cor_threshold absolute correlation threshold (default 0.8)
#' @return character vector of pruned predictor names
prune_by_correlation <- function(data, predictor_vars, target_var, cor_threshold = 0.8) {
  # Defensive: ensure target_var is a single character string
  if (is.list(target_var) || length(target_var) != 1) {
    target_var <- as.character(target_var)[1]
  }
  dat <- data[, c(target_var, predictor_vars), drop = FALSE]
  dat <- dat[complete.cases(dat), ]
  if (nrow(dat) < 3 || length(predictor_vars) < 2) {
    return(predictor_vars)
  }
  cor_with_target <- sapply(predictor_vars, function(v) {
    if (v %in% names(dat)) {
      stats::cor(dat[[v]], dat[[target_var]], use = "complete.obs")
    } else {
      NA_real_
    }
  })
  cor_df <- data.frame(variable = predictor_vars, correlation = cor_with_target, abs_correlation = abs(cor_with_target), stringsAsFactors = FALSE)
  cor_df <- cor_df[!is.na(cor_df$correlation), ]
  vars_for_cor <- cor_df$variable
  if (length(vars_for_cor) < 2) {
    return(vars_for_cor)
  }
  covar_data <- dat[, vars_for_cor, drop = FALSE]
  cor_matrix <- stats::cor(covar_data, use = "complete.obs")
  high_cor_pairs <- which(abs(cor_matrix) > cor_threshold & abs(cor_matrix) < 1, arr.ind = TRUE)
  vars_removed <- character()
  if (nrow(high_cor_pairs) > 0) {
    pair_df <- data.frame(
      var1 = rownames(cor_matrix)[high_cor_pairs[, 1]],
      var2 = colnames(cor_matrix)[high_cor_pairs[, 2]],
      correlation = cor_matrix[high_cor_pairs],
      stringsAsFactors = FALSE
    )
    pair_df <- pair_df[pair_df$var1 < pair_df$var2, ]
    pair_df <- pair_df[order(abs(pair_df$correlation), decreasing = TRUE), ]
    for (i in seq_len(nrow(pair_df))) {
      v1 <- pair_df$var1[i]
      v2 <- pair_df$var2[i]
      if (v1 %in% vars_removed || v2 %in% vars_removed) next
      c1 <- cor_df$abs_correlation[cor_df$variable == v1]
      c2 <- cor_df$abs_correlation[cor_df$variable == v2]
      if (length(c1) == 0 || length(c2) == 0) next
      if (round(c1, 3) == round(c2, 3)) next
      if (c1 < c2) {
        vars_removed <- c(vars_removed, v1)
      } else {
        vars_removed <- c(vars_removed, v2)
      }
    }
  }
  setdiff(vars_for_cor, unique(vars_removed))
}
# ================================ UTILITY FUNCTIONS ================================

#' Install packages if not already installed
#'
#' @param packages character vector of package names to install
#' @return NULL but prints message to console when packages are installed
install_packages <- function(packages) {
  # if string, convert to list
  if (is.character(packages)) {
    packages <- list(packages)
  }
  for (package in packages) {
    if (!requireNamespace(package, quietly = TRUE)) {
      install.packages(package)
    }
  }
  print(paste("Installed packages:", paste(packages, collapse = ", ")))
}

#' Load packages, install if necessary
#'
#' @param packages character vector of package names to load
#' @param optional Optional character vector of package names that are not required (warn if missing but continue)
#' @return NULL but prints message to console when packages are loaded
load_packages <- function(packages, optional = NULL) {
  loaded <- character()
  failed <- character()

  for (package in packages) {
    if (package %in% optional) {
      # Optional package: try to load but don't fail if unavailable
      if (suppressWarnings(require(package, character.only = TRUE, quietly = TRUE))) {
        loaded <- c(loaded, package)
      } else {
        # Try to install
        if (!requireNamespace(package, quietly = TRUE)) {
          tryCatch(
            {
              install.packages(package, quiet = TRUE)
            },
            error = function(e) NULL
          )
        }
        # Try to load again
        if (suppressWarnings(require(package, character.only = TRUE, quietly = TRUE))) {
          loaded <- c(loaded, package)
        } else {
          warning("Optional package '", package, "' could not be loaded. Continuing without it.")
        }
      }
    } else {
      # Required package: must load successfully
      if (!suppressWarnings(require(package, character.only = TRUE, quietly = TRUE))) {
        # Try to install
        if (!requireNamespace(package, quietly = TRUE)) {
          tryCatch(
            {
              install.packages(package, quiet = TRUE)
            },
            error = function(e) {
              stop("Failed to install required package '", package, "': ", conditionMessage(e))
            }
          )
        }
        # Try to load
        if (!suppressWarnings(require(package, character.only = TRUE, quietly = TRUE))) {
          stop("Required package '", package, "' could not be loaded after installation attempt.")
        }
      }
      loaded <- c(loaded, package)
    }
  }

  if (length(loaded) > 0) {
    cat("Loaded packages:", paste(loaded, collapse = ", "), "\n")
  }
  if (length(failed) > 0) {
    warning("Failed to load packages:", paste(failed, collapse = ", "))
  }
}

# ================================ PREDICTOR SETS FROM MODEL PERMUTATION PRUNING ================================
# Generic model-wise permutation pruning (model_permutation_pruning.R) writes a CSV with
# model, variable, rmse_increase, cum_rmse, total_rmse, cum_frac, keep, permutation_keep_frac.
# Use this for "best covariate set" per model instead of model-specific pruning CSVs.

#' Get pruned predictor names for a model from the generic permutation-pruning CSV
#'
#' Looks for \code{cv_model_permutation_pruning_*.csv} in \code{figures/cv_pipeline_output}
#' and \code{modelling/cv_pipeline_output}, filters to the given model and \code{keep == TRUE},
#' and returns the variable names (optionally intersected with \code{data_cols}).
#'
#' @param model_name Character, e.g. \code{"GPR"}.
#' @param data_cols Optional character vector of column names; if provided, returned
#'   variables are intersected with this (so only names present in your data are returned).
#' @param search_dirs Directories to search for the CSV (default: figures then modelling output).
#' @return Character vector of variable names, or \code{NULL} if no CSV found or no kept variables.
get_pruned_predictors_for_model <- function(model_name,
                                            data_cols = NULL,
                                            search_dirs = c("figures/cv_pipeline_output", "modelling/cv_pipeline_output")) {
  pattern <- "cv_model_permutation_pruning_.*\\.csv"
  csv_path <- NULL
  for (d in search_dirs) {
    if (!dir.exists(d)) next
    f <- list.files(d, pattern = pattern, full.names = TRUE)
    if (length(f) > 0L) {
      csv_path <- f[1L]
      break
    }
  }
  if (is.null(csv_path)) {
    return(NULL)
  }
  df <- tryCatch(
    read.csv(csv_path, stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  if (is.null(df) || !all(c("model", "variable", "keep") %in% names(df))) {
    return(NULL)
  }
  keep_ok <- if (is.logical(df$keep)) df$keep else (toupper(as.character(df$keep)) %in% c("TRUE", "T", "1"))
  vars <- unique(df$variable[df$model == model_name & keep_ok])
  if (length(vars) == 0L) {
    return(NULL)
  }
  if (!is.null(data_cols)) vars <- intersect(vars, data_cols)
  vars
}

#' Combine per-model pruned variable CSVs into a single file.
#'
#' Reads all \code{pruned_variables_to_include_<model>.csv} files in
#' \code{covariate_dir} and writes \code{pruned_model_variables.csv} with
#' columns \code{model} and \code{variable}. Call after writing any
#' pruned_variables_to_include_*.csv so the combined file stays up to date.
#'
#' @param covariate_dir Directory containing the per-model CSVs (default:
#'   \code{figures/covariate_selection}).
#' @return Invisibly the path to the combined CSV, or \code{NULL} if no
#'   per-model files were found.
combine_pruned_model_variables <- function(covariate_dir = "figures/covariate_selection") {
  pruned_files <- list.files(covariate_dir, pattern = "^pruned_variables_to_include_.+\\.csv$", full.names = TRUE)
  pruned_files <- setdiff(pruned_files, file.path(covariate_dir, "pruned_model_variables.csv"))
  if (length(pruned_files) == 0L) return(invisible(NULL))
  combined_list <- lapply(pruned_files, function(f) {
    model_name <- sub("^pruned_variables_to_include_(.+)\\.csv$", "\\1", basename(f))
    d <- read.csv(f, stringsAsFactors = FALSE)
    if (!"variable" %in% names(d)) return(NULL)
    data.frame(model = model_name, variable = d$variable, stringsAsFactors = FALSE)
  })
  combined_list <- combined_list[!vapply(combined_list, is.null, logical(1L))]
  if (length(combined_list) == 0L) return(invisible(NULL))
  combined <- do.call(rbind, combined_list)
  out_combined <- file.path(covariate_dir, "pruned_model_variables.csv")
  write.csv(combined, out_combined, row.names = FALSE)
  out_combined
}

# ================================ DATA PREPARATION ================================

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

# ================================ INITIALISATION ================================

# load packages if not already loaded
load_packages(c("here", "mgcv", "tidyverse", "ggplot2", "visreg", "randomForest", "corrplot", "blockCV", "GauPro", "xgboost", "e1071", "neuralnet", "purrr", "sp", "viridisLite", "patchwork", "maps"))
world <- map_data("world")


#' Create a progress bar for hyperparameter tuning
#'
#' @param total total number of iterations for the progress bar
#' @param format format string for the progress bar display
#' @param clear whether to clear the progress bar when done
#' @param width width of the progress bar in characters
#' @return a progress bar object from the progress package
progress_bar <- function(total, format = "[:bar] :current/:total (:percent)", clear = FALSE, width = 60) {
  if (!requireNamespace("progress", quietly = TRUE)) {
    install.packages("progress")
  }
  progress::progress_bar$new(total = total, format = format, clear = clear, width = width)
}

# ================================ MODEL FITTING ================================

#' Calculate the aggregate carbon density for each core at a given depth
#'
#' @param core_data data frame containing the core data
#' @param aggregation_level the depth level to which to aggregate
#' @param aggregation_function the function to use to aggregate the carbon density (e.g. median, mean, sum)
#' @return data frame containing the aggregate carbon density for each core
calculate_core_aggregate <- function(core_data, aggregation_level = 5, aggregation_function = median) {
  core_data %>%
    filter(sediment_mean_depth_cm <= aggregation_level) %>%
    group_by(random_core_variable) %>%
    summarise(
      aggregate_carbon_density = aggregation_function(carbon_density_g_c_cm3, na.rm = TRUE),
      .groups = "drop"
    )
}


#' Calculate performance metrics for model predictions
#'
#' @param observed vector of observed (true) values
#' @param predicted vector of predicted values
#' @return data frame with r2, rmse, mae, and bias metrics
calculate_metrics <- function(observed, predicted) {
  r2 <- cor(observed, predicted, use = "complete.obs")^2
  rmse <- sqrt(mean((observed - predicted)^2, na.rm = TRUE))
  mae <- mean(abs(observed - predicted), na.rm = TRUE)
  bias <- mean(predicted - observed, na.rm = TRUE)

  return(data.frame(r2 = r2, rmse = rmse, mae = mae, bias = bias))
}


#' Fit a random forest model
#'
#' @param train_data training data frame containing predictor variables and response
#' @param test_data test data frame for making predictions
#' @param predictor_vars character vector of predictor variable names
#' @param hyperparams optional list of hyperparameters (ntree, mtry, nodesize). if null, uses defaults
#' @return list containing the fitted model and predictions on test data
fit_rf <- function(train_data, test_data, predictor_vars, hyperparams = NULL) {
  formula_str <- paste("median_carbon_density ~", paste(predictor_vars, collapse = " + "))
  formula_obj <- as.formula(formula_str)

  # use provided hyperparameters or defaults
  ntree <- if (!is.null(hyperparams$ntree)) hyperparams$ntree else 500
  mtry <- if (!is.null(hyperparams$mtry)) hyperparams$mtry else floor(sqrt(length(predictor_vars)))
  nodesize <- if (!is.null(hyperparams$nodesize)) hyperparams$nodesize else NULL

  model <- randomForest(formula_obj,
    data = train_data,
    ntree = ntree,
    mtry = mtry,
    nodesize = nodesize,
    importance = TRUE
  )

  predictions <- predict(model, newdata = test_data)

  return(list(model = model, predictions = predictions))
}

#' Fit an xgboost (gradient boosted trees) model
#'
#' @param train_data training data frame containing predictor variables and response
#' @param test_data test data frame for making predictions
#' @param predictor_vars character vector of predictor variable names
#' @param hyperparams optional list of hyperparameters (nrounds, max_depth, learning_rate, subsample, colsample_bytree). if null, uses defaults
#' @return list containing the fitted model and predictions on test data
fit_xgboost <- function(train_data, test_data, predictor_vars, hyperparams = NULL) {
  # prepare data matrices
  X_train <- as.matrix(train_data[, predictor_vars])
  y_train <- train_data$median_carbon_density
  X_test <- as.matrix(test_data[, predictor_vars])

  # use provided hyperparameters or defaults
  nrounds <- if (!is.null(hyperparams$nrounds)) hyperparams$nrounds else 100
  max_depth <- if (!is.null(hyperparams$max_depth)) hyperparams$max_depth else 6
  learning_rate <- if (!is.null(hyperparams$learning_rate)) hyperparams$learning_rate else 0.3
  subsample <- if (!is.null(hyperparams$subsample)) hyperparams$subsample else 0.8
  colsample_bytree <- if (!is.null(hyperparams$colsample_bytree)) hyperparams$colsample_bytree else 0.8
  nthreads <- if (!is.null(hyperparams$nthreads) && !is.na(hyperparams$nthreads)) as.integer(hyperparams$nthreads)[1L] else 1L

  # fit model
  model <- xgboost(
    x = X_train,
    y = y_train,
    nrounds = nrounds,
    max_depth = max_depth,
    learning_rate = learning_rate,
    subsample = subsample,
    colsample_bytree = colsample_bytree,
    nthread = nthreads,
    objective = "reg:squarederror"
  )

  predictions <- predict(model, newdata = X_test)

  return(list(model = model, predictions = predictions))
}

#' Fit a support vector machine (svm) model with automatic parameter tuning
#'
#' @param train_data training data frame containing predictor variables and response
#' @param test_data test data frame for making predictions
#' @param predictor_vars character vector of predictor variable names
#' @return list containing the fitted model and predictions on test data
fit_svm <- function(train_data, test_data, predictor_vars) {
  formula_str <- paste("median_carbon_density ~", paste(predictor_vars, collapse = " + "))
  formula_obj <- as.formula(formula_str)

  # tune svm parameters (cost, gamma, and kernel type) using grid search
  # tune.svm doesn't accept a vector for kernel, so we need to tune each kernel separately
  # and then select the best overall model
  kernels_to_try <- c("linear", "polynomial", "sigmoid")
  best_result <- NULL
  best_performance <- Inf

  for (kernel_type in kernels_to_try) {
    # tune parameters for this kernel type
    # note: linear kernel doesn't use gamma, so we skip it for non-linear kernels
    # polynomial and sigmoid kernels use both cost and gamma
    tune_result <- tryCatch(
      {
        tune.svm(formula_obj,
          data = train_data,
          gamma = 10^(-3:1),
          cost = 10^(-1:2),
          kernel = kernel_type,
          tunecontrol = tune.control(sampling = "cross", cross = 3)
        )
      },
      error = function(e) {
        cat("    Warning: Failed to tune", kernel_type, "kernel:", e$message, "\n")
        return(NULL)
      }
    )

    if (!is.null(tune_result)) {
      # extract performance (error) - lower is better
      current_performance <- tune_result$best.performance
      if (current_performance < best_performance) {
        best_performance <- current_performance
        best_result <- tune_result
      }
    }
  }

  # if no kernel worked, fall back to default svm
  if (is.null(best_result)) {
    cat("    Warning: All kernel tuning failed, using default svm\n")
    model <- svm(formula_obj, data = train_data)
  } else {
    model <- best_result$best.model
  }

  predictions <- predict(model, newdata = test_data)

  return(list(model = model, predictions = predictions))
}

#' Fit a neural network model
#'
#' @param train_data training data frame containing predictor variables and response
#' @param test_data test data frame for making predictions
#' @param predictor_vars character vector of predictor variable names
#' @param hyperparams optional list of hyperparameters (hidden, learningrate). if null, uses defaults
#' @return list containing the fitted model and predictions on test data (unscaled)
fit_nn <- function(train_data, test_data, predictor_vars, hyperparams = NULL) {
  # neuralnet requires all numeric data and works better with scaled data
  # ensure data is a data.frame (not tibble) to avoid namespace conflicts
  train_data <- as.data.frame(train_data)
  test_data <- as.data.frame(test_data)

  # prepare training data: select only numeric predictors and response
  # use explicit dplyr namespace to avoid conflicts
  train_numeric <- train_data %>%
    dplyr::select(dplyr::all_of(predictor_vars), median_carbon_density) %>%
    dplyr::mutate_all(as.numeric) %>%
    as.data.frame()

  # remove any rows with missing values
  train_numeric <- train_numeric[complete.cases(train_numeric), ]

  # scale the data (neural networks work better with scaled inputs)
  # store scaling parameters to apply to test data
  scale_params <- list(
    means = colMeans(train_numeric[, predictor_vars, drop = FALSE]),
    sds = apply(train_numeric[, predictor_vars, drop = FALSE], 2, sd)
  )

  # scale training data
  train_scaled <- train_numeric
  for (var in predictor_vars) {
    if (scale_params$sds[var] > 0) { # avoid division by zero
      train_scaled[[var]] <- (train_numeric[[var]] - scale_params$means[var]) / scale_params$sds[var]
    } else {
      train_scaled[[var]] <- 0 # if no variance, set to 0
    }
  }

  # scale response variable separately
  response_mean <- mean(train_numeric$median_carbon_density)
  response_sd <- sd(train_numeric$median_carbon_density)
  train_scaled$median_carbon_density <- (train_numeric$median_carbon_density - response_mean) / response_sd

  # prepare test data
  test_numeric <- test_data %>%
    dplyr::select(dplyr::all_of(predictor_vars), median_carbon_density) %>%
    dplyr::mutate_all(as.numeric) %>%
    as.data.frame()

  # scale test data using training scaling parameters
  test_scaled <- test_numeric
  for (var in predictor_vars) {
    if (scale_params$sds[var] > 0) {
      test_scaled[[var]] <- (test_numeric[[var]] - scale_params$means[var]) / scale_params$sds[var]
    } else {
      test_scaled[[var]] <- 0
    }
  }

  # create formula
  formula_str <- paste("median_carbon_density ~", paste(predictor_vars, collapse = " + "))
  formula_obj <- as.formula(formula_str)

  # use provided hyperparameters or defaults
  hidden <- if (!is.null(hyperparams$hidden)) hyperparams$hidden else c(10, 5)
  learningrate <- if (!is.null(hyperparams$learningrate)) hyperparams$learningrate else NULL

  # fit neural network
  # linear.output = true for regression (not classification)
  # act.fct = "logistic" for hidden layers, but output is linear for regression
  model_args <- list(
    formula = formula_obj,
    data = train_scaled,
    hidden = hidden,
    act.fct = "logistic",
    linear.output = TRUE,
    threshold = 0.01,
    stepmax = 1e+05
  )
  if (!is.null(learningrate)) model_args$learningrate <- learningrate

  model <- do.call(neuralnet, model_args)

  # predict on scaled test data
  # neuralnet predict returns a matrix, extract first column
  predictions_scaled <- predict(model, newdata = test_scaled[, predictor_vars, drop = FALSE])

  # convert predictions_scaled from matrix to vector if needed
  if (is.matrix(predictions_scaled)) {
    predictions_scaled <- predictions_scaled[, 1]
  }

  # unscale predictions back to original scale
  predictions <- predictions_scaled * response_sd + response_mean

  # store scaling parameters in model object for later use
  model$scale_params <- scale_params
  model$response_mean <- response_mean
  model$response_sd <- response_sd

  return(list(model = model, predictions = predictions))
}

fit_gpr <- function(train_data, test_data, predictor_vars) {
  formula_str <- paste("median_carbon_density ~", paste(predictor_vars, collapse = " + "))
  formula_obj <- as.formula(formula_str)

  model <- gpkm(formula_obj, data = train_data, kernel = "matern52")
  XX <- as.matrix(test_data[, predictor_vars, drop = FALSE])
  predictions <- predict(model, XX = XX)
  return(list(model = model, predictions = predictions))
}

# ================================ SPATIAL INTERPOLATION (KRIGING) ================================

#' Project points to a metric CRS for variogram/kriging (UTM zone from centroid)
#' @param sf_pts sf object with point geometry (WGS84)
#' @return sf object in UTM (metres)
project_for_kriging <- function(sf_pts) {
  cen <- sf::st_union(sf::st_geometry(sf_pts)) %>% sf::st_centroid()
  xy <- sf::st_coordinates(cen)
  zone <- floor((xy[1, "X"] + 180) / 6) + 1
  hem <- if (xy[1, "Y"] >= 0) " +north" else " +south"
  crs_utm <- paste0("+proj=utm +zone=", zone, hem, " +datum=WGS84 +units=m")
  sf::st_transform(sf_pts, crs_utm)
}

#' Fit ordinary kriging (spatial correlation only) and predict at test locations
#' Uses gstat; aggregates to unique locations and local kriging (nmax) to avoid singular matrix.
#' @param train_data data frame with median_carbon_density, longitude, latitude
#' @param test_data data frame with longitude, latitude
#' @param nmax max nearest neighbours per prediction (default 25); avoids singular covariance.
#' @return list(model = gstat object or NULL, predictions = vector)
fit_ok <- function(train_data, test_data, predictor_vars = NULL, nmax = 25) {
  if (!requireNamespace("gstat", quietly = TRUE)) {
    return(list(model = NULL, predictions = rep(NA_real_, nrow(test_data))))
  }
  # Aggregate to unique locations (duplicate coords cause singular covariance in krige)
  train_agg <- train_data %>%
    dplyr::group_by(longitude, latitude) %>%
    dplyr::summarise(median_carbon_density = mean(median_carbon_density, na.rm = TRUE), .groups = "drop") %>%
    as.data.frame()
  if (nrow(train_agg) < 4) {
    return(list(model = NULL, predictions = rep(NA_real_, nrow(test_data))))
  }
  train_sf <- sf::st_as_sf(train_agg, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)
  test_sf <- sf::st_as_sf(test_data, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)
  train_sf <- project_for_kriging(train_sf)
  test_sf <- project_for_kriging(test_sf)
  # Variogram: use cutoff up to ~1/3 of max distance so we have enough bins
  d_obs <- sf::st_coordinates(train_sf)
  max_d <- max(sqrt(outer(d_obs[, 1], d_obs[, 1], "-")^2 + outer(d_obs[, 2], d_obs[, 2], "-")^2), na.rm = TRUE)
  cutoff <- min(500e3, max(1.5 * max_d, 10e3))
  v <- tryCatch(
    gstat::variogram(median_carbon_density ~ 1, data = train_sf, cutoff = cutoff),
    error = function(e) NULL
  )
  if (is.null(v) || nrow(v) < 3) {
    return(list(model = NULL, predictions = rep(NA_real_, nrow(test_data))))
  }
  v0 <- var(train_agg$median_carbon_density, na.rm = TRUE)
  if (!is.finite(v0) || v0 <= 0) v0 <- 0.1
  fit <- tryCatch(
    gstat::fit.variogram(v, gstat::vgm(NA, "Sph", NA, NA), fit.method = 7),
    error = function(e) {
      tryCatch(
        gstat::fit.variogram(v, gstat::vgm(v0, "Sph", 50e3, 0.1), fit.method = 7),
        error = function(e2) {
          tryCatch(
            gstat::fit.variogram(v, gstat::vgm(v0, "Exp", 50e3, 0.1), fit.method = 7),
            error = function(e3) NULL
          )
        }
      )
    }
  )
  if (is.null(fit) || !is.data.frame(fit) || nrow(fit) < 2 ||
    (length(fit$range) >= 2 && fit$range[2] <= 0)) {
    return(list(model = NULL, predictions = rep(NA_real_, nrow(test_data))))
  }
  # Minimum nugget for numerical stability (avoids singular kriging matrix)
  sill_tot <- sum(fit$psill, na.rm = TRUE)
  if (sill_tot > 0 && fit$psill[1] < 1e-6 * sill_tot) fit$psill[1] <- 1e-6 * sill_tot
  pred <- tryCatch(
    gstat::krige(median_carbon_density ~ 1, locations = train_sf, newdata = test_sf, model = fit, nmax = nmax),
    error = function(e) NULL
  )
  if (is.null(pred)) {
    return(list(model = list(variogram = fit), predictions = rep(NA_real_, nrow(test_data))))
  }
  list(model = list(variogram = fit), predictions = as.numeric(pred$var1.pred))
}

#' Fit universal kriging with environmental covariates as drift and predict at test locations
#' Aggregates to unique locations and uses local kriging (nmax) to avoid singular matrix.
#' @param train_data data frame with median_carbon_density, longitude, latitude, and predictor_vars
#' @param test_data data frame with longitude, latitude, and predictor_vars
#' @param predictor_vars character vector of covariate names (drift terms)
#' @param nmax max nearest neighbours per prediction (default 25).
#' @return list(model = gstat object or NULL, predictions = vector)
fit_uk <- function(train_data, test_data, predictor_vars, nmax = 25) {
  if (!requireNamespace("gstat", quietly = TRUE) || length(predictor_vars) == 0) {
    return(list(model = NULL, predictions = rep(NA_real_, nrow(test_data))))
  }
  pv <- intersect(predictor_vars, names(train_data))
  pv <- pv[pv %in% names(test_data)]
  if (length(pv) == 0) {
    return(list(model = NULL, predictions = rep(NA_real_, nrow(test_data))))
  }
  # Aggregate to unique locations (duplicate coords cause singular covariance)
  train_agg <- train_data %>%
    dplyr::group_by(longitude, latitude) %>%
    dplyr::summarise(
      median_carbon_density = mean(median_carbon_density, na.rm = TRUE),
      dplyr::across(dplyr::all_of(pv), ~ if (is.numeric(.)) mean(., na.rm = TRUE) else .[1]),
      .groups = "drop"
    ) %>%
    as.data.frame()
  if (nrow(train_agg) < 4) {
    return(list(model = NULL, predictions = rep(NA_real_, nrow(test_data))))
  }
  train_sf <- sf::st_as_sf(train_agg, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)
  test_sf <- sf::st_as_sf(test_data, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)
  train_sf <- project_for_kriging(train_sf)
  test_sf <- project_for_kriging(test_sf)
  form <- as.formula(paste("median_carbon_density ~", paste(pv, collapse = " + ")))
  max_d <- max(sqrt(outer(sf::st_coordinates(train_sf)[, 1], sf::st_coordinates(train_sf)[, 1], "-")^2 +
    outer(sf::st_coordinates(train_sf)[, 2], sf::st_coordinates(train_sf)[, 2], "-")^2), na.rm = TRUE)
  cutoff <- min(500e3, max(1.5 * max_d, 10e3))
  v <- tryCatch(
    gstat::variogram(form, data = train_sf, cutoff = cutoff),
    error = function(e) NULL
  )
  if (is.null(v) || nrow(v) < 3) {
    return(list(model = NULL, predictions = rep(NA_real_, nrow(test_data))))
  }
  fit <- tryCatch(
    gstat::fit.variogram(v, gstat::vgm(NA, "Sph", NA, NA), fit.method = 7),
    error = function(e) {
      tryCatch(
        gstat::fit.variogram(v, gstat::vgm(0.1, "Sph", 50e3, 0.01), fit.method = 7),
        error = function(e2) {
          tryCatch(
            gstat::fit.variogram(v, gstat::vgm(0.1, "Exp", 50e3, 0.01), fit.method = 7),
            error = function(e3) NULL
          )
        }
      )
    }
  )
  if (is.null(fit) || !is.data.frame(fit) || nrow(fit) < 2 ||
    (length(fit$range) >= 2 && fit$range[2] <= 0)) {
    return(list(model = NULL, predictions = rep(NA_real_, nrow(test_data))))
  }
  sill_tot <- sum(fit$psill, na.rm = TRUE)
  if (sill_tot > 0 && fit$psill[1] < 1e-6 * sill_tot) fit$psill[1] <- 1e-6 * sill_tot
  pred <- tryCatch(
    gstat::krige(form, locations = train_sf, newdata = test_sf, model = fit, nmax = nmax),
    error = function(e) NULL
  )
  if (is.null(pred)) {
    return(list(model = list(variogram = fit), predictions = rep(NA_real_, nrow(test_data))))
  }
  list(model = list(variogram = fit), predictions = as.numeric(pred$var1.pred))
}

#' Fit INLA SPDE (Matérn field) on train data and predict at test locations (for CV).
#' Requires INLA; returns NA predictions if INLA not installed or fit fails.
#' @param train_data data frame with median_carbon_density, longitude, latitude, and predictor_vars
#' @param test_data data frame with longitude, latitude, and predictor_vars
#' @param predictor_vars character vector of covariate names (fixed effects)
#' @return list(model = INLA fit or NULL, predictions = vector)
fit_inla_cv <- function(train_data, test_data, predictor_vars) {
  if (!requireNamespace("INLA", quietly = TRUE)) {
    return(list(model = NULL, predictions = rep(NA_real_, nrow(test_data))))
  }
  if (nrow(train_data) < 5 || nrow(test_data) < 1) {
    return(list(model = NULL, predictions = rep(NA_real_, nrow(test_data))))
  }
  coords <- c("longitude", "latitude")
  value_var <- "median_carbon_density"
  loc_obs <- as.matrix(train_data[, coords])
  loc_pred <- as.matrix(test_data[, coords])
  y <- as.numeric(train_data[[value_var]])
  pv <- intersect(if (length(predictor_vars) == 0) character(0) else predictor_vars, names(train_data))
  pv <- pv[pv %in% names(test_data)]
  bnd <- try(INLA::inla.nonconvex.hull(loc_obs, convex = -0.05), silent = TRUE)
  if (inherits(bnd, "try-error")) {
    return(list(model = NULL, predictions = rep(NA_real_, nrow(test_data))))
  }
  mesh <- try(INLA::inla.mesh.2d(boundary = bnd, max.edge = c(0.5, 2), cutoff = 0.2), silent = TRUE)
  if (inherits(mesh, "try-error")) {
    return(list(model = NULL, predictions = rep(NA_real_, nrow(test_data))))
  }
  spde <- INLA::inla.spde2.matern(mesh, alpha = 2)
  A_obs <- INLA::inla.spde.make.A(mesh, loc = loc_obs)
  A_pred <- INLA::inla.spde.make.A(mesh, loc = loc_pred)
  field_idx <- INLA::inla.spde.make.index("field", n.spde = spde$n.spde)
  n_obs <- length(y)
  n_pred <- nrow(loc_pred)
  fixed_est <- list(intercept = rep(1, n_obs))
  fixed_pred <- list(intercept = rep(1, n_pred))
  for (v in pv) {
    fixed_est[[v]] <- as.numeric(train_data[[v]])
    fixed_pred[[v]] <- as.numeric(test_data[[v]])
  }
  s_est <- INLA::inla.stack(data = list(y = y), A = list(A_obs, 1), effects = list(field_idx, fixed_est), tag = "est")
  s_pred <- INLA::inla.stack(data = list(y = NA), A = list(A_pred, 1), effects = list(field_idx, fixed_pred), tag = "pred")
  stack <- INLA::inla.stack(s_est, s_pred)
  formula_str <- "y ~ -1 + intercept"
  if (length(pv) > 0) formula_str <- paste(formula_str, "+", paste(pv, collapse = " + "))
  formula_str <- paste(formula_str, "+ f(field, model = spde)")
  formula <- as.formula(formula_str)
  fit <- try(INLA::inla(formula, data = INLA::inla.stack.data(stack), control.predictor = list(A = INLA::inla.stack.A(stack), compute = TRUE)), silent = TRUE)
  if (inherits(fit, "try-error")) {
    return(list(model = NULL, predictions = rep(NA_real_, nrow(test_data))))
  }
  idx_pred <- INLA::inla.stack.index(stack, "pred")$data
  pred_vec <- as.numeric(fit$summary.fitted.values$mean[idx_pred])
  list(model = fit, predictions = pred_vec)
}

#' Fit GAM with spatial smooth + covariates (same specification as sensitivity analysis)
#' @param train_data data frame with median_carbon_density, longitude, latitude, and predictor_vars
#' @param test_data data frame with longitude, latitude, and predictor_vars
#' @param predictor_vars character vector of covariate names (linear terms)
#' @param k_spatial basis dimension for spatial smooth (default 80)
#' @return list(model = gam object or NULL, predictions = vector)
fit_gam <- function(train_data, test_data, predictor_vars, k_spatial = 80) {
  if (!requireNamespace("mgcv", quietly = TRUE)) {
    return(list(model = NULL, predictions = rep(NA_real_, nrow(test_data))))
  }
  y_train <- train_data$median_carbon_density
  family_used <- if (any(y_train <= 0, na.rm = TRUE)) {
    gaussian(link = "identity")
  } else {
    Gamma(link = "log")
  }
  covar_terms <- paste(predictor_vars, collapse = " + ")
  spatial_term <- paste0("s(longitude, latitude, k = ", k_spatial, ")")
  form <- as.formula(paste("median_carbon_density ~", spatial_term, "+", covar_terms))
  fit <- try(
    mgcv::gam(form, data = train_data, family = family_used, method = "REML"),
    silent = TRUE
  )
  if (inherits(fit, "try-error")) {
    return(list(model = NULL, predictions = rep(NA_real_, nrow(test_data))))
  }
  pred <- try(
    as.numeric(mgcv::predict.gam(fit, newdata = test_data, type = "response")),
    silent = TRUE
  )
  if (inherits(pred, "try-error") || any(is.na(pred))) {
    pred <- rep(mean(y_train, na.rm = TRUE), nrow(test_data))
  }
  list(model = fit, predictions = pred)
}

#' Create distance-buffered CV folds: test sets chosen so each test point is > buffer_m from training
#' Tests "Can I fill in gaps between my samples?"
#' @param core_sf sf object with point geometry (will use st_distance, so CRS in metres or WGS84)
#' @param n_folds number of folds (test fraction ≈ 1/n_folds)
#' @param buffer_m minimum distance (m) between test and training points
#' @return list of length n_folds; folds[[k]] = vector of row indices for test set for fold k
create_distance_buffered_folds <- function(core_sf, n_folds, buffer_m = 1000) {
  n <- nrow(core_sf)
  target_per_fold <- max(1L, floor(n / n_folds))
  assigned <- rep(FALSE, n)
  folds <- vector("list", n_folds)
  order_idx <- sample.int(n)
  for (f in seq_len(n_folds)) {
    while (length(folds[[f]]) < target_per_fold) {
      cand <- which(!assigned)
      if (length(cand) == 0) break
      if (length(folds[[f]]) == 0) {
        pick <- order_idx[order_idx %in% cand][1]
        folds[[f]] <- c(folds[[f]], pick)
        assigned[pick] <- TRUE
        next
      }
      dist_to_test <- sf::st_distance(core_sf[cand, ], core_sf[folds[[f]], ])
      min_dist <- apply(dist_to_test, 1, min)
      valid <- cand[as.numeric(min_dist) > buffer_m]
      if (length(valid) == 0) break
      pick <- order_idx[order_idx %in% valid][1]
      if (is.na(pick)) break
      folds[[f]] <- c(folds[[f]], pick)
      assigned[pick] <- TRUE
    }
  }
  remaining <- which(!assigned)
  for (i in remaining) {
    d <- sapply(seq_len(n_folds), function(k) {
      if (length(folds[[k]]) == 0) {
        return(Inf)
      }
      min(as.numeric(sf::st_distance(core_sf[i, ], core_sf[folds[[k]], ])))
    })
    f <- which.max(d)
    folds[[f]] <- c(folds[[f]], i)
    assigned[i] <- TRUE
  }
  folds
}


#' Clean database in preparation for modelling. Ensures that the carbon density and sediment_mean_depth are both above zero and not NA.
#'
#' @param dat data frame containing the data (including the carbon_density_g_c_cm3 and sediment_mean_depth_cm columns)
#' @return data frame containing the cleaned data
clean_data <- function(dat) {
  dat_clean <- dat %>%
    filter(!is.na(carbon_density_g_c_cm3) & !is.na(sediment_mean_depth_cm) &
      carbon_density_g_c_cm3 > 0 & sediment_mean_depth_cm >= 0)
  return(dat_clean)
}

#' Fit an exponential carbon decay model
#'
#' @param dat data frame containing the data
#' @return list containing the fitted model, cleaned data, prediction data frame, and residuals standard deviation
fit_exponential_carbon_decay <- function(dat) {
  # clean data
  dat_clean <- clean_data(dat) # belt and braces (data passed should already be cleaned)
  if (nrow(dat_clean) < 3) {
    warning("Insufficient data points after cleaning for NLS fitting.")
    return(NULL)
  }

  # Fit nls model: carbon_density = a * exp(b * depth)
  nls_model <- tryCatch(
    {
      nls(carbon_density_g_c_cm3 ~ a * exp(b * sediment_mean_depth_cm),
        data = dat_clean,
        start = list(
          a = max(dat_clean$carbon_density_g_c_cm3, na.rm = TRUE),
          b = -0.01
        )
      )
    },
    error = function(e) {
      cat("nls fitting failed, trying alternative starting values\n")
      nls(carbon_density_g_c_cm3 ~ a * exp(b * sediment_mean_depth_cm),
        data = dat_clean,
        start = list(
          a = mean(dat_clean$carbon_density_g_c_cm3, na.rm = TRUE),
          b = -0.001
        ),
        control = nls.control(maxiter = 500, warnOnly = TRUE)
      )
    }
  )

  # Create prediction data frame
  depth_seq <- seq(min(dat_clean$sediment_mean_depth_cm, na.rm = TRUE),
    max(dat_clean$sediment_mean_depth_cm, na.rm = TRUE),
    length.out = 200
  )
  pred_df <- data.frame(sediment_mean_depth_cm = depth_seq)
  pred_df$carbon_density_fit <- predict(nls_model, newdata = pred_df)

  # Calculate approximate confidence intervals based on residual standard error
  residuals_sd <- sd(residuals(nls_model), na.rm = TRUE)
  pred_df$carbon_density_se <- residuals_sd
  pred_df$carbon_density_lower <- pred_df$carbon_density_fit - 1.96 * pred_df$carbon_density_se
  pred_df$carbon_density_upper <- pred_df$carbon_density_fit + 1.96 * pred_df$carbon_density_se

  return(list(
    nls_model = nls_model,
    dat_clean = dat_clean,
    pred_df = pred_df,
    residuals_sd = residuals_sd
  ))
}

#' Fit an exponential carbon decay model per species
#'
#' @param dat_clean data frame containing the data
#' @param spec character string specifying the species
#' @return list containing the fitted model and prediction data frame
fit_species <- function(dat, spec) {
  species_dat <- dat %>% filter(seagrass_species == spec)
  if (nrow(species_dat) < 3) {
    return(NULL)
  }
  fit <- tryCatch(
    {
      nls(carbon_density_g_c_cm3 ~ a * exp(b * sediment_mean_depth_cm),
        data = species_dat,
        start = list(
          a = max(species_dat$carbon_density_g_c_cm3, na.rm = TRUE),
          b = -0.01
        )
      )
    },
    error = function(e) {
      nls(carbon_density_g_c_cm3 ~ a * exp(b * sediment_mean_depth_cm),
        data = species_dat,
        start = list(
          a = mean(species_dat$carbon_density_g_c_cm3, na.rm = TRUE),
          b = -0.001
        ),
        control = nls.control(maxiter = 500, warnOnly = TRUE)
      )
    }
  )
  # predictions
  depth_seq <- seq(
    min(species_dat$sediment_mean_depth_cm, na.rm = TRUE),
    max(species_dat$sediment_mean_depth_cm, na.rm = TRUE),
    length.out = 200
  )
  pred_df <- data.frame(
    seagrass_species = spec,
    sediment_mean_depth_cm = depth_seq
  )
  pred_df$carbon_density_fit <- predict(fit, newdata = pred_df)
  # approx CI
  residuals_sd <- sd(residuals(fit), na.rm = TRUE)
  pred_df$carbon_density_se <- residuals_sd
  pred_df$carbon_density_lower <- pred_df$carbon_density_fit - 1.96 * residuals_sd
  pred_df$carbon_density_upper <- pred_df$carbon_density_fit + 1.96 * residuals_sd
  list(fit = fit, pred_df = pred_df)
}


#' Compare different depth profile models
#'
#' Tests multiple functional forms: exponential, power law, linear, GAM
#' with and without species interactions. Also tests normalization by core aggregate.
#'
#' @param dat_clean data frame containing the data
#' @param normalise_by_core logical, whether to normalize carbon density by core aggregate
#' @param aggregate_fun function to calculate per-core aggregate (default: median of top 5 cm)
#' @return a list containing fitted models and comparison table
compare_depth_profile_models <- function(dat_clean,
                                         normalise_by_core = FALSE,
                                         aggregate_fun = calculate_core_aggregate) {
  # prepare data
  dat_sub <- dat_clean %>%
    filter(!is.na(carbon_density_g_c_cm3) & !is.na(sediment_mean_depth_cm) &
      !is.na(seagrass_species))

  # normalize if requested
  if (normalise_by_core) {
    core_aggregates <- aggregate_fun(dat_sub)
    dat_sub <- dat_sub %>%
      left_join(core_aggregates, by = "random_core_variable") %>%
      mutate(carbon_density_g_c_cm3 = carbon_density_g_c_cm3 / aggregate_carbon_density) %>%
      select(-aggregate_carbon_density)
  }

  species_fac <- as.factor(dat_sub$seagrass_species)
  depth <- dat_sub$sediment_mean_depth_cm
  carbon <- dat_sub$carbon_density_g_c_cm3

  models <- list()
  model_comparison <- data.frame()

  print("Exponential decay: c = a * exp(b * depth) - no interaction")
  # 1. exponential decay: c = a * exp(b * depth) - no interaction
  models$exp_null <- tryCatch(
    {
      nls(carbon ~ a * exp(b * depth),
        data = dat_sub,
        start = list(a = max(carbon, na.rm = TRUE), b = -0.01),
        control = nls.control(maxiter = 500, warnOnly = TRUE)
      )
    },
    error = function(e) NULL
  )

  if (!is.null(models$exp_null)) {
    aic_val <- AIC(models$exp_null)
    bic_val <- BIC(models$exp_null)
    model_comparison <- rbind(
      model_comparison,
      data.frame(
        model = "exp_null", aic = aic_val, bic = bic_val,
        n_params = 2, has_species_interaction = FALSE
      )
    )
  }

  print("Exponential decay: c = a * exp(b * depth) - species interaction")
  # 2. exponential decay with species interaction
  species_list <- unique(species_fac)
  exp_species_models <- list()
  aic_total <- 0
  n_params_total <- 0
  n_successful <- 0

  for (sp in species_list) {
    idx <- species_fac == sp
    if (sum(idx) >= 5) {
      carbon_subset <- carbon[idx]
      depth_subset <- depth[idx]
      valid_data <- !is.na(carbon_subset) & !is.na(depth_subset) & carbon_subset > 0

      if (sum(valid_data) >= 5) {
        m <- tryCatch(
          {
            nls(carbon_subset[valid_data] ~ a * exp(b * depth_subset[valid_data]),
              start = list(a = max(carbon_subset[valid_data], na.rm = TRUE), b = -0.01),
              control = nls.control(maxiter = 500, warnOnly = TRUE)
            )
          },
          error = function(e) NULL
        )

        if (!is.null(m)) {
          exp_species_models[[as.character(sp)]] <- m
          aic_total <- aic_total + AIC(m)
          n_params_total <- n_params_total + 2
          n_successful <- n_successful + 1
        }
      }
    }
  }

  if (n_successful > 0) {
    models$exp_species <- exp_species_models
    # approximate bic: aic + k * log(n) where k is number of parameters
    n_obs <- sum(sapply(exp_species_models, function(m) length(residuals(m))))
    bic_total <- aic_total + n_params_total * log(n_obs)
    model_comparison <- rbind(
      model_comparison,
      data.frame(
        model = "exp_species", aic = aic_total, bic = bic_total,
        n_params = n_params_total, has_species_interaction = TRUE
      )
    )
  }

  # 3. power law: c = a * depth^b - no interaction
  print("Power law: c = a * depth^b - species interaction")
  models$power_null <- tryCatch(
    {
      nls(carbon ~ a * depth^b,
        data = dat_sub,
        start = list(a = max(carbon, na.rm = TRUE), b = -0.5),
        control = nls.control(maxiter = 500, warnOnly = TRUE)
      )
    },
    error = function(e) NULL
  )

  if (!is.null(models$power_null)) {
    aic_val <- AIC(models$power_null)
    bic_val <- BIC(models$power_null)
    model_comparison <- rbind(
      model_comparison,
      data.frame(
        model = "power_null", aic = aic_val, bic = bic_val,
        n_params = 2, has_species_interaction = FALSE
      )
    )
  }

  # 4. power law with species interaction
  print("Power law: c = a * depth^b - species interaction")
  power_species_models <- list()
  aic_total <- 0
  n_params_total <- 0
  n_successful <- 0

  for (sp in species_list) {
    idx <- species_fac == sp
    if (sum(idx) >= 5) {
      carbon_subset <- carbon[idx]
      depth_subset <- depth[idx]
      valid_data <- !is.na(carbon_subset) & !is.na(depth_subset) & carbon_subset > 0

      if (sum(valid_data) >= 5) {
        m <- tryCatch(
          {
            nls(carbon_subset[valid_data] ~ a * depth_subset[valid_data]^b,
              start = list(a = max(carbon_subset[valid_data], na.rm = TRUE), b = -0.5),
              control = nls.control(maxiter = 500, warnOnly = TRUE)
            )
          },
          error = function(e) NULL
        )

        if (!is.null(m)) {
          power_species_models[[as.character(sp)]] <- m
          aic_total <- aic_total + AIC(m)
          n_params_total <- n_params_total + 2
          n_successful <- n_successful + 1
        }
      }
    }
  }

  if (n_successful > 0) {
    models$power_species <- power_species_models
    n_obs <- sum(sapply(power_species_models, function(m) length(residuals(m))))
    bic_total <- aic_total + n_params_total * log(n_obs)
    model_comparison <- rbind(
      model_comparison,
      data.frame(
        model = "power_species", aic = aic_total, bic = bic_total,
        n_params = n_params_total, has_species_interaction = TRUE
      )
    )
  }

  # 5. linear: c = a + b * depth - no interaction
  print("Linear: c = a + b * depth - no interaction")
  models$linear_null <- tryCatch(
    {
      lm(carbon ~ depth, data = dat_sub)
    },
    error = function(e) NULL
  )

  if (!is.null(models$linear_null)) {
    aic_val <- AIC(models$linear_null)
    bic_val <- BIC(models$linear_null)
    model_comparison <- rbind(
      model_comparison,
      data.frame(
        model = "linear_null", aic = aic_val, bic = bic_val,
        n_params = 2, has_species_interaction = FALSE
      )
    )
  }

  # 6. linear with species interaction
  print("Linear: c = a + b * depth - species interaction")
  models$linear_species <- tryCatch(
    {
      lm(carbon ~ depth * species_fac, data = dat_sub)
    },
    error = function(e) NULL
  )

  if (!is.null(models$linear_species)) {
    aic_val <- AIC(models$linear_species)
    bic_val <- BIC(models$linear_species)
    n_params <- length(coef(models$linear_species))
    model_comparison <- rbind(
      model_comparison,
      data.frame(
        model = "linear_species", aic = aic_val, bic = bic_val,
        n_params = n_params, has_species_interaction = TRUE
      )
    )
  }

  # 7. GAM without species interaction
  print("GAM: c = a * exp(b * depth) - no interaction")
  models$gam_null <- tryCatch(
    {
      gam(carbon ~ s(depth, k = 5), data = dat_sub, method = "REML")
    },
    error = function(e) NULL
  )

  if (!is.null(models$gam_null)) {
    aic_val <- AIC(models$gam_null)
    bic_val <- BIC(models$gam_null)
    n_params <- length(coef(models$gam_null))
    model_comparison <- rbind(
      model_comparison,
      data.frame(
        model = "gam_null", aic = aic_val, bic = bic_val,
        n_params = n_params, has_species_interaction = FALSE
      )
    )
  }

  # 8. GAM with species interaction (by smooth)
  print("GAM: c = a * exp(b * depth) - species interaction (by smooth)")
  models$gam_species_by <- tryCatch(
    {
      gam(carbon ~ s(depth, by = species_fac, k = 5) + species_fac,
        data = dat_sub, method = "REML"
      )
    },
    error = function(e) NULL
  )

  if (!is.null(models$gam_species_by)) {
    aic_val <- AIC(models$gam_species_by)
    bic_val <- BIC(models$gam_species_by)
    n_params <- length(coef(models$gam_species_by))
    model_comparison <- rbind(
      model_comparison,
      data.frame(
        model = "gam_species_by", aic = aic_val, bic = bic_val,
        n_params = n_params, has_species_interaction = TRUE
      )
    )
  }

  # 9. GAM with tensor product (allows interaction)
  print("GAM: c = a * exp(b * depth) - species interaction (by tensor product)")
  models$gam_tensor <- tryCatch(
    {
      gam(carbon ~ te(depth, species_fac, k = c(5, 3)), data = dat_sub, method = "REML")
    },
    error = function(e) NULL
  )

  if (!is.null(models$gam_tensor)) {
    aic_val <- AIC(models$gam_tensor)
    bic_val <- BIC(models$gam_tensor)
    n_params <- length(coef(models$gam_tensor))
    model_comparison <- rbind(
      model_comparison,
      data.frame(
        model = "gam_tensor", aic = aic_val, bic = bic_val,
        n_params = n_params, has_species_interaction = TRUE
      )
    )
  }

  # 10. GAM with unique core random effect
  # print("GAM: c = a * exp(b * depth) - no interaction (with unique core random effect)")
  # models$gam_core_re <- tryCatch({
  #   gam(carbon ~ s(depth) + s(random_core_variable, bs = "re"), data = dat_sub, method = "REML")
  # }, error = function(e) NULL)

  # if (!is.null(models$gam_core_re)) {
  #   aic_val <- AIC(models$gam_core_re)
  #   bic_val <- BIC(models$gam_core_re)
  #   n_params <- length(coef(models$gam_core_re))
  #   model_comparison <- rbind(model_comparison,
  #     data.frame(model = "gam_core_re", aic = aic_val, bic = bic_val,
  #                n_params = n_params, has_species_interaction = FALSE))
  # }

  # # 11. GAM with unique core random effect and tensor product interaction
  # print("GAM: c = a * exp(b * depth) - species interaction (by tensor product) (with unique core random effect)")
  # models$gam_core_re_tensor <- tryCatch({
  #   gam(carbon ~ te(depth, random_core_variable, bs = "re", k = c(5, 3)), data = dat_sub, method = "REML")
  # }, error = function(e) NULL)

  # if (!is.null(models$gam_core_re_tensor)) {
  #   aic_val <- AIC(models$gam_core_re_tensor)
  #   bic_val <- BIC(models$gam_core_re_tensor)
  #   n_params <- length(coef(models$gam_core_re_tensor))
  #   model_comparison <- rbind(model_comparison,
  #     data.frame(model = "gam_core_re_tensor", aic = aic_val, bic = bic_val,
  #                n_params = n_params, has_species_interaction = TRUE))
  # }

  # # 12. GAM with unique core random effect and species interaction
  # print("GAM: c = a * exp(b * depth) - species interaction (with unique core random effect)")
  # models$gam_core_re_species <- tryCatch({
  #   gam(carbon ~ s(depth) + s(random_core_variable, bs = "re") + species_fac, data = dat_sub, method = "REML")
  # }, error = function(e) NULL)

  # if (!is.null(models$gam_core_re_species)) {
  #   aic_val <- AIC(models$gam_core_re_species)
  #   bic_val <- BIC(models$gam_core_re_species)
  #   n_params <- length(coef(models$gam_core_re_species))
  #   model_comparison <- rbind(model_comparison,
  #     data.frame(model = "gam_core_re_species", aic = aic_val, bic = bic_val,
  #                n_params = n_params, has_species_interaction = TRUE))
  # }

  # calculate delta aic and bic
  if (nrow(model_comparison) > 0) {
    model_comparison <- model_comparison %>%
      arrange(aic) %>%
      mutate(
        delta_aic = aic - min(aic, na.rm = TRUE),
        delta_bic = bic - min(bic, na.rm = TRUE)
      )
  }

  return(list(
    models = models,
    comparison = model_comparison,
    dat_sub = dat_sub,
    normalised = normalise_by_core
  ))
}

#' Run comprehensive depth profile model investigation
#'
#' Compares models with and without normalization, and with/without species interactions
#'
#' @param dat_clean data frame containing the data
#' @param aggregate_fun function to calculate per-core aggregate
#' @return a list containing all comparison results
investigate_depth_profile_models <- function(dat_clean,
                                             aggregate_fun = calculate_core_aggregate) {
  cat("\n=== COMPREHENSIVE DEPTH PROFILE MODEL INVESTIGATION ===\n\n")

  # 1. models without normalization
  cat("1. MODELS WITHOUT NORMALIZATION\n")
  cat("   -----------------------------\n")
  results_not_norm <- compare_depth_profile_models(
    dat_clean,
    normalise_by_core = FALSE,
    aggregate_fun = aggregate_fun
  )

  cat("\n   Model comparison (AIC, lower is better):\n")
  print(results_not_norm$comparison %>%
    select(model, aic, delta_aic, has_species_interaction) %>%
    arrange(delta_aic))

  # 2. models with normalization
  cat("\n2. MODELS WITH NORMALIZATION (by core aggregate)\n")
  cat("   ------------------------------------------------\n")
  results_norm <- compare_depth_profile_models(
    dat_clean,
    normalise_by_core = TRUE,
    aggregate_fun = aggregate_fun
  )

  cat("\n   Model comparison (AIC, lower is better):\n")
  print(results_norm$comparison %>%
    select(model, aic, delta_aic, has_species_interaction) %>%
    arrange(delta_aic))

  # 3. compare best models with and without normalization
  cat("\n3. COMPARING NORMALIZATION EFFECT\n")
  cat("   --------------------------------\n")

  best_not_norm <- results_not_norm$comparison$model[1]
  best_norm <- results_norm$comparison$model[1]

  cat("   Best model without normalization:", best_not_norm, "\n")
  cat("   Best model with normalization:", best_norm, "\n")

  aic_not_norm <- results_not_norm$comparison$aic[1]
  aic_norm <- results_norm$comparison$aic[1]

  cat("   AIC (not normalized):", round(aic_not_norm, 2), "\n")
  cat("   AIC (normalized):", round(aic_norm, 2), "\n")
  cat("   Delta AIC:", round(aic_norm - aic_not_norm, 2), "\n")

  if (aic_norm < aic_not_norm) {
    cat("   -> Normalization improves model fit\n")
  } else {
    cat("   -> No normalization provides better fit\n")
  }

  # 4. test species interaction significance
  cat("\n4. TESTING SPECIES INTERACTION SIGNIFICANCE\n")
  cat("   ------------------------------------------\n")

  # compare best model with vs without species interaction (for normalized data)
  if (nrow(results_norm$comparison) > 0) {
    best_with_interaction <- results_norm$comparison %>%
      filter(has_species_interaction == TRUE) %>%
      arrange(delta_aic) %>%
      slice(1)

    best_without_interaction <- results_norm$comparison %>%
      filter(has_species_interaction == FALSE) %>%
      arrange(delta_aic) %>%
      slice(1)

    if (nrow(best_with_interaction) > 0 && nrow(best_without_interaction) > 0) {
      delta_aic_interaction <- best_with_interaction$aic - best_without_interaction$aic

      cat(
        "   Best model without species interaction:", best_without_interaction$model,
        "(AIC =", round(best_without_interaction$aic, 2), ")\n"
      )
      cat(
        "   Best model with species interaction:", best_with_interaction$model,
        "(AIC =", round(best_with_interaction$aic, 2), ")\n"
      )
      cat("   Delta AIC:", round(delta_aic_interaction, 2), "\n")

      if (delta_aic_interaction < -10) {
        cat("   -> Strong evidence for species-specific depth profiles\n")
      } else if (delta_aic_interaction < -5) {
        cat("   -> Moderate evidence for species-specific depth profiles\n")
      } else if (delta_aic_interaction < 0) {
        cat("   -> Weak evidence for species-specific depth profiles\n")
      } else {
        cat("   -> No evidence for species-specific depth profiles\n")
      }
    }
  }

  # create comparison plots
  p1 <- plot_model_comparison(results_not_norm)
  p2 <- plot_model_comparison(results_norm)

  return(list(
    results_not_normalised = results_not_norm,
    results_normalised = results_norm,
    plot_not_normalised = p1,
    plot_normalised = p2
  ))
}


# ================================ HYPERPARAMETER TUNING ================================

#' Tune random forest hyperparameters using cross-validation
#'
#' @param train_data training data frame containing predictor variables and response
#' @param predictor_vars character vector of predictor variable names
#' @param n_folds number of folds for cross-validation during tuning (default: 3)
#' @param verbose whether to print detailed tuning progress (default: false)
#' @return list of best hyperparameters (mtry, ntree, nodesize)
tune_rf <- function(train_data, predictor_vars, n_folds = 3, verbose = FALSE) {
  formula_str <- paste("median_carbon_density ~", paste(predictor_vars, collapse = " + "))
  formula_obj <- as.formula(formula_str)

  # parameter grid
  mtry_vals <- c(2, 4, 6, floor(sqrt(length(predictor_vars))), length(predictor_vars))
  ntree_vals <- c(300, 500, 700)
  nodesize_vals <- c(5, 10, 20)

  best_rmse <- Inf
  best_params <- NULL

  pb <- progress_bar(
    total = length(mtry_vals) * length(ntree_vals) * length(nodesize_vals)
  )
  for (mtry in mtry_vals) {
    for (ntree in ntree_vals) {
      for (nodesize in nodesize_vals) {
        # simple cv
        cv_rmse <- c()
        folds <- sample(rep(1:n_folds, length.out = nrow(train_data)))

        for (fold in 1:n_folds) {
          train_fold <- train_data[folds != fold, ]
          val_fold <- train_data[folds == fold, ]

          model <- randomForest(formula_obj,
            data = train_fold,
            mtry = mtry, ntree = ntree, nodesize = nodesize,
            importance = FALSE
          )
          pred <- predict(model, val_fold)
          cv_rmse <- c(cv_rmse, sqrt(mean((val_fold$median_carbon_density - pred)^2)))
        }

        mean_rmse <- mean(cv_rmse)
        if (verbose) cat(sprintf("  mtry=%d, ntree=%d, nodesize=%d: RMSE=%.4f\n", mtry, ntree, nodesize, mean_rmse))

        if (mean_rmse < best_rmse) {
          best_rmse <- mean_rmse
          best_params <- list(mtry = mtry, ntree = ntree, nodesize = nodesize)
        }
        pb$tick()
      }
    }
  }

  return(best_params)
}

#' Tune xgboost hyperparameters using cross-validation
#'
#' @param train_data training data frame containing predictor variables and response
#' @param predictor_vars character vector of predictor variable names
#' @param n_folds number of folds for cross-validation during tuning (default: 3)
#' @param verbose whether to print detailed tuning progress (default: false)
#' @return list of best hyperparameters (max_depth, learning_rate, subsample, colsample_bytree, nrounds)
tune_xgboost <- function(train_data, predictor_vars, n_folds = 3, verbose = FALSE) {
  X_train <- as.matrix(train_data[, predictor_vars])
  y_train <- train_data$median_carbon_density

  # parameter grid (reduced for speed)
  max_depth_vals <- c(3, 6, 9)
  learning_rate_vals <- c(0.1, 0.3)
  subsample_vals <- c(0.8, 1.0)
  colsample_bytree_vals <- c(0.8, 1.0)
  nrounds_vals <- c(100, 200)

  best_rmse <- Inf
  best_params <- NULL

  folds <- sample(rep(1:n_folds, length.out = nrow(train_data)))

  total_combos <- length(max_depth_vals) * length(learning_rate_vals) * length(subsample_vals) *
    length(colsample_bytree_vals) * length(nrounds_vals)
  pb <- progress_bar(
    total = total_combos
  )
  for (max_depth in max_depth_vals) {
    for (learning_rate in learning_rate_vals) {
      for (subsample in subsample_vals) {
        for (colsample_bytree in colsample_bytree_vals) {
          for (nrounds in nrounds_vals) {
            cv_rmse <- c()

            for (fold in 1:n_folds) {
              X_fold_train <- X_train[folds != fold, ]
              y_fold_train <- y_train[folds != fold]
              X_fold_val <- X_train[folds == fold, ]
              y_fold_val <- y_train[folds == fold]

              model <- xgboost(
                x = X_fold_train, y = y_fold_train,
                nrounds = nrounds, max_depth = max_depth,
                learning_rate = learning_rate, subsample = subsample,
                colsample_bytree = colsample_bytree,
                nthread = 1L,
                objective = "reg:squarederror"
              )
              pred <- predict(model, X_fold_val)
              cv_rmse <- c(cv_rmse, sqrt(mean((y_fold_val - pred)^2)))
            }

            mean_rmse <- mean(cv_rmse)
            if (mean_rmse < best_rmse) {
              best_rmse <- mean_rmse
              best_params <- list(
                max_depth = max_depth, learning_rate = learning_rate,
                subsample = subsample, colsample_bytree = colsample_bytree,
                nrounds = nrounds
              )
            }
            pb$tick()
          }
        }
      }
    }
  }

  return(best_params)
}

#' Tune neural network hyperparameters using cross-validation
#'
#' @param train_data training data frame containing predictor variables and response
#' @param predictor_vars character vector of predictor variable names
#' @param n_folds number of folds for cross-validation during tuning (default: 3)
#' @param verbose whether to print detailed tuning progress (default: false)
#' @return list of best hyperparameters (hidden, learningrate)
tune_nn <- function(train_data, predictor_vars, n_folds = 3, verbose = FALSE) {
  train_data <- as.data.frame(train_data)
  train_numeric <- train_data %>%
    dplyr::select(dplyr::all_of(predictor_vars), median_carbon_density) %>%
    dplyr::mutate_all(as.numeric) %>%
    as.data.frame()
  train_numeric <- train_numeric[complete.cases(train_numeric), ]

  # scale data
  scale_params <- list(
    means = colMeans(train_numeric[, predictor_vars, drop = FALSE]),
    sds = apply(train_numeric[, predictor_vars, drop = FALSE], 2, sd)
  )
  train_scaled <- train_numeric
  for (var in predictor_vars) {
    if (scale_params$sds[var] > 0) {
      train_scaled[[var]] <- (train_numeric[[var]] - scale_params$means[var]) / scale_params$sds[var]
    } else {
      train_scaled[[var]] <- 0
    }
  }
  response_mean <- mean(train_numeric$median_carbon_density)
  response_sd <- sd(train_numeric$median_carbon_density)
  train_scaled$median_carbon_density <- (train_numeric$median_carbon_density - response_mean) / response_sd

  formula_str <- paste("median_carbon_density ~", paste(predictor_vars, collapse = " + "))
  formula_obj <- as.formula(formula_str)

  # parameter grid
  hidden_layers <- list(c(5), c(10), c(5, 5), c(10, 5))
  learningrate_vals <- c(0.01, 0.1)

  best_rmse <- Inf
  best_params <- NULL
  folds <- sample(rep(1:n_folds, length.out = nrow(train_scaled)))
  pb <- progress_bar(
    total = length(hidden_layers) * length(learningrate_vals)
  )
  for (hidden in hidden_layers) {
    for (learningrate in learningrate_vals) {
      cv_rmse <- c()

      for (fold in 1:n_folds) {
        train_fold <- train_scaled[folds != fold, ]
        val_fold <- train_scaled[folds == fold, ]

        tryCatch(
          {
            model <- neuralnet(formula_obj,
              data = train_fold,
              hidden = hidden, learningrate = learningrate,
              linear.output = TRUE, threshold = 0.01, stepmax = 1e+05
            )
            pred_scaled <- predict(model, val_fold[, predictor_vars, drop = FALSE])
            if (is.matrix(pred_scaled)) pred_scaled <- pred_scaled[, 1]
            pred <- pred_scaled * response_sd + response_mean
            obs <- val_fold$median_carbon_density * response_sd + response_mean
            cv_rmse <- c(cv_rmse, sqrt(mean((obs - pred)^2)))
          },
          error = function(e) {
            cv_rmse <<- c(cv_rmse, Inf) # skip if model fails
          }
        )
      }

      mean_rmse <- mean(cv_rmse)
      if (verbose) {
        cat(sprintf(
          "  hidden=%s, lr=%.2f: RMSE=%.4f\n",
          paste(hidden, collapse = ","), learningrate, mean_rmse
        ))
      }

      if (mean_rmse < best_rmse && !is.infinite(mean_rmse)) {
        best_rmse <- mean_rmse
        best_params <- list(hidden = hidden, learningrate = learningrate)
      }
      pb$tick()
    }
  }

  return(best_params)
}

#' Unified function to tune hyperparameters for different model types
#'
#' @param model_type character string specifying model type: "rf", "xgb", or "nn"
#' @param train_data training data frame containing predictor variables and response
#' @param predictor_vars character vector of predictor variable names
#' @param n_folds number of folds for cross-validation during tuning (default: 3)
#' @param verbose whether to print detailed tuning progress (default: false)
#' @return list of best hyperparameters for the specified model type
tune_hyperparameters <- function(model_type, train_data, predictor_vars, n_folds = 3, verbose = FALSE) {
  cat("Tuning", model_type, "hyperparameters...\n")

  if (model_type == "rf" || model_type == "random_forest") {
    return(tune_rf(train_data, predictor_vars, n_folds, verbose))
  } else if (model_type == "xgb" || model_type == "xgboost") {
    return(tune_xgboost(train_data, predictor_vars, n_folds, verbose))
  } else if (model_type == "nn" || model_type == "neural_network") {
    return(tune_nn(train_data, predictor_vars, n_folds, verbose))
  } else {
    stop("Unknown model type. Use 'rf', 'xgb', or 'nn'")
  }
}

# ================================ CROSS-VALIDATION ================================
#' Run cross-validation for multiple models with optional hyperparameter tuning
#'
#' @param cv_method_name character string describing the cv method (e.g., "random_split", "spatial_split")
#' @param fold_indices vector of fold assignments for each observation (or list of test indices per fold)
#' @param core_data data frame containing predictor variables and response (median_carbon_density)
#' @param predictor_vars character vector of predictor variable names
#' @param tune_hyperparams logical, whether to tune hyperparameters (default: false)
#' @param nested_tuning logical, if true uses nested cv (tune per fold), if false tunes once on first fold (default: true)
#' @param verbose logical, whether to print detailed progress during tuning (default: false)
#' @param buffer_m optional; if set with core_sf, training set = points not in test and > buffer_m from all test points (distance-buffered CV)
#' @param core_sf optional sf object of core locations (required when buffer_m is set; used for kriging if provided)
#' @param return_predictions if TRUE, return list(metrics = ..., predictions = ...) with observed/predicted/longitude/latitude per test point
#' @param skip_inla if TRUE, do not run INLA (avoids slow/hanging fits when INLA is slow or problematic)
#' @param models optional character vector of model names to run; NULL = all models. Names: "Random Forest", "XGBoost", "SVM", "Neural Network", "GPR", "Ordinary Kriging", "Universal Kriging (env drift)", "GAM", "INLA".
#' @return data frame with performance metrics (r2, rmse, mae, bias) for each model and fold; or list(metrics, predictions) if return_predictions
run_cv <- function(cv_method_name, fold_indices, core_data, predictor_vars, tune_hyperparams = FALSE, nested_tuning = TRUE, verbose = FALSE, buffer_m = NA_real_, core_sf = NULL, return_predictions = FALSE, models = NULL) {
  all_models <- c("RF", "XGB", "SVM", "NN", "GPR", "OK", "UK", "GAM", "INLA")
  if (is.null(models)) models <- all_models
  models <- intersect(models, all_models)
  if (length(models) == 0) stop("run_cv: no valid models selected")
  cat("Selected models for CV:", paste(models, collapse = ", "), "\n")

  cat("\n=== Running CV:", cv_method_name, if (length(models) < length(all_models)) paste0("(", paste(models, collapse = ", "), ")") else "", "===\n")

  # normalize fold_indices to a simple vector format
  # blockCV returns folds_ids as a list (one element per fold), so convert if needed
  if (is.list(fold_indices)) {
    # convert list format to vector: create a vector where each row index gets its fold number
    n_rows <- nrow(core_data)
    fold_vector <- integer(n_rows)
    for (fold_num in seq_along(fold_indices)) {
      fold_vector[fold_indices[[fold_num]]] <- fold_num
    }
    fold_indices <- fold_vector
  } else if (!is.vector(fold_indices) || !is.numeric(fold_indices)) {
    # ensure it's a numeric vector
    fold_indices <- as.numeric(fold_indices)
  }

  # ensure fold_indices has the right length
  if (length(fold_indices) != nrow(core_data)) {
    stop("fold_indices length (", length(fold_indices), ") does not match core_data nrow (", nrow(core_data), ")")
  }

  use_buffer <- !is.na(buffer_m) && buffer_m > 0 && !is.null(core_sf) && inherits(core_sf, "sf")
  if (use_buffer) {
    cat("  Distance-buffered CV: training points must be >", buffer_m / 1000, "km from test points\n")
  }

  if (tune_hyperparams && nested_tuning) {
    cat("  Using NESTED cross-validation for hyperparameter tuning\n")
    cat("  (Hyperparameters tuned separately for each fold)\n\n")
  } else if (tune_hyperparams && !nested_tuning) {
    cat("  Using SINGLE-FOLD hyperparameter tuning\n")
    cat("  (Hyperparameters tuned once on first fold's training data, used for all folds)\n\n")
  } else {
    cat("  Using default hyperparameters\n\n")
  }

  results_list <- list()
  preds_list <- list()

  # non-nested tuning: tune hyperparameters once on first fold's training data
  # # these will be used for all folds if nested_tuning = false
  rf_params <- xgb_params <- nn_params <- NULL
  # if (tune_hyperparams && !nested_tuning) {
  #   cat("  Tuning hyperparameters on first fold's training data...\n")
  #   first_train_indices <- which(fold_indices != 1)
  #   first_train_data <- core_data[first_train_indices, ] %>% as.data.frame()
  #   if ("RF" %in% models) rf_params <- tune_hyperparameters("rf", first_train_data, predictor_vars, verbose = verbose)
  #   if ("XGB" %in% models) xgb_params <- tune_hyperparameters("xgb", first_train_data, predictor_vars, verbose = verbose)
  #   if ("NN" %in% models) nn_params <- tune_hyperparameters("nn", first_train_data, predictor_vars, verbose = verbose)
  #   cat("  Using these tuned hyperparameters for all folds\n\n")
  # }

  num_folds <- max(fold_indices)
  for (fold in 1:num_folds) {
    cat("  Fold", fold, "/", num_folds, "...\n")

    # get train/test split for this fold
    test_indices <- which(fold_indices == fold)
    train_indices <- which(fold_indices != fold)

    # distance-buffered: exclude from training any point within buffer_m of a test point
    if (use_buffer && length(test_indices) > 0 && length(train_indices) > 0) {
      dist_train_to_test <- sf::st_distance(core_sf[train_indices, ], core_sf[test_indices, ])
      min_dist_to_test <- apply(dist_train_to_test, 1, min)
      train_indices <- train_indices[as.numeric(min_dist_to_test) > buffer_m]
    }

    # skip if no test data
    if (length(test_indices) == 0) {
      cat("    Warning: No test data in fold", fold, "- skipping\n")
      next
    }

    if (length(train_indices) == 0) {
      cat("    Warning: No training data after buffer in fold", fold, "- skipping\n")
      next
    }

    train_data <- core_data[train_indices, ] %>% as.data.frame()
    test_data <- core_data[test_indices, ] %>% as.data.frame()

    # nested tuning: tune hyperparameters separately for each fold
    # this ensures hyperparameters are optimized on each fold's training data only
    if (tune_hyperparams && nested_tuning) {
      cat("    Tuning hyperparameters for this fold...\n")
      if ("RF" %in% models) rf_params <- tune_hyperparameters("rf", train_data, predictor_vars, verbose = verbose)
      if ("XGB" %in% models) xgb_params <- tune_hyperparameters("xgb", train_data, predictor_vars, verbose = verbose)
      if ("NN" %in% models) nn_params <- tune_hyperparameters("nn", train_data, predictor_vars, verbose = verbose)
    }
    # if non-nested tuning, rf_params, xgb_params, nn_params are already set above
    # if no tuning, they remain null (defaults will be used)

    # fit models with tuned hyperparameters (or defaults if tuning disabled)
    # each model in its own tryCatch so one failure does not skip others
    if ("RF" %in% models) {
      tryCatch(
        {
          rf_result <- fit_rf(train_data, test_data, predictor_vars, hyperparams = rf_params)
          rf_metrics <- calculate_metrics(test_data$median_carbon_density, rf_result$predictions)
          results_list[[paste0("fold_", fold, "_rf")]] <- data.frame(
            method = cv_method_name,
            fold = fold,
            model = "RF",
            rf_metrics
          )
          if (return_predictions) preds_list[[length(preds_list) + 1]] <- data.frame(method = cv_method_name, fold = fold, model = "RF", observed = test_data$median_carbon_density, predicted = rf_result$predictions, longitude = test_data$longitude, latitude = test_data$latitude, stringsAsFactors = FALSE)
        },
        error = function(e) {
          cat("    RF fold", fold, ":", e$message, "\n")
        }
      )
    }
    if ("XGB" %in% models) {
      tryCatch(
        {
          xgb_result <- fit_xgboost(train_data, test_data, predictor_vars, hyperparams = xgb_params)
          xgb_metrics <- calculate_metrics(test_data$median_carbon_density, xgb_result$predictions)
          results_list[[paste0("fold_", fold, "_xgb")]] <- data.frame(
            method = cv_method_name,
            fold = fold,
            model = "XGB",
            xgb_metrics
          )
          if (return_predictions) preds_list[[length(preds_list) + 1]] <- data.frame(method = cv_method_name, fold = fold, model = "XGB", observed = test_data$median_carbon_density, predicted = xgb_result$predictions, longitude = test_data$longitude, latitude = test_data$latitude, stringsAsFactors = FALSE)
        },
        error = function(e) {
          cat("    XGB fold", fold, ":", e$message, "\n")
        }
      )
    }
    if ("SVM" %in% models) {
      tryCatch(
        {
          svm_result <- fit_svm(train_data, test_data, predictor_vars)
          svm_metrics <- calculate_metrics(test_data$median_carbon_density, svm_result$predictions)
          results_list[[paste0("fold_", fold, "_svm")]] <- data.frame(
            method = cv_method_name,
            fold = fold,
            model = "SVM",
            svm_metrics
          )
          if (return_predictions) preds_list[[length(preds_list) + 1]] <- data.frame(method = cv_method_name, fold = fold, model = "SVM", observed = test_data$median_carbon_density, predicted = svm_result$predictions, longitude = test_data$longitude, latitude = test_data$latitude, stringsAsFactors = FALSE)
        },
        error = function(e) {
          cat("    SVM fold", fold, ":", e$message, "\n")
        }
      )
    }
    if ("NN" %in% models) {
      tryCatch(
        {
          nn_result <- fit_nn(train_data, test_data, predictor_vars, hyperparams = nn_params)
          nn_metrics <- calculate_metrics(test_data$median_carbon_density, nn_result$predictions)
          results_list[[paste0("fold_", fold, "_nn")]] <- data.frame(
            method = cv_method_name,
            fold = fold,
            model = "NN",
            nn_metrics
          )
          if (return_predictions) preds_list[[length(preds_list) + 1]] <- data.frame(method = cv_method_name, fold = fold, model = "NN", observed = test_data$median_carbon_density, predicted = nn_result$predictions, longitude = test_data$longitude, latitude = test_data$latitude, stringsAsFactors = FALSE)
        },
        error = function(e) {
          cat("    NN fold", fold, ":", e$message, "\n")
        }
      )
    }
    if ("GPR" %in% models) {
      tryCatch(
        {
          gpr_result <- fit_gpr(train_data, test_data, predictor_vars)
          gpr_metrics <- calculate_metrics(test_data$median_carbon_density, gpr_result$predictions)
          results_list[[paste0("fold_", fold, "_gpr")]] <- data.frame(
            method = cv_method_name,
            fold = fold,
            model = "GPR",
            gpr_metrics
          )
          if (return_predictions) preds_list[[length(preds_list) + 1]] <- data.frame(method = cv_method_name, fold = fold, model = "GPR", observed = test_data$median_carbon_density, predicted = gpr_result$predictions, longitude = test_data$longitude, latitude = test_data$latitude, stringsAsFactors = FALSE)
        },
        error = function(e) {
          cat("    GPR fold", fold, ":", e$message, "\n")
        }
      )
    }
    if ("OK" %in% models) {
      tryCatch(
        {
          ok_result <- fit_ok(train_data, test_data)
          ok_metrics <- if (!all(is.na(ok_result$predictions))) {
            calculate_metrics(test_data$median_carbon_density, ok_result$predictions)
          } else {
            data.frame(r2 = NA_real_, rmse = NA_real_, mae = NA_real_, bias = NA_real_)
          }
          results_list[[paste0("fold_", fold, "_ok")]] <- data.frame(
            method = cv_method_name,
            fold = fold,
            model = "OK",
            ok_metrics
          )
          if (return_predictions && !all(is.na(ok_result$predictions))) preds_list[[length(preds_list) + 1]] <- data.frame(method = cv_method_name, fold = fold, model = "OK", observed = test_data$median_carbon_density, predicted = ok_result$predictions, longitude = test_data$longitude, latitude = test_data$latitude, stringsAsFactors = FALSE)
        },
        error = function(e) {
          cat("    OK fold", fold, ":", e$message, "\n")
          results_list[[paste0("fold_", fold, "_ok")]] <<- data.frame(
            method = cv_method_name, fold = fold, model = "OK",
            r2 = NA_real_, rmse = NA_real_, mae = NA_real_, bias = NA_real_
          )
        }
      )
    }
    if ("UK" %in% models) {
      tryCatch(
        {
          uk_result <- fit_uk(train_data, test_data, predictor_vars)
          uk_metrics <- if (!all(is.na(uk_result$predictions))) {
            calculate_metrics(test_data$median_carbon_density, uk_result$predictions)
          } else {
            data.frame(r2 = NA_real_, rmse = NA_real_, mae = NA_real_, bias = NA_real_)
          }
          results_list[[paste0("fold_", fold, "_uk")]] <- data.frame(
            method = cv_method_name,
            fold = fold,
            model = "UK",
            uk_metrics
          )
          if (return_predictions && !all(is.na(uk_result$predictions))) preds_list[[length(preds_list) + 1]] <- data.frame(method = cv_method_name, fold = fold, model = "UK", observed = test_data$median_carbon_density, predicted = uk_result$predictions, longitude = test_data$longitude, latitude = test_data$latitude, stringsAsFactors = FALSE)
        },
        error = function(e) {
          cat("    UK fold", fold, ":", e$message, "\n")
          results_list[[paste0("fold_", fold, "_uk")]] <<- data.frame(
            method = cv_method_name, fold = fold, model = "Universal Kriging (env drift)",
            r2 = NA_real_, rmse = NA_real_, mae = NA_real_, bias = NA_real_
          )
        }
      )
    }
    if ("GAM" %in% models) {
      tryCatch(
        {
          gam_result <- fit_gam(train_data, test_data, predictor_vars)
          if (!all(is.na(gam_result$predictions))) {
            gam_metrics <- calculate_metrics(test_data$median_carbon_density, gam_result$predictions)
            results_list[[paste0("fold_", fold, "_gam")]] <- data.frame(
              method = cv_method_name,
              fold = fold,
              model = "GAM",
              gam_metrics
            )
            if (return_predictions) preds_list[[length(preds_list) + 1]] <- data.frame(method = cv_method_name, fold = fold, model = "GAM", observed = test_data$median_carbon_density, predicted = gam_result$predictions, longitude = test_data$longitude, latitude = test_data$latitude, stringsAsFactors = FALSE)
          }
        },
        error = function(e) {
          cat("    GAM fold", fold, ":", e$message, "\n")
        }
      )
    }
    if ("INLA" %in% models) {
      tryCatch(
        {
          inla_result <- fit_inla_cv(train_data, test_data, predictor_vars)
          if (!all(is.na(inla_result$predictions))) {
            inla_metrics <- calculate_metrics(test_data$median_carbon_density, inla_result$predictions)
            results_list[[paste0("fold_", fold, "_inla")]] <- data.frame(
              method = cv_method_name,
              fold = fold,
              model = "INLA",
              inla_metrics
            )
            if (return_predictions) preds_list[[length(preds_list) + 1]] <- data.frame(method = cv_method_name, fold = fold, model = "INLA", observed = test_data$median_carbon_density, predicted = inla_result$predictions, longitude = test_data$longitude, latitude = test_data$latitude, stringsAsFactors = FALSE)
          } else {
            results_list[[paste0("fold_", fold, "_inla")]] <- data.frame(
              method = cv_method_name, fold = fold, model = "INLA",
              r2 = NA_real_, rmse = NA_real_, mae = NA_real_, bias = NA_real_
            )
          }
        },
        error = function(e) {
          cat("    INLA fold", fold, ":", e$message, "\n")
          results_list[[paste0("fold_", fold, "_inla")]] <<- data.frame(
            method = cv_method_name, fold = fold, model = "INLA",
            r2 = NA_real_, rmse = NA_real_, mae = NA_real_, bias = NA_real_
          )
        }
      )
    }
    if ("INLA" %in% models) {
      results_list[[paste0("fold_", fold, "_inla")]] <- data.frame(
        method = cv_method_name, fold = fold, model = "INLA",
        r2 = NA_real_, rmse = NA_real_, mae = NA_real_, bias = NA_real_
      )
    }
  }

  # combine results
  results_df <- dplyr::bind_rows(results_list)
  if (return_predictions && length(preds_list) > 0) {
    return(list(metrics = results_df, predictions = dplyr::bind_rows(preds_list)))
  }
  return(results_df)
}


# ================================ SPATIAL AUTOCORRELATION ANALYSIS ================================

#' Calculate spatial autocorrelation range for a raster
#'
#' Uses blockCV's cv_spatial_autocor function to estimate the range of spatial
#' autocorrelation using variogram analysis.
#'
#' @param raster_path Path to the NetCDF raster file
#' @param varname Variable name within the NetCDF file
#' @param name Descriptive name for the variable
#' @param sample_points Optional sf object with sample points
#' @param plot_variogram Whether to plot the variogram
#' @return List with autocorrelation range and variogram plot
calculate_spatial_autocorrelation <- function(raster_path,
                                              varname = NULL,
                                              name = "variable",
                                              num_sample = 5000,
                                              sample_points = NULL,
                                              plot_variogram = TRUE) {
  cat("Processing:", name, "\n")
  cat("  Loading raster:", raster_path, "\n")

  # Check if file exists

  if (!file.exists(raster_path)) {
    cat("  ERROR: File not found!\n")
    return(list(name = name, range_m = NA, error = "File not found"))
  }

  # Load raster
  tryCatch(
    {
      if (!is.null(varname)) {
        r <- terra::rast(raster_path, subds = varname)
      } else {
        r <- terra::rast(raster_path)
      }

      # If multiple layers, take first one
      if (terra::nlyr(r) > 1) {
        cat("  Multiple layers detected, using first layer\n")
        r <- r[[1]]
      }
      # Display raster information
      cat("  Raster information:\n")
      print(r)
      # Ensure raster has CRS
      if (is.na(terra::crs(r))) {
        cat("  Setting CRS to EPSG:4326\n")
        terra::crs(r) <- "EPSG:4326"
      }

      # Use cv_spatial_autocor from blockCV to estimate autocorrelation range
      # This function fits an exponential variogram and estimates the range
      cat("  Calculating spatial autocorrelation...\n")

      autocor_result <- blockCV::cv_spatial_autocor(
        r = r,
        num_sample = min(num_sample, terra::ncell(r)), # Sample points for variogram
        plot = FALSE,
        progress = FALSE
      )

      range_m <- autocor_result$range
      cat("  Estimated autocorrelation range:", round(range_m / 1000, 2), "km\n")

      # Create variogram plot
      if (plot_variogram) {
        variogram_plot <- autocor_result$plot +
          labs(
            title = name,
            subtitle = paste("Range:", round(range_m / 1000, 2), "km")
          ) +
          theme_minimal() +
          theme(
            plot.title = element_text(size = 10, face = "bold"),
            plot.subtitle = element_text(size = 9)
          )
      } else {
        variogram_plot <- NULL
      }

      return(list(
        name = name,
        range_m = range_m,
        range_km = range_m / 1000,
        variogram_plot = variogram_plot,
        autocor_object = autocor_result
      ))
    },
    error = function(e) {
      cat("  ERROR:", e$message, "\n")
      return(list(name = name, range_m = NA, error = e$message))
    }
  )
}


# ================================ PLOTTING ================================
#' Custom function for scatter plots with continuous color using turbo colormap
#' turbo colormap is from viridisLite package (dependency of ggplot2)
#'
#' @param data data frame containing the data
#' @param mapping ggplot mapping
#' @param ... additional arguments to pass to ggplot
#' @return a ggplot object
scatter_with_color <- function(data, mapping, ...) {
  # add color aesthetic to the mapping by combining with existing mapping
  # use aes() to create new mapping that includes color
  new_mapping <- aes(
    x = !!mapping$x,
    y = !!mapping$y,
    color = carbon_density_g_c_cm3
  )

  ggplot(data = data, mapping = new_mapping) +
    geom_point(alpha = 0.5, ...) +
    scale_color_gradientn(
      colors = viridisLite::turbo(256),
      name = "Carbon density\n(gC cm^-3)",
      guide = guide_colorbar(
        barwidth = 1,
        barheight = 10,
        title.position = "right",
        title.hjust = 0.5
      )
    ) +
    theme(legend.position = "right")
}

#' Define turbo colormap for categorical or integer variables and return a named vector of colors
#'
#' If the supplied variable is a factor or character, one color is assigned per unique value.
#' If it is an integer or can be converted to integer, assigns a color for each integer from min to max,
#' and names the vector with those integer values (normalized so that the largest integer is mapped to the last color).
#'
#' @param data Data frame containing the variable to color by
#' @param color_var Name of the variable (column) to use for coloring (as string)
#' @return Named vector of hex color codes with names corresponding to unique values (factor/char) or integer labels
define_turbo_colors <- function(data, color_var) {
  values <- data[[color_var]]

  # Check for integer type or integer-like variable
  if (is.integer(values) || (is.numeric(values) && all(!is.na(values)) && all(values == as.integer(values)))) {
    values_int <- as.integer(values)
    # Normalise to the largest integer
    int_range <- sort(unique(values_int))
    if (length(int_range) == 1) {
      turbo_colors <- viridisLite::turbo(1)
      names(turbo_colors) <- int_range
      return(turbo_colors)
    }
    max_int <- max(int_range)
    # assign one color for each integer from smallest to max
    color_vals <- min(int_range):max_int
    n_colors <- length(color_vals)
    turbo_colors <- viridisLite::turbo(n_colors)
    names(turbo_colors) <- color_vals
    return(turbo_colors)
  } else {
    # treat as categorical (factor or character)
    categories <- unique(values)
    n_colors <- length(categories)
    turbo_colors <- viridisLite::turbo(n_colors)
    names(turbo_colors) <- categories
    return(turbo_colors)
  }
}

#' Get consistent color mapping for seagrass species
#'
#' Returns a named vector of colors for seagrass species that can be used
#' consistently across all plots. Uses a curated palette for visual distinction.
#'
#' @return Named vector of hex color codes for each species
#' @export
get_species_colors <- function() {
  c(
    "Cymodocea nodosa" = "#3366CC", # blue
    "Halophila stipulacea" = "#FF9900", # orange
    "Posidonia oceanica" = "#00CED1", # dark cyan
    "Zostera marina" = "#32CD32", # lime green
    "Zostera noltei" = "#CCCC00", # yellow-green
    "Zostera marina and Zostera noltei" = "#006400", # dark green
    "Zostera marina and Cymodocea nodosa" = "#8B0000", # dark red
    "Unspecified" = "#808080" # gray
  )
}

#' Get consistent color mapping for regions
#'
#' Returns a named vector of colors for regions that can be used
#' consistently across all plots. Uses a curated palette for visual distinction.
#'
#' @return Named vector of hex color codes for each region
#' @export
get_region_colors <- function() {
  c(
    "Mediterranean Sea"       = "#1E90FF", # dodger blue
    "North European Atlantic" = "#228B22", # forest green
    "South European Atlantic" = "#60cc32", # lime green
    "Baltic Sea"              = "#FFD700", # gold
    "Black Sea"               = "#2F4F4F" # dark slate gray
  )
}

#' Get colors for a specific grouping variable
#'
#' Convenience function that returns appropriate colors based on the group type.
#' Falls back to turbo palette if group type is unknown.
#'
#' @param group_type Either "species" or "region"
#' @param data Optional data frame to extract unique values from (for fallback)
#' @param group_var Optional column name (for fallback)
#' @return Named vector of hex color codes
#' @export
get_group_colors <- function(group_type, data = NULL, group_var = NULL) {
  if (group_type == "species") {
    return(get_species_colors())
  } else if (group_type == "region") {
    return(get_region_colors())
  } else if (!is.null(data) && !is.null(group_var)) {
    return(define_turbo_colors(data, group_var))
  } else {
    stop("Unknown group_type. Use 'species' or 'region', or provide data and group_var for fallback.")
  }
}

#' Plot a pairs plot of specified variables with a continuous variable as the color
#' and save it with a combined colorbar legend.
#'
#' @param data The input data frame
#' @param env_vars Character vector of variable names to plot on axes (must be columns in data)
#' @param color_var The variable name to use as the color aesthetic (must be in data)
#' @param output_file Filename to save the figure (png)
#' @param width Width of output image in inches (default 15)
#' @param height Height of output image in inches (default 15)
#' @param dpi Dots per inch for png output (default 300)
#' @return NULL (saves output to file; displays plot)
plot_env_pairs <- function(data, env_vars, color_var, output_file = "environmental_variables_pairs.png",
                           width = 15, height = 15, dpi = 300) {
  require(GGally)
  require(ggplot2)
  require(dplyr)
  require(cowplot)

  # Prepare data so color_var is present
  dat_pairs <- data %>%
    dplyr::select(dplyr::all_of(c(env_vars, color_var)))

  # Use a custom wrapper for scatter plots that knows to color by color_var
  pairs_plot <- ggpairs(
    dat_pairs,
    columns = env_vars,
    lower = list(continuous = wrap(scatter_with_color)), # assumes scatter_with_color uses color_var
    upper = list(continuous = wrap("cor", size = 3)),
    diag = list(continuous = "densityDiag")
  ) +
    theme_minimal()

  # Extract the colorbar legend from one of the lower triangle plots
  first_lower_plot <- getPlot(pairs_plot, 2, 1)
  colorbar_legend <- cowplot::get_legend(first_lower_plot)

  # Create combined plot with pairs plot and colorbar
  if (!is.null(colorbar_legend)) {
    library(grid)
    library(gridExtra)

    png(output_file, width = width + 1, height = height, units = "in", res = dpi)
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(1, 2, widths = unit(c(10, 1), "null"))))

    pushViewport(viewport(layout.pos.col = 1))
    print(pairs_plot, newpage = FALSE)
    popViewport()

    pushViewport(viewport(layout.pos.col = 2))
    grid.draw(colorbar_legend)
    popViewport()

    dev.off()

    # also print to screen
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(1, 2, widths = unit(c(10, 1), "null"))))
    pushViewport(viewport(layout.pos.col = 1))
    print(pairs_plot, newpage = FALSE)
    popViewport()
    pushViewport(viewport(layout.pos.col = 2))
    grid.draw(colorbar_legend)
    popViewport()
  } else {
    # Just plot and save if legend extraction fails
    print(pairs_plot)
    ggsave(output_file, plot = pairs_plot, width = width, height = height, dpi = dpi)
  }
  invisible(NULL)
}

#' Plot depth relationships by group (rotated axes: depth on y, C-density on x)
#'
#' @param results list containing the results of the model fitting
#' @param group_name the name of the group variable
#' @param xlim optional limits for the x-axis
#' @param use_facets whether to use faceted plots (default TRUE, auto-enabled if >6 groups)
#' @param colors optional named vector of colors (use get_species_colors() or get_region_colors())
#' @return a ggplot object
plot_depth_by_group <- function(results, group_name, xlim = NULL, use_facets = TRUE, colors = NULL) {
  if (is.null(results)) {
    return(NULL)
  }

  dat_sub <- results$dat_sub
  group_var <- results$group_var
  group_coefs <- results$group_coefs
  group_models <- results$group_models

  # depth sequence for predictions
  depth_seq <- seq(min(dat_sub$sediment_mean_depth_cm, na.rm = TRUE),
    max(dat_sub$sediment_mean_depth_cm, na.rm = TRUE),
    length.out = 100
  )

  # auto-enable facets for many groups
  n_groups <- length(unique(dat_sub[[group_var]]))
  use_facets <- use_facets || n_groups > 6

  # build predictions for each group using delta method for CIs
  build_group_predictions <- function(i) {
    g <- group_coefs$group[i]
    g_char <- as.character(g)

    if (!is.null(group_models) && g_char %in% names(group_models)) {
      model_g <- group_models[[g_char]]

      pred_result <- tryCatch(
        {
          pred_fit <- predict(model_g, newdata = data.frame(sediment_mean_depth_cm = depth_seq))
          vcov_mat <- vcov(model_g)
          coefs <- coef(model_g)
          a <- coefs[["a"]]
          b <- coefs[["b"]]

          # delta method: gradient = [exp(b*d), a*d*exp(b*d)]
          se_vec <- sapply(depth_seq, function(d) {
            grad <- c(exp(b * d), a * d * exp(b * d))
            sqrt(t(grad) %*% vcov_mat %*% grad)
          })

          list(fit = pred_fit, lower = pred_fit - 1.96 * se_vec, upper = pred_fit + 1.96 * se_vec)
        },
        error = function(e) NULL
      )

      if (!is.null(pred_result)) {
        return(data.frame(
          sediment_mean_depth_cm = depth_seq,
          group = g,
          carbon_density_fit = pred_result$fit,
          carbon_density_lower = pred_result$lower,
          carbon_density_upper = pred_result$upper
        ))
      }
    }

    # fallback: simple prediction without CIs
    data.frame(
      sediment_mean_depth_cm = depth_seq,
      group = g,
      carbon_density_fit = group_coefs$a[i] * exp(group_coefs$b[i] * depth_seq),
      carbon_density_lower = NA,
      carbon_density_upper = NA
    )
  }

  # build predictions and filter for valid CIs
  preds_df <- bind_rows(lapply(seq_len(nrow(group_coefs)), build_group_predictions))
  preds_df[[group_var]] <- preds_df$group

  # set up colors: use provided colors or fall back to turbo palette
  if (is.null(colors)) {
    group_colors <- define_turbo_colors(dat_sub, group_var)
  } else {
    group_colors <- colors
  }

  # build plot: points -> ribbons -> lines
  p <- ggplot() +
    geom_point(
      data = dat_sub,
      aes(x = carbon_density_g_c_cm3, y = sediment_mean_depth_cm, color = !!sym(group_var)),
      alpha = 0.5
    )

  # add confidence ribbons if available
  ribbon_df <- preds_df %>% filter(!is.na(carbon_density_lower))
  if (nrow(ribbon_df) > 0) {
    p <- p +
      geom_ribbon(
        data = ribbon_df,
        aes(
          y = sediment_mean_depth_cm,
          xmin = carbon_density_lower,
          xmax = carbon_density_upper,
          fill = !!sym(group_var)
        ),
        alpha = 0.18,
        show.legend = FALSE
      )
  }

  # add fitted lines
  p <- p +
    geom_line(
      data = preds_df,
      aes(x = carbon_density_fit, y = sediment_mean_depth_cm, color = !!sym(group_var)),
      linewidth = 1.5
    ) +
    scale_y_reverse() +
    scale_color_manual(name = group_name, values = group_colors) +
    scale_fill_manual(values = group_colors) +
    labs(
      y = "Depth (cm)",
      x = expression("Carbon density (gC cm"^
        {
          -3
        } * ")"),
      title = paste("Depth-Carbon density relationship by", group_name),
      subtitle = paste("Delta AIC:", round(results$delta_aic, 2))
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")

  if (!is.null(xlim)) {
    p <- p + scale_x_continuous(limits = xlim)
  }

  if (use_facets) {
    p <- p + facet_wrap(as.formula(paste("~", group_var)), scales = "free_y")
  }

  print(p)
  return(p)
}


#' Plot model comparison results
#'
#' @param comparison_results list from compare_depth_profile_models
#' @return a ggplot object
plot_model_comparison <- function(comparison_results) {
  comp <- comparison_results$comparison

  if (nrow(comp) == 0) {
    cat("No models to compare\n")
    return(NULL)
  }

  # create a long format for plotting
  comp_long <- comp %>%
    select(model, delta_aic, delta_bic, has_species_interaction) %>%
    pivot_longer(
      cols = c(delta_aic, delta_bic),
      names_to = "criterion",
      values_to = "delta"
    ) %>%
    mutate(criterion = ifelse(criterion == "delta_aic", "Delta AIC", "Delta BIC"))

  p <- ggplot(comp_long, aes(
    x = reorder(model, delta), y = delta,
    fill = has_species_interaction
  )) +
    geom_bar(stat = "identity") +
    facet_wrap(~criterion, scales = "free_y") +
    coord_flip() +
    scale_fill_manual(
      values = c("FALSE" = "lightblue", "TRUE" = "lightcoral"),
      name = "Species\ninteraction"
    ) +
    labs(
      x = "Model", y = "Delta (relative to best model)",
      title = "Model comparison: depth profile functional forms"
    ) +
    theme_minimal()

  return(p)
}


#' Plot raster stack as grid of subplots with individual colorbars
#'
#' @param raster_stack SpatRaster object with multiple layers
#' @param sample_points Optional data.frame with longitude, latitude columns
#' @param point_color Optional: either a color string (e.g., "red") or column name from sample_points
#' @param raster_info Optional list with descriptions for layer names (for better labels)
#' @param ncol Number of columns in the grid
#' @param color_scale Color scale for rasters (default: "viridis")
#' @return Combined plot object (patchwork)
plot_raster_stack <- function(raster_stack,
                              sample_points = NULL,
                              point_color = NULL,
                              raster_info = NULL,
                              ncol = 3,
                              color_scale = "viridis") {
  # Convert raster stack to data frame
  raster_df <- terra::as.data.frame(raster_stack, xy = TRUE)

  # Get layer names
  layer_names <- names(raster_stack)
  n_layers <- length(layer_names)

  # Determine point color range if using a column
  point_color_range <- NULL
  if (!is.null(sample_points) && !is.null(point_color) && point_color %in% names(sample_points)) {
    point_color_range <- range(sample_points[[point_color]], na.rm = TRUE)
  }

  # Create individual plots for each layer
  plot_list <- list()

  for (i in seq_along(layer_names)) {
    layer_name <- layer_names[i]

    # Get description if available
    if (!is.null(raster_info) && layer_name %in% names(raster_info)) {
      title <- paste0(raster_info[[layer_name]]$description, "\n(", layer_name, ")")
    } else {
      title <- layer_name
    }

    # Prepare data for this layer
    layer_data <- raster_df %>%
      dplyr::select(x, y, !!sym(layer_name)) %>%
      dplyr::rename(value = !!sym(layer_name))

    # Create base plot for this layer
    p <- ggplot(layer_data, aes(x = x, y = y, fill = value)) +
      geom_raster() +
      scale_fill_viridis_c(
        option = color_scale,
        na.value = "transparent",
        name = NULL,
        guide = guide_colorbar(
          title.position = "right",
          title.hjust = 0.5,
          barwidth = 0.5,
          barheight = 3
        )
      ) +
      labs(
        title = title,
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 6),
        legend.position = "right",
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
      ) +
      coord_fixed()
    # label entire x axis Latitude, entire y axis Longitude
    p <- p + labs(x = "Latitude", y = "Longitude") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))


    # Add sample points if provided
    if (!is.null(sample_points) && !is.null(point_color)) {
      # Check if point_color is a column name or a color string
      if (point_color %in% names(sample_points)) {
        # It's a column - use continuous color scale with shared range
        points_df <- sample_points %>%
          dplyr::select(longitude, latitude, !!sym(point_color)) %>%
          dplyr::rename(point_value = !!sym(point_color))

        # Only add colorbar to first plot if using column
        if (i == 1) {
          p <- p +
            geom_point(
              data = points_df,
              aes(x = longitude, y = latitude, color = point_value),
              inherit.aes = FALSE,
              size = 0.8,
              alpha = 0.7,
              stroke = 0
            ) +
            scale_color_viridis_c(
              option = "plasma",
              name = point_color,
              limits = point_color_range,
              guide = guide_colorbar(
                title.position = "right",
                title.hjust = 0.5,
                barwidth = 0.5,
                barheight = 8,
                order = 2
              )
            )
        } else {
          # Other plots: just points, no colorbar
          p <- p +
            geom_point(
              data = points_df,
              aes(x = longitude, y = latitude, color = point_value),
              inherit.aes = FALSE,
              size = 0.8,
              alpha = 0.7,
              stroke = 0,
              show.legend = FALSE
            ) +
            scale_color_viridis_c(
              option = "plasma",
              limits = point_color_range,
              guide = "none"
            )
        }
      } else {
        # It's a color string - same for all plots
        points_df <- sample_points %>%
          dplyr::select(longitude, latitude)

        p <- p +
          geom_point(
            data = points_df,
            aes(x = longitude, y = latitude),
            inherit.aes = FALSE,
            color = point_color,
            size = 0.8,
            alpha = 0.7
          )
      }
    }

    # add on geographical features (country outlines)
    p <- p +
      geom_polygon(
        data = world, aes(x = long, y = lat, group = group),
        fill = "#eeeeee", color = "#a5a5a5", linewidth = 0.1
      ) +
      coord_fixed()
    # set extent to the extent of the plot
    p <- p + coord_cartesian(
      xlim = c(min(layer_data$x), max(layer_data$x)),
      ylim = c(min(layer_data$y), max(layer_data$y))
    )

    plot_list[[i]] <- p
  }

  # Combine plots using patchwork
  combined_plot <- wrap_plots(plot_list, ncol = ncol)

  return(combined_plot)
}


## ================================ DEPRECATED FUNCTIONS ================================
# #' Function to fit and compare models with/without interaction
# #'
# #' @param dat_clean data frame containing the data
# #' @param group_var the variable to group by
# #' @param group_name the name of the group variable
# #' @param normalise_by_core logical, whether to normalize carbon density by core aggregate
# #' @param aggregate_fun function to calculate per-core aggregate (default: median of top 5 cm)
# #' @return a list containing the results of the model fitting
# # function to fit and compare models with/without interaction
# compare_depth_relationships <- function(dat_clean, group_var, group_name,
#                                         normalise_by_core = FALSE,
#                                         aggregate_fun = calculate_core_aggregate) {
#   cat("\n=== testing variation by", group_name, "===\n")

#   # remove missing values for this grouping variable
#   dat_sub <- dat_clean %>%
#     filter(!is.na(!!sym(group_var)) & !is.na(carbon_density_g_c_cm3) &
#            !is.na(sediment_mean_depth_cm))

#   if (normalise_by_core) {
#     # calculate aggregate value per core using the provided function
#     core_aggregates <- aggregate_fun(dat_sub)

#     # normalize by dividing by the aggregate value
#     dat_sub <- dat_sub %>%
#       left_join(core_aggregates, by = "random_core_variable") %>%
#       mutate(carbon_density_g_c_cm3 = carbon_density_g_c_cm3 / aggregate_carbon_density) %>%
#       select(-aggregate_carbon_density)  # remove the aggregate column
#   }

#   # check if enough data
#   n_groups <- length(unique(dat_sub[[group_var]]))
#   if (n_groups < 2) {
#     cat("insufficient groups (", n_groups, "), skipping\n")
#     return(NULL)
#   }

#   # null model: exponential decay with no interaction (same relationship for all groups)
#   model_null <- tryCatch({
#     nls(carbon_density_g_c_cm3 ~ a * exp(b * sediment_mean_depth_cm),
#         data = dat_sub,
#         start = list(a = max(dat_sub$carbon_density_g_c_cm3, na.rm = TRUE),
#                      b = -0.01),
#         control = nls.control(maxiter = 500, warnOnly = TRUE))
#   }, error = function(e) NULL)

#   if (is.null(model_null)) {
#     cat("failed to fit null model\n")
#     return(NULL)
#   }

#   # interaction model: exponential decay with different relationship per group
#   # fit separate models per group and combine
#   groups <- unique(dat_sub[[group_var]])
#   group_models <- list()
#   group_coefs <- data.frame()

#   for (g in groups) {
#     dat_g <- dat_sub %>% filter(!!sym(group_var) == g)
#     if (nrow(dat_g) < 5) next  # need minimum data points

#     model_g <- tryCatch({
#       nls(carbon_density_g_c_cm3 ~ a * exp(b * sediment_mean_depth_cm),
#           data = dat_g,
#           start = list(a = max(dat_g$carbon_density_g_c_cm3, na.rm = TRUE),
#                        b = -0.01),
#           control = nls.control(maxiter = 500, warnOnly = TRUE))
#     }, error = function(e) NULL)

#     if (!is.null(model_g)) {
#       group_models[[as.character(g)]] <- model_g
#       coefs <- coef(model_g)
#       coef_se <- summary(model_g)$coefficients[, "Std. Error"]
#       group_coefs <- rbind(group_coefs,
#                            data.frame(group = g, a = coefs[["a"]], b = coefs[["b"]],
#                                       a_se = coef_se[["a"]], b_se = coef_se[["b"]],
#                                       aic = AIC(model_g), n = nrow(dat_g)))
#     }
#   }

#   if (nrow(group_coefs) == 0) {
#     cat("failed to fit group-specific models\n")
#     return(NULL)
#   }

#   # compare models using aic
#   aic_null <- AIC(model_null)
#   aic_interaction <- sum(group_coefs$aic)
#   delta_aic <- aic_interaction - aic_null

#   cat("null model (no interaction) aic:", round(aic_null, 2), "\n")
#   cat("interaction model (group-specific) aic:", round(aic_interaction, 2), "\n")
#   cat("delta aic:", round(delta_aic, 2), "\n")
#   if (delta_aic < -10) {
#     cat("strong evidence for group-specific relationships\n")
#   } else if (delta_aic < -5) {
#     cat("moderate evidence for group-specific relationships\n")
#   } else if (delta_aic < 0) {
#     cat("weak evidence for group-specific relationships\n")
#   } else {
#     cat("no evidence for group-specific relationships (null model preferred)\n")
#   }

#   # also test using gam with smooth interaction
#   gam_null <- tryCatch({
#     gam(carbon_density_g_c_cm3 ~ s(sediment_mean_depth_cm, k = 5),
#         data = dat_sub, method = "REML")
#   }, error = function(e) NULL)

#   # create formula for gam interaction model
#   group_fac <- as.factor(dat_sub[[group_var]])
#   gam_interaction <- tryCatch({
#     gam(carbon_density_g_c_cm3 ~ s(sediment_mean_depth_cm, by = group_fac, k = 5) + group_fac,
#         data = dat_sub, method = "REML")
#   }, error = function(e) NULL)

#   if (!is.null(gam_null) && !is.null(gam_interaction)) {
#     gam_compare <- anova(gam_null, gam_interaction, test = "F")
#     cat("\ngam comparison (f-test):\n")
#     print(gam_compare)
#     if (gam_compare$`Pr(>F)`[2] < 0.05) {
#       cat("gam: significant interaction (p < 0.05)\n")
#     } else {
#       cat("gam: no significant interaction (p >= 0.05)\n")
#     }
#   }

#   # return results for plotting
#   return(list(
#     group_var = group_var,
#     group_name = group_name,
#     model_null = model_null,
#     group_models = group_models,
#     group_coefs = group_coefs,
#     dat_sub = dat_sub,
#     delta_aic = delta_aic,
#     gam_null = gam_null,
#     gam_interaction = gam_interaction
#   ))
# }
