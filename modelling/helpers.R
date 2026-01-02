# ================================ UTILITY FUNCTIONS ================================

#' Install packages if not already installed
#' 
#' @param packages character vector of package names to install
#' @return NULL but prints message to console when packages are installed
install_packages <- function(packages) {
  # if string, convert to list
  if(is.character(packages)) {
    packages <- list(packages)
  }
  for(package in packages) {
    if(!requireNamespace(package, quietly = TRUE)) {
      install.packages(package)
    }
  }
  print(paste("Installed packages:", paste(packages, collapse = ", ")))
}

#' Load packages, install if necessary
#' 
#' @param packages character vector of package names to load
#' @return NULL but prints message to console when packages are loaded
load_packages <- function(packages) {
  for (package in packages) {
    if (!suppressWarnings(require(package, character.only = TRUE))) {
      install_packages(package)
      library(package, character.only = TRUE)
    }
  }
  print(paste("Loaded packages:", paste(packages, collapse = ", ")))
}

# ================================ INITIALISATION ================================

# load packages if not already loaded
load_packages(c("here", "mgcv", "tidyverse", "ggplot2", "visreg", "randomForest", "corrplot", "blockCV", "GauPro", "xgboost", "e1071", "neuralnet", "purrr", "sp"))

#' Create a progress bar for hyperparameter tuning
#'
#' @param total total number of iterations for the progress bar
#' @param format format string for the progress bar display
#' @param clear whether to clear the progress bar when done
#' @param width width of the progress bar in characters
#' @return a progress bar object from the progress package
progress_bar <- function(total, format = "[:bar] :current/:total (:percent)", clear = FALSE, width = 60) {
  if(!requireNamespace("progress", quietly = TRUE)) {
    install.packages("progress")
  }
  progress::progress_bar$new(total = total, format = format, clear = clear, width = width)
}

# ================================ MODEL FITTING ================================

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
  ntree <- if(!is.null(hyperparams$ntree)) hyperparams$ntree else 500
  mtry <- if(!is.null(hyperparams$mtry)) hyperparams$mtry else floor(sqrt(length(predictor_vars)))
  nodesize <- if(!is.null(hyperparams$nodesize)) hyperparams$nodesize else NULL
  
  model <- randomForest(formula_obj, 
                       data = train_data,
                       ntree = ntree,
                       mtry = mtry,
                       nodesize = nodesize,
                       importance = TRUE)
  
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
  nrounds <- if(!is.null(hyperparams$nrounds)) hyperparams$nrounds else 100
  max_depth <- if(!is.null(hyperparams$max_depth)) hyperparams$max_depth else 6
  learning_rate <- if(!is.null(hyperparams$learning_rate)) hyperparams$learning_rate else 0.3
  subsample <- if(!is.null(hyperparams$subsample)) hyperparams$subsample else 0.8
  colsample_bytree <- if(!is.null(hyperparams$colsample_bytree)) hyperparams$colsample_bytree else 0.8
  
  # fit model
  model <- xgboost(x = X_train,
                   y = y_train,
                   nrounds = nrounds,
                   max_depth = max_depth,
                   learning_rate = learning_rate,
                   subsample = subsample,
                   colsample_bytree = colsample_bytree,
                   objective = "reg:squarederror")
  
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
  
  for(kernel_type in kernels_to_try) {
    # tune parameters for this kernel type
    # note: linear kernel doesn't use gamma, so we skip it for non-linear kernels
    # polynomial and sigmoid kernels use both cost and gamma
    tune_result <- tryCatch({
      tune.svm(formula_obj,
               data = train_data,
               gamma = 10^(-3:1),
               cost = 10^(-1:2),
               kernel = kernel_type,
               tunecontrol = tune.control(sampling = "cross", cross = 3))
    }, error = function(e) {
      cat("    Warning: Failed to tune", kernel_type, "kernel:", e$message, "\n")
      return(NULL)
    })
    
    if(!is.null(tune_result)) {
      # extract performance (error) - lower is better
      current_performance <- tune_result$best.performance
      if(current_performance < best_performance) {
        best_performance <- current_performance
        best_result <- tune_result
      }
    }
  }
  
  # if no kernel worked, fall back to default svm
  if(is.null(best_result)) {
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
  for(var in predictor_vars) {
    if(scale_params$sds[var] > 0) {  # avoid division by zero
      train_scaled[[var]] <- (train_numeric[[var]] - scale_params$means[var]) / scale_params$sds[var]
    } else {
      train_scaled[[var]] <- 0  # if no variance, set to 0
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
  for(var in predictor_vars) {
    if(scale_params$sds[var] > 0) {
      test_scaled[[var]] <- (test_numeric[[var]] - scale_params$means[var]) / scale_params$sds[var]
    } else {
      test_scaled[[var]] <- 0
    }
  }
  
  # create formula
  formula_str <- paste("median_carbon_density ~", paste(predictor_vars, collapse = " + "))
  formula_obj <- as.formula(formula_str)
  
  # use provided hyperparameters or defaults
  hidden <- if(!is.null(hyperparams$hidden)) hyperparams$hidden else c(10, 5)
  learningrate <- if(!is.null(hyperparams$learningrate)) hyperparams$learningrate else NULL
  
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
  if(!is.null(learningrate)) model_args$learningrate <- learningrate
  
  model <- do.call(neuralnet, model_args)
  
  # predict on scaled test data
  # neuralnet predict returns a matrix, extract first column
  predictions_scaled <- predict(model, newdata = test_scaled[, predictor_vars, drop = FALSE])
  
  # convert predictions_scaled from matrix to vector if needed
  if(is.matrix(predictions_scaled)) {
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
  
  model <- gpkm(formula_obj, data = train_data, kernel)
  predictions <- predict(model, newdata = test_data)
  return(list(model = model, predictions = predictions))
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
  for(mtry in mtry_vals) {
    for(ntree in ntree_vals) {
      for(nodesize in nodesize_vals) {
        # simple cv
        cv_rmse <- c()
        folds <- sample(rep(1:n_folds, length.out = nrow(train_data)))
        
        for(fold in 1:n_folds) {
          train_fold <- train_data[folds != fold, ]
          val_fold <- train_data[folds == fold, ]
          
          model <- randomForest(formula_obj, data = train_fold,
                               mtry = mtry, ntree = ntree, nodesize = nodesize,
                               importance = FALSE)
          pred <- predict(model, val_fold)
          cv_rmse <- c(cv_rmse, sqrt(mean((val_fold$median_carbon_density - pred)^2)))
        }
        
        mean_rmse <- mean(cv_rmse)
        if(verbose) cat(sprintf("  mtry=%d, ntree=%d, nodesize=%d: RMSE=%.4f\n", mtry, ntree, nodesize, mean_rmse))
        
        if(mean_rmse < best_rmse) {
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
  for(max_depth in max_depth_vals) {
    for(learning_rate in learning_rate_vals) {
      for(subsample in subsample_vals) {
        for(colsample_bytree in colsample_bytree_vals) {
          for(nrounds in nrounds_vals) {
            cv_rmse <- c()
            
            for(fold in 1:n_folds) {
              X_fold_train <- X_train[folds != fold, ]
              y_fold_train <- y_train[folds != fold]
              X_fold_val <- X_train[folds == fold, ]
              y_fold_val <- y_train[folds == fold]
              
              model <- xgboost(x = X_fold_train, y = y_fold_train,
                              nrounds = nrounds, max_depth = max_depth,
                              learning_rate = learning_rate, subsample = subsample,
                              colsample_bytree = colsample_bytree,
                              objective = "reg:squarederror")
              pred <- predict(model, X_fold_val)
              cv_rmse <- c(cv_rmse, sqrt(mean((y_fold_val - pred)^2)))
            }
            
            mean_rmse <- mean(cv_rmse)
            if(mean_rmse < best_rmse) {
              best_rmse <- mean_rmse
              best_params <- list(max_depth = max_depth, learning_rate = learning_rate,
                                 subsample = subsample, colsample_bytree = colsample_bytree,
                                 nrounds = nrounds)
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
  for(var in predictor_vars) {
    if(scale_params$sds[var] > 0) {
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
  for(hidden in hidden_layers) {
    for(learningrate in learningrate_vals) {
      cv_rmse <- c()
      
      for(fold in 1:n_folds) {
        train_fold <- train_scaled[folds != fold, ]
        val_fold <- train_scaled[folds == fold, ]
        
        tryCatch({
          model <- neuralnet(formula_obj, data = train_fold,
                            hidden = hidden, learningrate = learningrate,
                            linear.output = TRUE, threshold = 0.01, stepmax = 1e+05)
          pred_scaled <- predict(model, val_fold[, predictor_vars, drop = FALSE])
          if(is.matrix(pred_scaled)) pred_scaled <- pred_scaled[, 1]
          pred <- pred_scaled * response_sd + response_mean
          obs <- val_fold$median_carbon_density * response_sd + response_mean
          cv_rmse <- c(cv_rmse, sqrt(mean((obs - pred)^2)))
        }, error = function(e) {
          cv_rmse <<- c(cv_rmse, Inf)  # skip if model fails
        })
      }
      
      mean_rmse <- mean(cv_rmse)
      if(verbose) cat(sprintf("  hidden=%s, lr=%.2f: RMSE=%.4f\n", 
                             paste(hidden, collapse=","), learningrate, mean_rmse))
      
      if(mean_rmse < best_rmse && !is.infinite(mean_rmse)) {
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
  
  if(model_type == "rf" || model_type == "random_forest") {
    return(tune_rf(train_data, predictor_vars, n_folds, verbose))
  } else if(model_type == "xgb" || model_type == "xgboost") {
    return(tune_xgboost(train_data, predictor_vars, n_folds, verbose))
  } else if(model_type == "nn" || model_type == "neural_network") {
    return(tune_nn(train_data, predictor_vars, n_folds, verbose))
  } else {
    stop("Unknown model type. Use 'rf', 'xgb', or 'nn'")
  }
}

# ================================ CROSS-VALIDATION ================================
#' Run cross-validation for multiple models with optional hyperparameter tuning
#'
#' @param cv_method_name character string describing the cv method (e.g., "random_split", "spatial_split")
#' @param fold_indices vector of fold assignments for each observation
#' @param core_data data frame containing predictor variables and response (median_carbon_density)
#' @param predictor_vars character vector of predictor variable names
#' @param tune_hyperparams logical, whether to tune hyperparameters (default: false)
#' @param nested_tuning logical, if true uses nested cv (tune per fold), if false tunes once on first fold (default: true)
#' @param verbose logical, whether to print detailed progress during tuning (default: false)
#' @return data frame with performance metrics (r2, rmse, mae, bias) for each model and fold
run_cv <- function(cv_method_name, fold_indices, core_data, predictor_vars, tune_hyperparams = FALSE, nested_tuning = TRUE, verbose = FALSE) {
  cat("\n=== Running CV:", cv_method_name, "===\n")
  
  # normalize fold_indices to a simple vector format
  # blockCV returns folds_ids as a list (one element per fold), so convert if needed
  if(is.list(fold_indices)) {
    # convert list format to vector: create a vector where each row index gets its fold number
    n_rows <- nrow(core_data)
    fold_vector <- integer(n_rows)
    for(fold_num in seq_along(fold_indices)) {
      fold_vector[fold_indices[[fold_num]]] <- fold_num
    }
    fold_indices <- fold_vector
  } else if(!is.vector(fold_indices) || !is.numeric(fold_indices)) {
    # ensure it's a numeric vector
    fold_indices <- as.numeric(fold_indices)
  }
  
  # ensure fold_indices has the right length
  if(length(fold_indices) != nrow(core_data)) {
    stop("fold_indices length (", length(fold_indices), ") does not match core_data nrow (", nrow(core_data), ")")
  }
  
  if(tune_hyperparams && nested_tuning) {
    cat("  Using NESTED cross-validation for hyperparameter tuning\n")
    cat("  (Hyperparameters tuned separately for each fold)\n\n")
  } else if(tune_hyperparams && !nested_tuning) {
    cat("  Using SINGLE-FOLD hyperparameter tuning\n")
    cat("  (Hyperparameters tuned once on first fold's training data, used for all folds)\n\n")
  } else {
    cat("  Using default hyperparameters\n\n")
  }
  
  results_list <- list()
  
  # non-nested tuning: tune hyperparameters once on first fold's training data
  # these will be used for all folds if nested_tuning = false
  rf_params <- xgb_params <- nn_params <- NULL
  if(tune_hyperparams && !nested_tuning) {
    cat("  Tuning hyperparameters on first fold's training data...\n")
    first_train_indices <- which(fold_indices != 1)
    first_train_data <- core_data[first_train_indices, ] %>% as.data.frame()
    
    rf_params <- tune_hyperparameters("rf", first_train_data, predictor_vars, verbose = verbose)
    xgb_params <- tune_hyperparameters("xgb", first_train_data, predictor_vars, verbose = verbose)
    nn_params <- tune_hyperparameters("nn", first_train_data, predictor_vars, verbose = verbose)
    cat("  Using these tuned hyperparameters for all folds\n\n")
  }
  
  num_folds <- max(fold_indices)
  for(fold in 1:num_folds) {
    cat("  Fold", fold, "/", num_folds, "...\n")
    
    # get train/test split for this fold
    train_indices <- which(fold_indices != fold)
    test_indices <- which(fold_indices == fold)

    # skip if no test data
    if(length(test_indices) == 0) {
      cat("    Warning: No test data in fold", fold, "- skipping\n")
      next
    }
    
    train_data <- core_data[train_indices, ] %>% as.data.frame()
    test_data <- core_data[test_indices, ] %>% as.data.frame()
    
    # nested tuning: tune hyperparameters separately for each fold
    # this ensures hyperparameters are optimized on each fold's training data only
    if(tune_hyperparams && nested_tuning) {
      cat("    Tuning hyperparameters for this fold...\n")
      rf_params <- tune_hyperparameters("rf", train_data, predictor_vars, verbose = verbose)
      xgb_params <- tune_hyperparameters("xgb", train_data, predictor_vars, verbose = verbose)
      nn_params <- tune_hyperparameters("nn", train_data, predictor_vars, verbose = verbose)
    }
    # if non-nested tuning, rf_params, xgb_params, nn_params are already set above
    # if no tuning, they remain null (defaults will be used)
    
    # fit models with tuned hyperparameters (or defaults if tuning disabled)
    tryCatch({
      # random forest
      # rf_result <- fit_rf(train_data, test_data, predictor_vars, hyperparams = rf_params)
      # rf_metrics <- calculate_metrics(test_data$median_carbon_density, rf_result$predictions)
      # results_list[[paste0("fold_", fold, "_rf")]] <- data.frame(
      #   method = cv_method_name,
      #   fold = fold,
      #   model = "Random Forest",
      #   rf_metrics
      # )
      
      # # xgboost
      # xgb_result <- fit_xgboost(train_data, test_data, predictor_vars, hyperparams = xgb_params)
      # xgb_metrics <- calculate_metrics(test_data$median_carbon_density, xgb_result$predictions)
      # results_list[[paste0("fold_", fold, "_xgb")]] <- data.frame(
      #   method = cv_method_name,
      #   fold = fold,
      #   model = "XGBoost",
      #   xgb_metrics
      # )
      
      # # svm
      # svm_result <- fit_svm(train_data, test_data, predictor_vars)
      # svm_metrics <- calculate_metrics(test_data$median_carbon_density, svm_result$predictions)
      # results_list[[paste0("fold_", fold, "_svm")]] <- data.frame(
      #   method = cv_method_name,
      #   fold = fold,
      #   model = "SVM",
      #   svm_metrics
      # )

      # # neural network
      # nn_result <- fit_nn(train_data, test_data, predictor_vars, hyperparams = nn_params)
      # nn_metrics <- calculate_metrics(test_data$median_carbon_density, nn_result$predictions)
      # results_list[[paste0("fold_", fold, "_nn")]] <- data.frame(
      #   method = cv_method_name,
      #   fold = fold,
      #   model = "Neural Network",
      #   nn_metrics
      # )

      # gpr
      gpr_result <- fit_gpr(train_data, test_data, predictor_vars)
      gpr_metrics <- calculate_metrics(test_data$median_carbon_density, gpr_result$predictions)
      results_list[[paste0("fold_", fold, "_gpr")]] <- data.frame(
        method = cv_method_name,
        fold = fold,
        model = "GPR",
        gpr_metrics
      )
    }, error = function(e) {
      cat("    Error in fold", fold, ":", e$message, "\n")
    })
  }
  
  # combine results
  results_df <- dplyr::bind_rows(results_list)
  return(results_df)
}


# ================================ PLOTTING ================================

#' Plot depth relationships by group (rotated axes: depth on y, C-density on x)
#' 
#' @param results list containing the results of the model fitting
#' @param group_name the name of the group variable
#' @param xlim optional limits for the x-axis
#' @return a ggplot object
plot_depth_by_group <- function(results, group_name, xlim = NULL) {
  if (is.null(results)) return(NULL)
  
  dat_sub <- results$dat_sub
  group_var <- results$group_var
  group_coefs <- results$group_coefs
  
  # create prediction data for each group
  depth_seq <- seq(min(dat_sub$sediment_mean_depth_cm, na.rm = TRUE),
                  max(dat_sub$sediment_mean_depth_cm, na.rm = TRUE),
                  length.out = 100)
  
  pred_list <- list()
  for (i in seq_len(nrow(group_coefs))) {
    g <- group_coefs$group[i]
    pred_df <- data.frame(
      sediment_mean_depth_cm = depth_seq,
      group = g,
      carbon_density_fit = group_coefs$a[i] * exp(group_coefs$b[i] * depth_seq)
    )
    pred_list[[i]] <- pred_df
  }
  pred_all <- bind_rows(pred_list)
  
  # plot (rotated like above: y=depth, x=carbon density, y reversed)
  p <- ggplot() +
    geom_point(data = dat_sub,
               aes(x = carbon_density_g_c_cm3, y = sediment_mean_depth_cm, color = !!sym(group_var)),
               alpha = 0.3, size = 0.5) +
    geom_line(data = pred_all,
              aes(x = carbon_density_fit, y = sediment_mean_depth_cm, color = group),
              linewidth = 1) +
    scale_y_reverse() +
    labs(
      y = "Depth (cm)",
      x = expression("Carbon density (gC cm"^{-3}*")"),
      title = paste("Depth-Carbon density relationship by", group_name),
      subtitle = paste("Delta AIC:", round(results$delta_aic, 2))
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # add xlim if provided
  if (!is.null(xlim)) {
    p <- p + scale_x_continuous(limits = xlim)
  }
  
  # if many groups, use facets (facet vertically: strips on side so rotate facet)
  n_groups <- length(unique(dat_sub[[group_var]]))
  if (n_groups > 6) {
    p <- p + facet_wrap(as.formula(paste("~", group_var)), scales = "free_y")
  }
  
  print(p)
  return(p)
}

#' Function to fit and compare models with/without interaction
#' 
#' @param dat_clean data frame containing the data
#' @param group_var the variable to group by
#' @param group_name the name of the group variable
#' @return a list containing the results of the model fitting
# function to fit and compare models with/without interaction
compare_depth_relationships <- function(dat_clean, group_var, group_name) {
  cat("\n=== testing variation by", group_name, "===\n")
  
  # remove missing values for this grouping variable
  dat_sub <- dat_clean %>%
    filter(!is.na(!!sym(group_var)) & !is.na(carbon_density_g_c_cm3) & 
           !is.na(sediment_mean_depth_cm))
  
  # check if enough data
  n_groups <- length(unique(dat_sub[[group_var]]))
  if (n_groups < 2) {
    cat("insufficient groups (", n_groups, "), skipping\n")
    return(NULL)
  }
  
  # null model: exponential decay with no interaction (same relationship for all groups)
  model_null <- tryCatch({
    nls(carbon_density_g_c_cm3 ~ a * exp(b * sediment_mean_depth_cm),
        data = dat_sub,
        start = list(a = max(dat_sub$carbon_density_g_c_cm3, na.rm = TRUE), 
                     b = -0.01),
        control = nls.control(maxiter = 500, warnOnly = TRUE))
  }, error = function(e) NULL)
  
  if (is.null(model_null)) {
    cat("failed to fit null model\n")
    return(NULL)
  }
  
  # interaction model: exponential decay with different relationship per group
  # fit separate models per group and combine
  groups <- unique(dat_sub[[group_var]])
  group_models <- list()
  group_coefs <- data.frame()
  
  for (g in groups) {
    dat_g <- dat_sub %>% filter(!!sym(group_var) == g)
    if (nrow(dat_g) < 5) next  # need minimum data points
    
    model_g <- tryCatch({
      nls(carbon_density_g_c_cm3 ~ a * exp(b * sediment_mean_depth_cm),
          data = dat_g,
          start = list(a = max(dat_g$carbon_density_g_c_cm3, na.rm = TRUE), 
                       b = -0.01),
          control = nls.control(maxiter = 500, warnOnly = TRUE))
    }, error = function(e) NULL)
    
    if (!is.null(model_g)) {
      group_models[[as.character(g)]] <- model_g
      coefs <- coef(model_g)
      group_coefs <- rbind(group_coefs, 
                           data.frame(group = g, a = coefs[["a"]], b = coefs[["b"]],
                                      aic = AIC(model_g), n = nrow(dat_g)))
    }
  }
  
  if (nrow(group_coefs) == 0) {
    cat("failed to fit group-specific models\n")
    return(NULL)
  }
  
  # compare models using aic
  aic_null <- AIC(model_null)
  aic_interaction <- sum(group_coefs$aic)
  delta_aic <- aic_interaction - aic_null
  
  cat("null model (no interaction) aic:", round(aic_null, 2), "\n")
  cat("interaction model (group-specific) aic:", round(aic_interaction, 2), "\n")
  cat("delta aic:", round(delta_aic, 2), "\n")
  if (delta_aic < -10) {
    cat("strong evidence for group-specific relationships\n")
  } else if (delta_aic < -5) {
    cat("moderate evidence for group-specific relationships\n")
  } else if (delta_aic < 0) {
    cat("weak evidence for group-specific relationships\n")
  } else {
    cat("no evidence for group-specific relationships (null model preferred)\n")
  }
  
  # also test using gam with smooth interaction
  gam_null <- tryCatch({
    gam(carbon_density_g_c_cm3 ~ s(sediment_mean_depth_cm, k = 5),
        data = dat_sub, method = "REML")
  }, error = function(e) NULL)
  
  # create formula for gam interaction model
  group_fac <- as.factor(dat_sub[[group_var]])
  gam_interaction <- tryCatch({
    gam(carbon_density_g_c_cm3 ~ s(sediment_mean_depth_cm, by = group_fac, k = 5) + group_fac,
        data = dat_sub, method = "REML")
  }, error = function(e) NULL)
  
  if (!is.null(gam_null) && !is.null(gam_interaction)) {
    gam_compare <- anova(gam_null, gam_interaction, test = "F")
    cat("\ngam comparison (f-test):\n")
    print(gam_compare)
    if (gam_compare$`Pr(>F)`[2] < 0.05) {
      cat("gam: significant interaction (p < 0.05)\n")
    } else {
      cat("gam: no significant interaction (p >= 0.05)\n")
    }
  }
  
  # return results for plotting
  return(list(
    group_var = group_var,
    group_name = group_name,
    model_null = model_null,
    group_coefs = group_coefs,
    dat_sub = dat_sub,
    delta_aic = delta_aic,
    gam_null = gam_null,
    gam_interaction = gam_interaction
  ))
}