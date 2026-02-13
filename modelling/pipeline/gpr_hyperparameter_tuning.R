# GPR Hyperparameter Tuning for 1km Spatial Cross-Validation
#
# This script tunes Gaussian Process Regression hyperparameters to maximize
# R² performance on 1km spatial block cross-validation. It addresses overfitting
# concerns by testing different regularization strategies.

## ================================ SETUP ================================
rm(list = ls())
setwd(here::here())
set.seed(42)

source("modelling/R/helpers.R")
source("modelling/R/extract_covariates_from_rasters.R")
source("modelling/R/gpr_funs.R")
load_packages(c("here", "mgcv", "tidyverse", "ggplot2", "GauPro", "blockCV", "sf", "viridisLite", "dplyr"))

## ================================ CONFIGURATION ================================
target_var <- "median_carbon_density_100cm"
block_size_km <- 1  # 1km spatial blocks for CV
n_folds <- 5
out_dir <- "figures/cv_pipeline_output"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Hyperparameter search space
# Note: GauPro uses different parameterization than standard GP
# We'll test different approaches to regularization

# Kernel options to test
kernels_to_test <- c("matern52", "matern32", "gaussian")

# Feature selection: number of top variables to test
n_features_to_test <- c(5, 10, 15, 20, 25, 30, "all")  # "all" uses all pruned variables

# Feature selection method
feature_selection_method <- "correlation"  # "correlation" (fast) or "gpr_importance" (slower but more accurate)

# For regularization, GauPro has several options:
# 1. nugget parameter (additive noise) - not directly accessible in gpkm
# 2. Use fewer predictors (feature selection) - THIS IS WHAT WE'RE TESTING
# 3. Use different kernels (matern32 is smoother than matern52)

## ================================ DATA LOADING ================================
cat("\n========================================\n")
cat("LOADING DATA\n")
cat("========================================\n\n")

raster_covariates_dat <- read_rds("data/all_extracted_new.rds")
raster_covariates_dat <- process_rs_covariates(raster_covariates_dat)

# Load pruned variables: prefer generic model permutation pruning, then GPR-specific, then GAM list
gam_pruned_vars <- get_pruned_predictors_for_model("GPR", data_cols = names(raster_covariates_dat))
if (is.null(gam_pruned_vars) || length(gam_pruned_vars) < 2L) {
  pruned_file_gpr <- "figures/covariate_selection/pruned_variables_to_include_gpr.csv"
  pruned_file <- "figures/covariate_selection/pruned_variables_to_include.csv"
  if (file.exists(pruned_file_gpr)) {
    gam_pruned_vars <- read_csv(pruned_file_gpr, show_col_types = FALSE)$variable
    cat("Loaded", length(gam_pruned_vars), "GPR-pruned covariates (gpr_covariate_pruning fallback).\n\n")
  } else if (file.exists(pruned_file)) {
    gam_pruned_vars <- read_csv(pruned_file, show_col_types = FALSE)$variable
    cat("Loaded", length(gam_pruned_vars), "pruned covariates (GAM list fallback).\n\n")
  } else {
    stop("Pruned variables not found. Run model_permutation_pruning.R or covariate_pruning.R / gpr_covariate_pruning.R first.")
  }
} else {
  cat("Loaded", length(gam_pruned_vars), "GPR covariates from model permutation pruning (1 km spatial CV).\n\n")
}
gam_pruned_vars <- gam_pruned_vars[gam_pruned_vars %in% names(raster_covariates_dat)]

# Prepare data - match cv_pipeline structure if using its cache
# cv_pipeline uses core-level aggregation, but we can use point-level data
# The folds will be matched by spatial location
predictor_vars <- gam_pruned_vars
need_vars <- c("longitude", "latitude", target_var, predictor_vars)
dat_complete <- raster_covariates_dat[, need_vars, drop = FALSE]
dat_complete <- dat_complete[complete.cases(dat_complete), ]

cat("Complete cases: ", nrow(dat_complete), "\n\n")

# Note: If using cv_pipeline cache, the folds are based on core-level aggregation
# but we can still use them by matching spatial locations
# The fold assignments are based on spatial blocks, so they should work with point-level data

## ================================ FEATURE SELECTION ================================
cat("========================================\n")
cat("FEATURE SELECTION/RANKING\n")
cat("========================================\n\n")

cat("Ranking variables by importance for feature selection...\n")
cat("Method:", feature_selection_method, "\n\n")

if (feature_selection_method == "correlation") {
  # Fast method: rank by absolute correlation with target
  cat("Computing correlations with target variable...\n")
  correlations <- sapply(predictor_vars, function(var) {
    if (var %in% names(dat_complete)) {
      cor(dat_complete[[var]], dat_complete[[target_var]], use = "complete.obs")
    } else {
      NA
    }
  })
  
  # Rank by absolute correlation
  var_importance <- data.frame(
    variable = predictor_vars,
    correlation = correlations,
    abs_correlation = abs(correlations),
    stringsAsFactors = FALSE
  )
  var_importance <- var_importance[order(-var_importance$abs_correlation), ]
  var_importance <- var_importance[!is.na(var_importance$correlation), ]
  
  cat("Top 10 variables by absolute correlation:\n")
  print(head(var_importance, 10))
  cat("\n")
  
} else if (feature_selection_method == "gpr_importance") {
  # Slower but more accurate: fit a quick GPR and use permutation importance
  cat("Fitting quick GPR model for importance estimation...\n")
  cat("(This may take a few minutes)\n\n")
  
  # Use a subset for speed
  n_sample <- min(500, nrow(dat_complete))
  sample_idx <- sample(nrow(dat_complete), n_sample)
  dat_sample <- dat_complete[sample_idx, ]
  
  # Fit baseline GPR
  formula_str <- paste(target_var, "~", paste(predictor_vars, collapse = " + "))
  formula_obj <- as.formula(formula_str)
  train_data_combined <- dat_sample[, c(target_var, predictor_vars), drop = FALSE]
  train_data_combined <- train_data_combined[complete.cases(train_data_combined), ]
  
  if (nrow(train_data_combined) < 10) {
    cat("WARNING: Insufficient data for GPR importance. Using correlation instead.\n")
    feature_selection_method <- "correlation"
    # Re-run correlation method
    correlations <- sapply(predictor_vars, function(var) {
      if (var %in% names(dat_complete)) {
        cor(dat_complete[[var]], dat_complete[[target_var]], use = "complete.obs")
      } else {
        NA
      }
    })
    var_importance <- data.frame(
      variable = predictor_vars,
      correlation = correlations,
      abs_correlation = abs(correlations),
      stringsAsFactors = FALSE
    )
    var_importance <- var_importance[order(-var_importance$abs_correlation), ]
    var_importance <- var_importance[!is.na(var_importance$correlation), ]
  } else {
    tryCatch({
      train_X <- train_data_combined[, predictor_vars, drop = FALSE]
      scale_params_imp <- gpr_compute_scale_params(train_X, predictor_vars)
      train_X_scaled <- as.data.frame(gpr_apply_scale(train_X, scale_params_imp))
      train_scaled <- data.frame(y = train_data_combined[[target_var]], train_X_scaled)
      names(train_scaled)[1] <- target_var
      baseline_model <- GauPro::gpkm(formula_obj, data = train_scaled, kernel = "matern52")
      X_scaled <- as.matrix(train_X_scaled)
      baseline_pred <- baseline_model$pred(X_scaled, se.fit = FALSE)
      baseline_rmse <- sqrt(mean((train_data_combined[[target_var]] - baseline_pred)^2, na.rm = TRUE))
      
      # Permutation importance (scale permuted data with same params)
      importance_scores <- numeric(length(predictor_vars))
      names(importance_scores) <- predictor_vars
      
      pb <- utils::txtProgressBar(min = 0, max = length(predictor_vars), style = 3)
      for (i in seq_along(predictor_vars)) {
        var <- predictor_vars[i]
        utils::setTxtProgressBar(pb, i)
        
        dat_permuted <- train_data_combined
        dat_permuted[[var]] <- sample(dat_permuted[[var]])
        X_perm_scaled <- gpr_apply_scale(as.matrix(dat_permuted[, predictor_vars, drop = FALSE]), scale_params_imp)
        pred_permuted <- baseline_model$pred(X_perm_scaled, se.fit = FALSE)
        rmse_permuted <- sqrt(mean((train_data_combined[[target_var]] - pred_permuted)^2, na.rm = TRUE))
        
        importance_scores[i] <- rmse_permuted - baseline_rmse
      }
      close(pb)
      cat("\n")
      
      var_importance <- data.frame(
        variable = predictor_vars,
        importance = importance_scores,
        stringsAsFactors = FALSE
      )
      var_importance <- var_importance[order(-var_importance$importance), ]
      
      cat("Top 10 variables by permutation importance:\n")
      print(head(var_importance, 10))
      cat("\n")
      
    }, error = function(e) {
      cat("ERROR in GPR importance calculation:", e$message, "\n")
      cat("Falling back to correlation method...\n")
      correlations <- sapply(predictor_vars, function(var) {
        if (var %in% names(dat_complete)) {
          cor(dat_complete[[var]], dat_complete[[target_var]], use = "complete.obs")
        } else {
          NA
        }
      })
      var_importance <- data.frame(
        variable = predictor_vars,
        correlation = correlations,
        abs_correlation = abs(correlations),
        stringsAsFactors = FALSE
      )
      var_importance <- var_importance[order(-var_importance$abs_correlation), ]
      var_importance <- var_importance[!is.na(var_importance$correlation), ]
    })
  }
}

# Create feature sets for testing
feature_sets <- list()
for (n_feat in n_features_to_test) {
  if (n_feat == "all") {
    feature_sets[["all"]] <- var_importance$variable
  } else {
    n_feat_num <- as.numeric(n_feat)
    if (n_feat_num <= nrow(var_importance)) {
      feature_sets[[as.character(n_feat)]] <- var_importance$variable[1:n_feat_num]
    }
  }
}

cat("Feature sets to test:\n")
for (name in names(feature_sets)) {
  cat("  ", name, " features: ", length(feature_sets[[name]]), " variables\n", sep = "")
}
cat("\n")

## ================================ CREATE SPATIAL FOLDS ================================
cat("========================================\n")
cat("LOADING/CREATING 1KM SPATIAL FOLDS\n")
cat("========================================\n\n")

block_size_m <- block_size_km * 1000

# Try to load from cv_pipeline cache first (shared cache)
# cv_pipeline.R writes to figures/cv_pipeline_output; also check modelling/cv_pipeline_output
cv_pipeline_cache_path <- NULL
for (base in c("figures/cv_pipeline_output", "modelling/cv_pipeline_output")) {
  candidate <- file.path(base, "spatial_folds_cache.rds")
  if (file.exists(candidate)) {
    cv_pipeline_cache_path <- candidate
    break
  }
}
spatial_folds_loaded <- FALSE
folds_ids <- NULL

if (!is.null(cv_pipeline_cache_path)) {
  cat("Checking cv_pipeline cache for 1km spatial folds (", cv_pipeline_cache_path, ")...\n", sep = "")
  cached <- tryCatch(readRDS(cv_pipeline_cache_path), error = function(e) NULL)
  
  if (!is.null(cached) && "spatial_strategies" %in% names(cached)) {
    # Look for 1km spatial block strategy
    target_method <- paste0("spatial_block_", block_size_m, "m")
    matching_strategy <- NULL
    
    for (strategy in cached$spatial_strategies) {
      if (strategy$method == target_method) {
        matching_strategy <- strategy
        break
      }
    }
    
    if (!is.null(matching_strategy)) {
      # cv_pipeline uses core-level aggregated data
      # We need to match our point-level data to the core-level folds
      # The folds are based on spatial blocks, so we can match by location
      
      # Check if we can directly use the folds (same number of rows)
      if (length(matching_strategy$folds) == nrow(dat_complete)) {
        # Direct match - same number of data points
        folds_ids <- matching_strategy$folds
        spatial_folds_loaded <- TRUE
        cat("  ✓ Loaded 1km spatial folds from cv_pipeline cache (direct match)\n")
        cat("  Method:", matching_strategy$method, "\n")
        cat("  Folds:", length(unique(folds_ids)), "\n\n")
      } else {
        # Different number of rows - need to match by spatial location
        cat("  Data structure differs: cv_pipeline has", length(matching_strategy$folds), 
            "rows, we have", nrow(dat_complete), "\n")
        cat("  Attempting to match by spatial location...\n")
        
        # Load the data used by cv_pipeline to get core locations
        dat_cv_pipeline <- read_rds("data/all_extracted_new.rds")
        dat_cv_pipeline <- process_rs_covariates(dat_cv_pipeline)
        
        # Check if cv_pipeline uses core aggregation
        if ("random_core_variable" %in% names(dat_cv_pipeline)) {
          # Aggregate to core level (matching cv_pipeline structure)
          core_data_cv <- dat_cv_pipeline %>%
            dplyr::group_by(random_core_variable) %>%
            dplyr::summarise(
              longitude = dplyr::first(longitude),
              latitude = dplyr::first(latitude),
              .groups = "drop"
            ) %>%
            dplyr::filter(!is.na(longitude) & !is.na(latitude)) %>%
            dplyr::arrange(random_core_variable)
          
          # Ensure we have the same number of cores as folds
          if (nrow(core_data_cv) == length(matching_strategy$folds)) {
            # Create lookup: core location -> fold assignment
            # Use rounded coordinates for matching (handles floating point differences)
            core_locations <- paste(round(core_data_cv$longitude, 8), 
                                    round(core_data_cv$latitude, 8), sep = "_")
            fold_lookup <- data.frame(
              location = core_locations,
              fold = matching_strategy$folds,
              longitude = core_data_cv$longitude,
              latitude = core_data_cv$latitude,
              stringsAsFactors = FALSE
            )
            
            # Match our point-level data to folds by location
            point_locations <- paste(round(dat_complete$longitude, 8), 
                                     round(dat_complete$latitude, 8), sep = "_")
            folds_ids <- fold_lookup$fold[match(point_locations, fold_lookup$location)]
            
            # For points not matching any core location exactly, assign to nearest fold
            unmatched <- is.na(folds_ids)
            if (sum(unmatched) > 0) {
              cat("  ", sum(unmatched), "points don't match core locations exactly.\n")
              cat("  Assigning to nearest fold by spatial distance...\n")
              
              # Find nearest core for unmatched points (vectorized for speed)
              unmatched_coords <- dat_complete[unmatched, c("longitude", "latitude"), drop = FALSE]
              
              for (i in seq_len(nrow(unmatched_coords))) {
                idx <- which(unmatched)[i]
                lon <- unmatched_coords$longitude[i]
                lat <- unmatched_coords$latitude[i]
                
                # Calculate distances to all cores
                dists <- sqrt((fold_lookup$longitude - lon)^2 + 
                             (fold_lookup$latitude - lat)^2)
                nearest_idx <- which.min(dists)
                folds_ids[idx] <- fold_lookup$fold[nearest_idx]
              }
            }
            
            spatial_folds_loaded <- TRUE
            cat("  ✓ Matched folds by spatial location\n")
            cat("  Matched", sum(!unmatched), "points directly,", sum(unmatched), "by nearest neighbor\n")
            cat("  Folds:", length(unique(folds_ids)), "\n\n")
          } else {
            cat("  Cannot match: core count mismatch\n")
            cat("  Will create new folds...\n\n")
          }
        } else {
          cat("  Cannot determine cv_pipeline data structure\n")
          cat("  Will create new folds...\n\n")
        }
      }
    } else {
      cat("  1km spatial folds not found in cv_pipeline cache\n")
      if (!is.null(cached$spatial_strategies) && length(cached$spatial_strategies) > 0) {
        cat("  Available methods:", paste(sapply(cached$spatial_strategies, function(s) s$method), collapse = ", "), "\n\n")
      }
    }
  }
} else {
  cat("No cv_pipeline spatial cache found (looked in figures/cv_pipeline_output and modelling/cv_pipeline_output).\n")
  cat("  Run cv_pipeline.R first to share 1km folds, or folds will be created locally.\n\n")
}

# If not loaded from cv_pipeline cache, try local cache or create new
if (!spatial_folds_loaded) {
  cat("Trying local cache or creating new folds...\n")
  local_cache_path <- file.path(out_dir, "spatial_folds_1km_cache.rds")
  
  if (file.exists(local_cache_path)) {
    cat("Loading from local cache...\n")
    spatial_cv <- tryCatch(readRDS(local_cache_path), error = function(e) NULL)
    if (!is.null(spatial_cv) && "folds_ids" %in% names(spatial_cv)) {
      if (length(spatial_cv$folds_ids) == nrow(dat_complete)) {
        folds_ids <- spatial_cv$folds_ids
        spatial_folds_loaded <- TRUE
        cat("  ✓ Loaded from local cache (", local_cache_path, ")\n", sep = "")
      } else {
        cat("  Local cache fold length (", length(spatial_cv$folds_ids), ") != data rows (", nrow(dat_complete), "), creating new folds.\n", sep = "")
      }
    }
  }
  
  if (!spatial_folds_loaded) {
    cat("Creating new spatial folds (this may take a moment)...\n")
    # Convert to sf object
    core_sf <- sf::st_as_sf(
      dat_complete,
      coords = c("longitude", "latitude"),
      crs = 4326
    )
    
    spatial_cv <- blockCV::cv_spatial(
      x = core_sf,
      k = n_folds,
      size = block_size_m,
      selection = "random",
      iteration = 50,
      progress = FALSE,
      biomod2 = FALSE,
      hexagon = FALSE,
      plot = FALSE
    )
    
    folds_ids <- spatial_cv$folds_ids
    
    # Save to local cache
    tryCatch({
      saveRDS(spatial_cv, local_cache_path)
      cat("  Saved to local cache:", local_cache_path, "\n")
    }, error = function(e) {
      cat("  Could not save local cache:", e$message, "\n")
    })
    
    cat("  Spatial folds created: ", length(unique(folds_ids)), " folds\n\n")
  }
}

if (is.null(folds_ids)) {
  stop("Failed to load or create spatial folds. Check data and blockCV installation.")
}

## ================================ HYPERPARAMETER TUNING ================================
cat("========================================\n")
cat("HYPERPARAMETER TUNING\n")
cat("========================================\n\n")

run_config <- list(
  target_var = target_var,
  block_size_km = block_size_km,
  n_folds = n_folds,
  kernels = sort(kernels_to_test),
  feature_sets_names = sort(names(feature_sets)),
  n_obs = nrow(dat_complete),
  n_fold_ids = length(unique(folds_ids))
)
cache_path <- file.path(out_dir, "gpr_tuning_results_cache.rds")
tuning_results <- list()
if (file.exists(cache_path)) {
  cached <- tryCatch(readRDS(cache_path), error = function(e) NULL)
  if (!is.null(cached) &&
      identical(cached$config$target_var, run_config$target_var) &&
      identical(cached$config$block_size_km, run_config$block_size_km) &&
      identical(cached$config$n_folds, run_config$n_folds) &&
      identical(cached$config$kernels, run_config$kernels) &&
      identical(cached$config$feature_sets_names, run_config$feature_sets_names) &&
      identical(cached$config$n_obs, run_config$n_obs)) {
    tuning_results <- cached$tuning_results
    cat("Loaded tuning cache:", length(tuning_results), "combinations already computed.\n")
  }
}

cat("Testing combinations of:\n")
cat("  - Kernels:", length(kernels_to_test), "options:", paste(kernels_to_test, collapse = ", "), "\n")
cat("  - Feature sets:", length(feature_sets), "options:", paste(names(feature_sets), collapse = ", "), "\n")
cat("  - Total combinations:", length(kernels_to_test) * length(feature_sets), "\n\n")

combination_num <- 0
total_combinations <- length(kernels_to_test) * length(feature_sets)

for (kernel in kernels_to_test) {
  for (feat_set_name in names(feature_sets)) {
    combination_num <- combination_num + 1
    result_key <- paste(kernel, feat_set_name, sep = "_")
    if (result_key %in% names(tuning_results)) {
      cat("\n--- Combination", combination_num, "/", total_combinations,
          ": Kernel =", kernel, ", Features =", feat_set_name, " (cached) ---\n")
      next
    }
    predictor_vars_subset <- feature_sets[[feat_set_name]]

    cat("\n--- Combination", combination_num, "/", total_combinations,
        ": Kernel =", kernel, ", Features =", feat_set_name, "(",
        length(predictor_vars_subset), "vars) ---\n")

    fold_results <- list()
    for (fold in 1:n_folds) {
      cat("  Fold", fold, "/", n_folds, "... ")
      train_idx <- folds_ids != fold
      test_idx <- folds_ids == fold
      train_data <- dat_complete[train_idx, ]
      test_data <- dat_complete[test_idx, ]
      result <- fit_gpr_cv(
        train_data = train_data,
        test_data = test_data,
        value_var = target_var,
        predictor_vars = predictor_vars_subset,
        kernel = kernel
      )
      fold_results[[fold]] <- result
      if (!is.null(result$error)) {
        cat("ERROR:", result$error, "\n")
      } else {
        cat("R² =", round(result$r2, 4), ", RMSE =", round(result$rmse, 4), "\n")
      }
    }

    r2_values <- sapply(fold_results, function(x) x$r2)
    rmse_values <- sapply(fold_results, function(x) x$rmse)
    valid_r2 <- r2_values[!is.na(r2_values)]
    valid_rmse <- rmse_values[!is.na(rmse_values)]

    if (length(valid_r2) > 0) {
      tuning_results[[result_key]] <- list(
        kernel = kernel,
        n_features = length(predictor_vars_subset),
        feature_set = feat_set_name,
        mean_r2 = mean(valid_r2),
        sd_r2 = sd(valid_r2),
        mean_rmse = mean(valid_rmse),
        sd_rmse = sd(valid_rmse),
        n_folds_valid = length(valid_r2),
        fold_r2 = r2_values,
        fold_rmse = rmse_values,
        features_used = predictor_vars_subset
      )
      cat("  Summary: Mean R² =", round(mean(valid_r2), 4),
          "±", round(sd(valid_r2), 4), "\n")
      cat("           Mean RMSE =", round(mean(valid_rmse), 4),
          "±", round(sd(valid_rmse), 4), "\n")
    } else {
      cat("  No valid results for this combination\n")
    }
    saveRDS(list(config = run_config, tuning_results = tuning_results), cache_path)
  }
}

## ================================ RESULTS SUMMARY ================================
cat("\n========================================\n")
cat("TUNING RESULTS SUMMARY\n")
cat("========================================\n\n")

if (length(tuning_results) == 0) {
  cat("ERROR: No valid results obtained. Check data and model fitting.\n")
} else {
  # Create results dataframe
  results_df <- do.call(rbind, lapply(tuning_results, function(x) {
    data.frame(
      kernel = x$kernel,
      n_features = x$n_features,
      feature_set = x$feature_set,
      mean_r2 = x$mean_r2,
      sd_r2 = x$sd_r2,
      mean_rmse = x$mean_rmse,
      sd_rmse = x$sd_rmse,
      n_folds_valid = x$n_folds_valid,
      stringsAsFactors = FALSE
    )
  }))
  
  # Sort by mean R²
  results_df <- results_df[order(-results_df$mean_r2), ]
  
  cat("Results (sorted by R²):\n")
  print(results_df)
  cat("\n")
  
  # Best configuration
  best_config <- results_df[1, ]
  best_result <- tuning_results[[paste(best_config$kernel, best_config$feature_set, sep = "_")]]
  
  cat("BEST CONFIGURATION:\n")
  cat("  Kernel:", best_config$kernel, "\n")
  cat("  Number of features:", best_config$n_features, "\n")
  cat("  Feature set:", best_config$feature_set, "\n")
  cat("  Mean R²:", round(best_config$mean_r2, 4), "±", round(best_config$sd_r2, 4), "\n")
  cat("  Mean RMSE:", round(best_config$mean_rmse, 4), "±", round(best_config$sd_rmse, 4), "\n")
  cat("\n  Features used:\n")
  cat("    ", paste(best_result$features_used, collapse = ", "), "\n\n")
  
  # Save results
  write.csv(results_df,
            file.path(out_dir, "gpr_hyperparameter_tuning_results.csv"),
            row.names = FALSE)
  cat("Results saved to: gpr_hyperparameter_tuning_results.csv\n\n")

  # Save best config for downstream scripts (gpr_fit_and_predict, etc.)
  best_config_list <- list(
    kernel = best_config$kernel,
    feature_set = best_config$feature_set,
    n_features = best_config$n_features,
    features_used = best_result$features_used
  )
  saveRDS(best_config_list, file.path(out_dir, "gpr_best_config.rds"))
  cat("Best config saved to: gpr_best_config.rds\n\n")

  # Plot comparison: Kernel × Feature set
  cat("Creating comparison plots...\n")
  
  # Plot 1: R² by kernel and feature set
  p_comparison <- ggplot(results_df, aes(x = factor(n_features), y = mean_r2, fill = kernel)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
    geom_errorbar(aes(ymin = mean_r2 - sd_r2, ymax = mean_r2 + sd_r2),
                  position = position_dodge(0.9), width = 0.2, color = "black") +
    labs(
      x = "Number of Features",
      y = "Mean R² (1km spatial CV)",
      fill = "Kernel",
      title = "GPR Hyperparameter Tuning: Kernel × Feature Selection"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(size = 12, face = "bold"),
          legend.position = "bottom")
  
  print(p_comparison)
  ggsave(file.path(out_dir, "gpr_hyperparameter_tuning_comparison.png"), 
         p_comparison, width = 12, height = 6)
  
  # Plot 2: Heatmap of R² by kernel and n_features
  p_heatmap <- ggplot(results_df, aes(x = factor(n_features), y = kernel, fill = mean_r2)) +
    geom_tile(color = "white", linewidth = 0.5) +
    scale_fill_viridis_c(name = "R²", option = "plasma") +
    labs(
      x = "Number of Features",
      y = "Kernel",
      title = "GPR Performance Heatmap: R² by Kernel and Feature Count"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(size = 12, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(p_heatmap)
  ggsave(file.path(out_dir, "gpr_hyperparameter_tuning_heatmap.png"), 
         p_heatmap, width = 10, height = 6)
  
  cat("Plots saved:\n")
  cat("  - gpr_hyperparameter_tuning_comparison.png\n")
  cat("  - gpr_hyperparameter_tuning_heatmap.png\n\n")
}

## ================================ ADVICE ON SPARSE DISTRIBUTIONS ================================
cat("========================================\n")
cat("ADVICE: SPARSE DISTRIBUTIONS & OVERFITTING\n")
cat("========================================\n\n")

cat("YOUR COVARIATE DISTRIBUTIONS:\n")
cat("----------------------------------------\n")
cat("Many covariates show sparse, gapped distributions:\n")
cat("  - bottomT_mean_daily_C: Clustered around 8-12°C and 17-19°C\n")
cat("  - No_sig_wave_heights_mean: Discrete peaks at 10, 25-30, 35-40\n")
cat("  - S3A_OLCI_ERRNT: Bimodal distribution with large gap\n")
cat("  - Surf_spco2_mean_uatm: Sparse peaks at 330, 360, 380, 400\n")
cat("  - wave_height_p95_m: Discrete values at 1.0, 1.5, 2.0, 2.5, 3.0, 3.5m\n\n")

cat("WHY THIS CAUSES JAGGED PARTIAL DEPENDENCE:\n")
cat("----------------------------------------\n")
cat("1. DATA GAPS:\n")
cat("   - Large regions with no observations\n")
cat("   - GPR must interpolate/extrapolate in these gaps\n")
cat("   - Model uncertainty is high in these regions\n")
cat("   - Predictions can vary wildly between nearby empty regions\n\n")

cat("2. OVERFITTING TO SPARSE CLUSTERS:\n")
cat("   - GPR fits closely to observed clusters\n")
cat("   - In gaps, model is unconstrained and can show high-frequency patterns\n")
cat("   - Matern 5/2 kernel allows rapid changes (high flexibility)\n")
cat("   - This creates the 'wiggly' appearance\n\n")

cat("3. CONDITIONAL EFFECTS:\n")
cat("   - Variables may only matter in specific ranges\n")
cat("   - Gaps hide the true relationship\n")
cat("   - Partial dependence averages over gaps, creating artifacts\n\n")

cat("RECOMMENDATIONS:\n")
cat("----------------------------------------\n")
cat("1. USE SMOOTHER KERNELS:\n")
cat("   - Matern 3/2 is smoother than Matern 5/2\n")
cat("   - Gaussian kernel is smoothest (but may underfit)\n")
cat("   - Test which gives best R² on spatial CV\n\n")

cat("2. REGULARIZATION:\n")
cat("   - GauPro handles regularization internally\n")
cat("   - Smoother kernels = implicit regularization\n")
cat("   - Consider reducing number of predictors (feature selection)\n\n")

cat("3. INTERPRETATION:\n")
cat("   - Focus on regions with data (not gaps)\n")
cat("   - Use smoothing for visualization (LOESS)\n")
cat("   - Check permutation importance (shows true importance)\n")
cat("   - Wide confidence intervals in gaps = high uncertainty\n\n")

cat("4. MODEL VALIDATION:\n")
cat("   - Good R² on spatial CV = model generalizes well\n")
cat("   - Jagged PDPs don't necessarily mean poor model\n")
cat("   - If R² is good, model is learning useful patterns\n")
cat("   - Jaggedness may reflect true conditional effects\n\n")

cat("5. DATA COLLECTION:\n")
cat("   - Consider collecting more data in gaps\n")
cat("   - Or acknowledge limitations in sparse regions\n")
cat("   - Focus interpretation on well-sampled ranges\n\n")

cat("\n========================================\n")
cat("TUNING COMPLETE\n")
cat("========================================\n\n")
