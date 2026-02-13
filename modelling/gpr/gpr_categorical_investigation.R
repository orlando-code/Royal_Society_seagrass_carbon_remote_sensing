# GPR paper: investigate inclusion of seagrass_species and Region as categorical predictors
#
# Runs 1km spatial CV for four predictor configurations (env only; +species; +region; +both)
# and writes a comparison table. Use the results to decide whether to include categoricals
# in the main model.

rm(list = ls())
setwd(here::here())
set.seed(42)

source("modelling/R/helpers.R")
source("modelling/R/extract_covariates_from_rasters.R")
source("modelling/R/gpr_funs.R")
load_packages(c("here", "tidyverse", "GauPro", "blockCV", "sf", "dplyr"))

target_var <- "median_carbon_density_100cm"
out_dir <- "figures/cv_pipeline_output"
n_folds <- 5
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("\n========================================\n")
cat("GPR CATEGORICAL INVESTIGATION\n")
cat("(species and region as predictors)\n")
cat("========================================\n\n")

# Load data
dat <- readr::read_rds("data/all_extracted_new.rds")
dat <- process_rs_covariates(dat)

# Pruned env vars: prefer generic model permutation pruning, then GPR-specific, then GAM list
env_vars <- get_pruned_predictors_for_model("GPR", data_cols = names(dat))
if (is.null(env_vars) || length(env_vars) < 2L) {
  pruned_file_gpr <- "figures/covariate_selection/pruned_variables_to_include_gpr.csv"
  pruned_file <- "figures/covariate_selection/pruned_variables_to_include.csv"
  if (file.exists(pruned_file_gpr)) {
    env_vars <- readr::read_csv(pruned_file_gpr, show_col_types = FALSE)$variable
  } else if (file.exists(pruned_file)) {
    env_vars <- readr::read_csv(pruned_file, show_col_types = FALSE)$variable
  } else {
    stop("Pruned variables not found. Run model_permutation_pruning.R or gpr_covariate_pruning.R first.")
  }
}
env_vars <- env_vars[env_vars %in% names(dat)]

kernel <- "matern52"
best_config_path <- file.path(out_dir, "gpr_best_config.rds")
if (file.exists(best_config_path)) {
  best <- readRDS(best_config_path)
  kernel <- best$kernel
  if (length(best$features_used) <= length(env_vars)) {
    env_vars <- intersect(best$features_used, names(dat))
  }
}

need_vars <- c("longitude", "latitude", target_var, env_vars)
has_species <- "seagrass_species" %in% names(dat)
has_region <- "Region" %in% names(dat)
if (has_species) need_vars <- c(need_vars, "seagrass_species")
if (has_region) need_vars <- c(need_vars, "Region")

dat_complete <- dat[, need_vars]
dat_complete <- dat_complete[complete.cases(dat_complete), ]

cat("Complete cases:", nrow(dat_complete), "\n")
cat("Environment predictors:", length(env_vars), "\n")
cat("seagrass_species in data:", has_species, "\n")
cat("Region in data:", has_region, "\n\n")

# Load 1km spatial folds (same logic as gpr_hyperparameter_tuning.R)
folds_ids <- NULL
cv_pipeline_cache_path <- NULL
for (base in c("figures/cv_pipeline_output", "modelling/cv_pipeline_output")) {
  candidate <- file.path(base, "spatial_folds_cache.rds")
  if (file.exists(candidate)) {
    cv_pipeline_cache_path <- candidate
    break
  }
}
if (!is.null(cv_pipeline_cache_path)) {
  cached <- tryCatch(readRDS(cv_pipeline_cache_path), error = function(e) NULL)
  if (!is.null(cached) && "spatial_strategies" %in% names(cached)) {
    matching_strategy <- NULL
    for (strategy in cached$spatial_strategies) {
      if (strategy$method == "spatial_block_1000m") {
        matching_strategy <- strategy
        break
      }
    }
    if (!is.null(matching_strategy)) {
      if (length(matching_strategy$folds) == nrow(dat_complete)) {
        folds_ids <- matching_strategy$folds
      } else if ("random_core_variable" %in% names(dat)) {
        dat_full <- readr::read_rds("data/all_extracted_new.rds")
        dat_full <- process_rs_covariates(dat_full)
        core_data_cv <- dat_full %>%
          dplyr::group_by(random_core_variable) %>%
          dplyr::summarise(longitude = dplyr::first(longitude), latitude = dplyr::first(latitude)) %>%
          dplyr::filter(!is.na(longitude) & !is.na(latitude)) %>%
          dplyr::arrange(random_core_variable)
        if (nrow(core_data_cv) == length(matching_strategy$folds)) {
          core_locations <- paste(round(core_data_cv$longitude, 8), round(core_data_cv$latitude, 8), sep = "_")
          fold_lookup <- data.frame(location = core_locations, fold = matching_strategy$folds,
                                   longitude = core_data_cv$longitude, latitude = core_data_cv$latitude, stringsAsFactors = FALSE)
          point_locations <- paste(round(dat_complete$longitude, 8), round(dat_complete$latitude, 8), sep = "_")
          folds_ids <- fold_lookup$fold[match(point_locations, fold_lookup$location)]
          unmatched <- is.na(folds_ids)
          if (sum(unmatched) > 0) {
            for (i in which(unmatched)) {
              lon <- dat_complete$longitude[i]; lat <- dat_complete$latitude[i]
              dists <- sqrt((fold_lookup$longitude - lon)^2 + (fold_lookup$latitude - lat)^2)
              folds_ids[i] <- fold_lookup$fold[which.min(dists)]
            }
          }
        }
      }
    }
  }
}
if (is.null(folds_ids)) {
  for (base in c("modelling/cv_pipeline_output", "figures/cv_pipeline_output")) {
    local_1km_path <- file.path(base, "spatial_folds_1km_cache.rds")
    if (file.exists(local_1km_path)) {
      spatial_cv <- tryCatch(readRDS(local_1km_path), error = function(e) NULL)
      if (!is.null(spatial_cv) && "folds_ids" %in% names(spatial_cv) && length(spatial_cv$folds_ids) == nrow(dat_complete)) {
        folds_ids <- spatial_cv$folds_ids
        break
      }
    }
  }
}
if (is.null(folds_ids)) {
  stop("1km spatial folds not found. Run gpr_hyperparameter_tuning.R first (it loads or creates folds and may save to spatial_folds_cache.rds or spatial_folds_1km_cache.rds).")
}
cat("Using", length(unique(folds_ids)), "spatial folds.\n\n")

# Define configurations (only include categoricals that exist)
configs <- list(
  env_only = env_vars
)
if (has_species) configs$env_plus_species <- c(env_vars, "seagrass_species")
if (has_region) configs$env_plus_region <- c(env_vars, "Region")
if (has_species && has_region) configs$env_plus_both <- c(env_vars, "seagrass_species", "Region")

results_list <- list()
for (config_name in names(configs)) {
  predictor_vars <- configs[[config_name]]
  cat("---", config_name, "(", length(predictor_vars), "predictors ) ---\n")
  r2_vec <- numeric(n_folds)
  rmse_vec <- numeric(n_folds)
  for (fold in 1:n_folds) {
    train_data <- dat_complete[folds_ids != fold, ]
    test_data <- dat_complete[folds_ids == fold, ]
    res <- fit_gpr_cv(
      train_data = train_data,
      test_data = test_data,
      value_var = target_var,
      predictor_vars = predictor_vars,
      kernel = kernel
    )
    r2_vec[fold] <- res$r2
    rmse_vec[fold] <- res$rmse
    if (!is.null(res$error)) cat("  Fold", fold, ":", res$error, "\n")
  }
  valid <- !is.na(r2_vec)
  results_list[[config_name]] <- data.frame(
    configuration = config_name,
    mean_r2 = mean(r2_vec[valid]),
    sd_r2 = if (sum(valid) > 1) sd(r2_vec[valid]) else NA_real_,
    mean_rmse = mean(rmse_vec[valid]),
    sd_rmse = if (sum(valid) > 1) sd(rmse_vec[valid]) else NA_real_,
    n_folds_valid = sum(valid),
    stringsAsFactors = FALSE
  )
  cat("  Mean R² =", round(mean(r2_vec[valid]), 4),
      ", Mean RMSE =", round(mean(rmse_vec[valid]), 4), "\n\n")
}

results_df <- do.call(rbind, results_list)
write.csv(results_df, file.path(out_dir, "gpr_categorical_investigation.csv"), row.names = FALSE)
cat("Saved:", file.path(out_dir, "gpr_categorical_investigation.csv"), "\n")

best_r2 <- which.max(results_df$mean_r2)
cat("\nBest configuration by mean R²:", results_df$configuration[best_r2],
    "(R² =", round(results_df$mean_r2[best_r2], 4), ")\n")
cat("Use this to set include_species / include_region in gpr_predictions.R if desired.\n")
