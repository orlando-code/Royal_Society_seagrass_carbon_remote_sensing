# Cross-Validation Pipeline: Random split vs spatial block CV
#
# Tests cross-fold random and spatial block CV; plots comparisons for each.
# Uses predefined spatial block sizes (1, 5, 25, 50, 100 km). No autocorrelation
# or distance-buffered CV (those remain commented out).

## ================================ SETUP ================================
# rm(list = ls())
setwd(here::here())
set.seed(42)

source("modelling/R/helpers.R")
source("modelling/R/plot_config.R")
source("modelling/R/extract_covariates_from_rasters.R")
load_packages(c(
  "here", "mgcv", "tidyverse", "ggplot2", "randomForest", "blockCV",
  "corrplot", "patchwork", "scales", "GauPro", "xgboost", "e1071", "neuralnet", "sf", "gstat"
))
# INLA optional (for INLA model in CV); install from https://www.r-inla.org if needed
if (requireNamespace("INLA", quietly = TRUE)) cat("INLA loaded for CV.\n") else cat("INLA not installed; INLA model will be skipped in CV.\n")

## ================================ CONFIG ================================
if (!exists("n_folds", envir = .GlobalEnv)) stop("n_folds must be set in the global environment!")
n_folds <- get("n_folds", envir = .GlobalEnv)
if (exists("dpi", envir = .GlobalEnv)) dpi <- get("dpi", envir = .GlobalEnv) else dpi <- 150
if (exists("show_titles", envir = .GlobalEnv)) show_titles <- get("show_titles", envir = .GlobalEnv) else show_titles <- TRUE
cat("NFOLDS:", n_folds, "\n")
use_hyperparameter_tuning <- FALSE
nested_tuning <- FALSE
# Use additional predictor sets based on external pruning/permutation importance (GAM/XGBoost/GPR).
# When FALSE, only the full predictor set ('all') is used in CV.
use_permutation_predictor_sets <- TRUE
# Set TRUE to run only random_split + first spatial block (no sensitivity) for a quick check/plot
quick_cv_check <- FALSE
# model_list <- c("GPR", "GAM", "XGB", "RF", "SVM", "NN", "OK", "UK", "INLA")
model_list <- c("GPR", "GAM", "XGB")

# Fixed spatial block sizes (km) for block CV – used for all spatial block strategies
# spatial_block_sizes_km <- c(1, 5, 10, 25, 50, 100)
spatial_block_sizes_km <- c(1, 5)
default_block_sizes_m <- spatial_block_sizes_km * 1000 # in metres

# Distance-buffered CV: hold out ~20% with test points > buffer_km from training (gap-filling test)
buffer_km <- 1
buffer_m <- buffer_km * 1000

# Threshold to flag variables with very large autocorrelation (excluded in sensitivity)
high_range_threshold_km <- 400
# Also exclude top N variables by range (so sensitivity always runs when we have autocorrelation data)
high_range_top_n <- 2L
# Maximum block size (m): blockCV fails with "Edge crosses" or "0 blocks" when blocks are too large
max_block_size_m <- 500000 # 500 km

out_dir <- "figures/cv_pipeline_output"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## ================================ DATA ================================
cat("\n========================================\n")
cat("LOADING CORE-LEVEL DATA\n")
cat("========================================\n\n")

dat <- read_rds("data/all_extracted_new.rds")
target_var <- "median_carbon_density_100cm"
# Predictors: raster-derived covariates that exist in dat (updated predictor set)
predictor_vars <- raster_covariates[raster_covariates %in% colnames(dat)]
if (length(predictor_vars) == 0) {
  stop("No predictor_vars found in dat. Check that all_extracted_new.rds contains columns matching raster config.")
}

# Build core-level data: one row per core with median_carbon_density (required by run_cv/fit_*)
meta_exclude <- c("Region", "number_id_final_version", "n_na_per_point", "seagrass_species", "random_core_variable", target_var, "longitude", "latitude")
if ("random_core_variable" %in% colnames(dat)) {
  core_data <- dat %>%
    dplyr::group_by(random_core_variable) %>%
    dplyr::summarise(
      median_carbon_density = median(.data[[target_var]], na.rm = TRUE),
      dplyr::across(dplyr::all_of(c(predictor_vars, "longitude", "latitude")), first),
      .groups = "drop"
    ) %>%
    dplyr::filter(!is.na(median_carbon_density))
} else {
  core_data <- dat %>%
    dplyr::mutate(median_carbon_density = .data[[target_var]]) %>%
    dplyr::select(longitude, latitude, median_carbon_density, dplyr::all_of(predictor_vars))
}

core_data_complete <- core_data %>%
  dplyr::select(longitude, latitude, median_carbon_density, dplyr::all_of(predictor_vars)) %>%
  dplyr::filter(complete.cases(.))
predictor_vars <- predictor_vars[predictor_vars %in% colnames(core_data_complete)]

n_cores <- nrow(core_data_complete)
cat("Cores (complete cases):", n_cores, "\n")
cat("Predictors (all):", length(predictor_vars), "\n")

# Load pruned predictor sets: GAM-based and XGBoost-based (used for additional CV runs)
pruned_file <- "figures/covariate_selection/pruned_variables_to_include_core_level.csv"
if (!file.exists(pruned_file)) pruned_file <- "figures/covariate_selection/pruned_variables_to_include.csv"
predictor_vars_pruned <- character(0)
if (file.exists(pruned_file)) {
  pruned_df <- read.csv(pruned_file, stringsAsFactors = FALSE)
  predictor_vars_pruned <- intersect(pruned_df$variable, colnames(core_data_complete))
  if (length(predictor_vars_pruned) < 2) {
    cat("GAM-pruned list has < 2 predictors in core data; skipping pruned CV.\n")
    predictor_vars_pruned <- character(0)
  } else {
    cat("Predictors (pruned, GAM):", length(predictor_vars_pruned), "\n")
  }
} else {
  cat("GAM-pruned variable list not found.\n")
}
# pruned_xgboost_file <- "figures/covariate_selection/pruned_variables_to_include_core_level_xgboost.csv"
# if (!file.exists(pruned_xgboost_file)) pruned_xgboost_file <- "figures/covariate_selection/pruned_variables_to_include_xgboost.csv"
# predictor_vars_pruned_xgboost <- character(0)
# if (file.exists(pruned_xgboost_file)) {
#   pruned_xgb_df <- read.csv(pruned_xgboost_file, stringsAsFactors = FALSE)
#   predictor_vars_pruned_xgboost <- intersect(pruned_xgb_df$variable, colnames(core_data_complete))
#   if (length(predictor_vars_pruned_xgboost) < 2) {
#     cat("XGBoost-pruned list has < 2 predictors in core data; skipping pruned_xgboost CV.\n")
#     predictor_vars_pruned_xgboost <- character(0)
#   } else {
#     cat("Predictors (pruned, XGBoost):", length(predictor_vars_pruned_xgboost), "\n")
#   }
# } else {
#   cat("XGBoost-pruned variable list not found; run covariate_pruning.R to generate it.\n")
# }
# # GPR pruned set from generic model permutation pruning (model_permutation_pruning.R)
# predictor_vars_pruned_gpr <- get_pruned_predictors_for_model("GPR", data_cols = colnames(core_data_complete))
# if (!is.null(predictor_vars_pruned_gpr) && length(predictor_vars_pruned_gpr) >= 2L) {
#   cat("Predictors (pruned, GPR from model permutation pruning):", length(predictor_vars_pruned_gpr), "\n")
# } else {
#   predictor_vars_pruned_gpr <- character(0)
#   if (!file.exists(file.path(out_dir, "cv_model_permutation_pruning_spatial_block_100000m.csv")) &&
#     length(list.files(out_dir, pattern = "cv_model_permutation_pruning_.*\\.csv")) == 0L) {
#     cat("GPR pruned set not found; run model_permutation_pruning.R to generate it.\n")
#   }
# }
# cat("\n")

block_sizes_m <- default_block_sizes_m
high_range_vars <- character(0)
autocor_note <- sprintf("Spatial block CV: fixed sizes %s km (autocorrelation not used).", paste(spatial_block_sizes_km, collapse = ", "))

cat("Block sizes (m) for spatial CV:", paste(block_sizes_m, collapse = ", "), "\n")
cat(autocor_note, "\n\n")

# Convert internal method name to display label (for plots; handles e.g. spatial_block_1e+05m)
method_to_display <- function(m) {
  if (grepl("^spatial_block_.+m$", m)) {
    size_str <- sub("^spatial_block_(.+)m$", "\\1", m)
    size_m <- as.numeric(size_str)
    size_km <- if (!is.na(size_m) && size_m >= 1000) round(size_m / 1000) else size_m
    if (!is.na(size_km)) {
      return(paste("Spatial block", size_km, "km"))
    }
  }
  if (grepl("^distance_buffer_", m)) {
    km <- sub("distance_buffer_|km", "", m)
    return(paste("Distance buffer", km, "km"))
  }
  switch(m,
    random_split = "Random split",
    env_cluster = "Env. cluster",
    paste0(toupper(substring(m, 1, 1)), substring(m, 2))
  )
}

## ================================ CV STRATEGIES ================================
cat("========================================\n")
cat("BUILDING CV STRATEGIES\n")
cat("========================================\n\n")

core_sf <- sf::st_as_sf(
  core_data_complete,
  coords = c("longitude", "latitude"),
  crs = 4326,
  remove = FALSE
)
light_core_sf <- sf::st_as_sf(
  core_data_complete %>% dplyr::select(-dplyr::all_of(predictor_vars)),
  coords = c("longitude", "latitude"),
  crs = 4326,
  remove = FALSE
) # create a version without all the variables to save memory when calculating spatial folds
cat("Core-level data:", nrow(core_data_complete), "cores,", ncol(core_data_complete), "variables.\n")

cat("\tCreating random folds\n")
# 1) Random split
random_folds <- sample(rep(seq_len(n_folds), length.out = n_cores))
cv_strategies <- list(
  list(method = "random_split", folds = random_folds, block_size_m = NA)
)


cat("\tCreating spatial folds\n")
# 4) Spatial block CV at each predefined block size (1, 5, 25, 50, 100 km) – use cache to avoid slow rebuilds
spatial_folds_cache_path <- file.path(out_dir, "spatial_folds_cache.rds")
cache_key <- list(n_cores = n_cores, block_sizes_m = block_sizes_m, n_folds = n_folds)
spatial_loaded_from_cache <- FALSE
if (file.exists(spatial_folds_cache_path)) {
  cached <- tryCatch(readRDS(spatial_folds_cache_path), error = function(e) NULL)
  if (!is.null(cached) && length(cached$spatial_strategies) > 0) {
    # Reuse cached strategies when present; warn if cache_key differs
    if (!identical(cached$cache_key$n_cores, cache_key$n_cores) ||
      !identical(cached$cache_key$block_sizes_m, cache_key$block_sizes_m) ||
      !identical(cached$cache_key$n_folds, cache_key$n_folds)) {
      cat("  Using cached spatial blocks (config differs from current).\n")
    }
    for (s in cached$spatial_strategies) {
      cv_strategies[[length(cv_strategies) + 1]] <- s
    }
    spatial_loaded_from_cache <- TRUE
    cat("  Loaded", length(cached$spatial_strategies), "spatial block strategies from cache.\n")
  }
}
if (!spatial_loaded_from_cache) {
  # get rid of all variables from 
  for (size_m in block_sizes_m) {
    cat("  Building spatial blocks ", size_m, " m ... ", sep = "")
    cv_spatial_out <- tryCatch(
      {
        blockCV::cv_spatial(
          x = light_core_sf[1:100,],
          k = n_folds,
          size = size_m,
          selection = "random",
          iteration = 50,
          progress = TRUE,
          biomod2 = FALSE,
          hexagon = FALSE,
          plot = FALSE
        ) # N.B. this is very resource intensive and slow.
      },
      error = function(e) {
        cat("failed:", e$message, "\n")
        return(NULL)
      }
    )
    if (!is.null(cv_spatial_out)) {
      cat("done.\n")
      fi <- cv_spatial_out$folds_ids
      n_actual_folds <- if (is.vector(fi) && !is.list(fi)) max(fi, na.rm = TRUE) else length(fi)
      if (n_actual_folds < 2) {
        cat("  spatial block size", size_m, "m: too few folds (", n_actual_folds, "), skipping.\n")
      } else {
        cat("  Added spatial_block_", size_m, "m (", n_actual_folds, " folds).\n", sep = "")
        cv_strategies[[length(cv_strategies) + 1]] <- list(
          method = paste0("spatial_block_", size_m, "m"),
          folds = cv_spatial_out$folds_ids,
          block_size_m = size_m,
          blocks = cv_spatial_out$blocks
        )
      }
    } else {
      cat("(no output).\n")
    }
  }
  # Save spatial strategies to cache for next run
  spatial_strategies <- cv_strategies[sapply(cv_strategies, function(s) grepl("^spatial_block_", s$method))]
  if (length(spatial_strategies) > 0) {
    tryCatch(
      {
        saveRDS(list(cache_key = cache_key, spatial_strategies = spatial_strategies), spatial_folds_cache_path)
        cat("  Saved spatial folds to cache:", spatial_folds_cache_path, "\n")
      },
      error = function(e) cat("  Could not save cache:", e$message, "\n")
    )
  }
}

if (exists("quick_cv_check") && isTRUE(quick_cv_check)) {
  idx_spatial <- which(sapply(cv_strategies, function(s) grepl("^spatial_block_", s$method)))
  keep <- if (length(idx_spatial) > 0) {
    c("random_split", cv_strategies[[idx_spatial[1]]]$method)
  } else {
    "random_split"
  }
  cv_strategies <- cv_strategies[sapply(cv_strategies, function(s) s$method %in% keep)]
  cat("Quick check: using only", length(cv_strategies), "CV strategies:", paste(keep, collapse = ", "), "\n")
}

cat("\nTotal CV strategies:", length(cv_strategies), "\n\n")

## ================================ FOLD ASSIGNMENTS AND ENV DISTRIBUTIONS ================================
# Build per-method fold vector (list -> vector for blockCV-style folds)
fold_vectors <- list()
for (i in seq_along(cv_strategies)) {
  s <- cv_strategies[[i]]
  if (is.list(s$folds)) {
    fv <- integer(n_cores)
    for (k in seq_along(s$folds)) fv[s$folds[[k]]] <- k
  } else {
    fv <- as.integer(s$folds)
  }
  fold_vectors[[s$method]] <- fv
}
# Long-format data: method, fold, and all predictor values
fold_env_list <- list()
for (m in names(fold_vectors)) {
  fold_env_list[[m]] <- core_data_complete %>%
    dplyr::mutate(method = m, fold = fold_vectors[[m]]) %>%
    dplyr::select(method, fold, dplyr::all_of(predictor_vars))
}
fold_env_df <- dplyr::bind_rows(fold_env_list)
fold_env_long <- fold_env_df %>%
  tidyr::pivot_longer(cols = dplyr::all_of(predictor_vars), names_to = "variable", values_to = "value") %>%
  mutate(method_label = vapply(method, method_to_display, character(1)))

# TODO: re-activate this
# # Plot: distributions of environmental variables across folds, by method
# cat("Building env distributions plot and saving ... ")
# p_env_by_fold <- ggplot(fold_env_long, aes(x = factor(fold), y = value, fill = factor(fold))) +
#   geom_boxplot(outlier.size = 0.6) +
#   facet_grid(variable ~ method_label, scales = "free_y") +
#   labs(x = "Fold", y = "Value", title = "Distribution of environmental variables across folds by CV method") +
#   theme_minimal() +
#   theme(
#     legend.position = "none",
#     axis.text.x = element_text(size = 8),
#     strip.text = element_text(size = 8),
#     plot.title = element_text(hjust = 0.5, size = 11)
#   )
# ggsave(file.path(out_dir, "cv_env_distributions_by_fold.png"), p_env_by_fold,
#   width = 2 + 2.5 * length(cv_strategies), height = 1.5 + 1.2 * length(predictor_vars),
#   dpi = 150, limitsize = FALSE
# )
# cat("Saved", file.path(out_dir, "cv_env_distributions_by_fold.png"), "\n")
# cat("Done.\n\n")

## ================================ RUN CV FOR EACH STRATEGY ================================
cat("========================================\n")
cat("RUNNING CROSS-VALIDATION\n")
cat("========================================\n\n")

predictor_sets <- list(all = predictor_vars)
# if (isTRUE(use_permutation_predictor_sets)) {
#   if (length(predictor_vars_pruned) >= 2) {
#     predictor_sets$pruned <- predictor_vars_pruned
#   }
#   if (length(predictor_vars_pruned_xgboost) >= 2) {
#     predictor_sets$pruned_xgboost <- predictor_vars_pruned_xgboost
#   }
#   if (length(predictor_vars_pruned_gpr) >= 2) {
#     predictor_sets$pruned_gpr <- predictor_vars_pruned_gpr
#   }
# }
n_ps <- length(predictor_sets)
cat("Running CV for ", n_ps, " predictor set(s): ", paste(names(predictor_sets), collapse = ", "), ".\n\n", sep = "")

all_results <- list()
all_predictions <- list()
n_runs <- length(cv_strategies) * length(predictor_sets)
run_idx <- 0
for (i in seq_along(cv_strategies)) {
  # TODO: run up to n_folds. Would this involve necessarily recalculating the fold assignment?
  s <- cv_strategies[[i]]
  buf_m <- if (!is.null(s$buffer_m)) s$buffer_m else NA_real_
  for (ps_name in names(predictor_sets)) {
    pvars <- predictor_sets[[ps_name]]
    run_idx <- run_idx + 1
    cat("  \n[", run_idx, "/", n_runs, "] ", s$method, ", predictor_set = ", ps_name, " ...\n", sep = "")
    res <- tryCatch(
      {
        run_cv(
          cv_method_name = s$method,
          fold_indices = s$folds,
          core_data = core_data_complete,
          predictor_vars = pvars,
          tune_hyperparams = use_hyperparameter_tuning,
          nested_tuning = nested_tuning,
          verbose = FALSE,
          buffer_m = buf_m,
          core_sf = core_sf,
          return_predictions = TRUE,
          models = model_list
        )
      },
      error = function(e) {
        cat("  CV failed for", s$method, " predictor_set=", ps_name, ": ", e$message, "\n", sep = "")
        return(NULL)
      }
    )
    if (is.null(res)) next
    # Handle list return (metrics + predictions) or legacy data frame
    if (is.list(res) && "metrics" %in% names(res)) {
      res_metrics <- res$metrics
      res_preds <- res$predictions
    } else {
      res_metrics <- res
      res_preds <- NULL
    }
    if (!is.null(res_metrics) && nrow(res_metrics) > 0) {
      res_metrics$block_size_m <- s$block_size_m
      res_metrics$predictor_set <- ps_name
      key <- paste0(s$method, "__", ps_name)
      all_results[[key]] <- res_metrics
    }
    if (!is.null(res_preds) && nrow(res_preds) > 0) {
      res_preds$block_size_m <- s$block_size_m
      res_preds$predictor_set <- ps_name
      all_predictions[[length(all_predictions) + 1]] <- res_preds
    }
  }
}

if (length(all_results) == 0) {
  stop("No CV results produced. Check data and fold generation.")
}

results_df <- dplyr::bind_rows(all_results)
predictions_df <- if (length(all_predictions) > 0) dplyr::bind_rows(all_predictions) else NULL
# Save predictions for plotting script (optional; only if return_predictions was TRUE)
if (!is.null(predictions_df) && nrow(predictions_df) > 0) {
  write.csv(predictions_df, file.path(out_dir, "cv_predictions.csv"), row.names = FALSE)
  cat("Saved predictions to", file.path(out_dir, "cv_predictions.csv"), "\n")
}

## ================================ SUMMARISE PERFORMANCE ================================
cat("\n========================================\n")
cat("SUMMARY: PERFORMANCE BY METHOD AND MODEL\n")
cat("========================================\n\n")

grp <- c("method", "model")
if ("predictor_set" %in% names(results_df)) grp <- c(grp, "predictor_set")

summary_by_method_model <- results_df %>%
  group_by(across(all_of(grp))) %>%
  summarise(
    n_folds = n(),
    mean_r2 = mean(r2, na.rm = TRUE),
    sd_r2 = sd(r2, na.rm = TRUE),
    mean_rmse = mean(rmse, na.rm = TRUE),
    sd_rmse = sd(rmse, na.rm = TRUE),
    mean_mae = mean(mae, na.rm = TRUE),
    sd_mae = sd(mae, na.rm = TRUE),
    mean_bias = mean(bias, na.rm = TRUE),
    sd_bias = sd(bias, na.rm = TRUE),
    .groups = "drop"
  )

print(summary_by_method_model)

## ================================ TABLES AND PLOTS ================================
safe_write_csv <- function(x, path) {
  tryCatch(
    {
      write.csv(x, path, row.names = FALSE)
      cat("Saved", path, "\n")
    },
    error = function(e) {
      cat("Warning: failed to write", path, ":", e$message, "\n")
    }
  )
}

safe_write_csv(results_df, file.path(out_dir, "cv_results_raw.csv"))
safe_write_csv(results_df, file.path(out_dir, "cv_results_detailed.csv"))
safe_write_csv(summary_by_method_model, file.path(out_dir, "cv_results_summary.csv"))
safe_write_csv(summary_by_method_model, file.path(out_dir, "cv_results_summary_by_predictor_set.csv"))


# Use only colours for models that appear in the data (keep order for legend)
models_present <- unique(results_df$model)
MODEL_COLOURS <- MODEL_COLOURS[names(MODEL_COLOURS) %in% models_present]

# Summary plot: R² and RMSE by method and model
metric_labeller <- c(r2 = "R²", rmse = "RMSE", mae = "MAE", bias = "Bias")
p_metrics <- results_df %>%
  dplyr::select(method, model, r2, rmse, mae, bias) %>%
  mutate(method_label = vapply(method, method_to_display, character(1))) %>%
  tidyr::pivot_longer(cols = c(r2, rmse, mae, bias), names_to = "metric", values_to = "value") %>%
  ggplot(aes(x = method_label, y = value, fill = model)) +
  geom_boxplot(outlier.size = 0.8, position = position_dodge(0.8), width = 0.7) +
  facet_wrap(~metric, scales = "free_y", ncol = 2, labeller = as_labeller(metric_labeller)) +
  scale_fill_manual(values = MODEL_COLOURS, na.value = "gray75") +
  labs(x = "CV method", y = NULL, title = if (show_titles) "Cross-validation performance by method and model" else NULL) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5)
  )
ggsave(file.path(out_dir, "cv_metrics_by_method_model.png"), p_metrics, width = 12, height = 8, dpi = dpi)
cat("\nSaved", file.path(out_dir, "cv_metrics_by_method_model.png"), "\n")

# Scatter: observed vs predicted for all methods (one panel per model)
if (!is.null(predictions_df) && nrow(predictions_df) > 0) {
  pred_all <- predictions_df %>% filter(if ("predictor_set" %in% names(predictions_df)) predictor_set == "all" else TRUE)
  if (nrow(pred_all) > 0) {
    pred_all <- pred_all %>% mutate(method_label = vapply(method, method_to_display, character(1)))
    p_scatter <- ggplot(pred_all, aes(x = predicted, y = observed)) +
      geom_point(alpha = 0.4, size = 1) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "red", linewidth = 0.5) +
      facet_wrap(~model, scales = "free", ncol = 3) +
      labs(x = "Predicted", y = "Observed", title = if (show_titles) "Observed vs predicted by model (all CV folds)" else NULL) +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5))
    ggsave(file.path(out_dir, "cv_observed_vs_predicted.png"), p_scatter, width = 12, height = 10, dpi = dpi)
    cat("Saved", file.path(out_dir, "cv_observed_vs_predicted.png"), "\n")
  }
}

# # Spatial GAM: predictions and errors side by side (inspired by interpolation.R)
# if (!is.null(predictions_df) && nrow(predictions_df) > 0) {
#   gam_pred <- predictions_df %>%
#     filter(model == "GAM", if ("predictor_set" %in% names(predictions_df)) predictor_set == "all" else TRUE) %>%
#     mutate(error = observed - predicted)
#   if (nrow(gam_pred) > 0) {
#     world <- ggplot2::map_data("world")
#     xlim <- range(gam_pred$longitude, na.rm = TRUE)
#     ylim <- range(gam_pred$latitude, na.rm = TRUE)
#     p_gam_pred <- ggplot() +
#       geom_polygon(data = world, aes(x = long, y = lat, group = group), fill = "#eeeeee", colour = "#a5a5a5") +
#       geom_point(data = gam_pred, aes(x = longitude, y = latitude, colour = predicted), size = 1.5, alpha = 0.8) +
#       scale_colour_viridis_c(option = "turbo", name = "Predicted") +
#       coord_cartesian(xlim = xlim, ylim = ylim) +
#       theme_minimal() +
#       labs(x = "Longitude", y = "Latitude", title = "GAM: spatial distribution of predictions")
#     p_gam_err <- ggplot() +
#       geom_polygon(data = world, aes(x = long, y = lat, group = group), fill = "#eeeeee", colour = "#a5a5a5") +
#       geom_point(data = gam_pred, aes(x = longitude, y = latitude, colour = error), size = 1.5, alpha = 0.8) +
#       scale_colour_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "Error (obs - pred)") +
#       coord_cartesian(xlim = xlim, ylim = ylim) +
#       theme_minimal() +
#       labs(x = "Longitude", y = "Latitude", title = "GAM: spatial distribution of errors")
#     p_gam_spatial <- p_gam_pred + p_gam_err + plot_layout(ncol = 2)
#     ggsave(file.path(out_dir, "cv_gam_spatial_predictions_and_errors.png"), p_gam_spatial, width = 12, height = 5, dpi = 150)
#     cat("Saved", file.path(out_dir, "cv_gam_spatial_predictions_and_errors.png"), "\n")
#   }
# }

# # Direct comparison: GAM vs ML methods under the same CV (mean R² by model and CV method)
# # Filter to full predictor set for clean comparison
# results_all <- results_df %>%
#   filter(if ("predictor_set" %in% names(results_df)) predictor_set == "all" else TRUE)
# summary_plot <- results_all %>%
#   group_by(method, model) %>%
#   summarise(mean_r2 = mean(r2, na.rm = TRUE), sd_r2 = sd(r2, na.rm = TRUE), .groups = "drop") %>%
#   mutate(method_label = vapply(method, method_to_display, character(1)))
# p_gam_vs_ml <- ggplot(summary_plot, aes(x = reorder(model, mean_r2), y = mean_r2, fill = model)) +
#   geom_col() +
#   geom_errorbar(aes(ymin = pmax(0, mean_r2 - sd_r2), ymax = pmin(1, mean_r2 + sd_r2)), width = 0.2, linewidth = 0.4) +
#   facet_wrap(~method_label, scales = "free_x", ncol = 2) +
#   scale_fill_manual(values = MODEL_COLOURS, na.value = "gray75") +
#   coord_flip() +
#   labs(x = NULL, y = "Mean R² (± SD across folds)", title = "GAM vs ML methods: same CV test") +
#   theme_minimal() +
#   theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
# ggsave(file.path(out_dir, "cv_gam_vs_ml_by_method.png"), p_gam_vs_ml, width = 10, height = 4 + 1.2 * length(unique(summary_plot$method)), dpi = 150, limitsize = FALSE)
# cat("Saved", file.path(out_dir, "cv_gam_vs_ml_by_method.png"), "\n")

# # Two-panel: Random split vs Spatial block (first available) – direct side-by-side
# spatial_methods <- unique(summary_plot$method[grepl("^spatial_block_", summary_plot$method)])
# compare_methods <- c("random_split", if (length(spatial_methods) > 0) spatial_methods[1] else NULL)
# if (length(compare_methods) >= 2) {
#   summary_compare <- summary_plot %>% filter(method %in% compare_methods)
#   summary_compare$method_label <- factor(
#     vapply(summary_compare$method, method_to_display, character(1)),
#     levels = c(method_to_display(compare_methods[1]), method_to_display(compare_methods[2]))
#   )
#   cv_method_colours <- setNames(c("#377eb8", "#e41a1c"), levels(summary_compare$method_label))
#   p_random_vs_spatial <- ggplot(summary_compare, aes(x = reorder(model, mean_r2), y = mean_r2, fill = method_label)) +
#     geom_col(position = position_dodge(0.9), width = 0.8) +
#     geom_errorbar(aes(ymin = pmax(0, mean_r2 - sd_r2), ymax = pmin(1, mean_r2 + sd_r2)),
#       position = position_dodge(0.9), width = 0.2, linewidth = 0.4
#     ) +
#     scale_fill_manual(values = cv_method_colours, name = NULL) +
#     coord_flip() +
#     labs(x = NULL, y = "Mean R² (± SD)", title = "Random split vs spatial block CV: GAM and ML methods") +
#     theme_minimal() +
#     theme(legend.title = element_blank(), legend.position = "bottom", plot.title = element_text(hjust = 0.5))
#   ggsave(file.path(out_dir, "cv_random_vs_spatial_gam_ml.png"), p_random_vs_spatial, width = 9, height = 5, dpi = 150)
#   cat("Saved", file.path(out_dir, "cv_random_vs_spatial_gam_ml.png"), "\n")
# }

# Bar summary: mean R² and mean RMSE by method (averaged across models)
grp_fold <- c("method", "fold")
if ("predictor_set" %in% names(results_df)) grp_fold <- c(grp_fold, "predictor_set")
grp_method <- c("method")
if ("predictor_set" %in% names(results_df)) grp_method <- c(grp_method, "predictor_set")

summary_method <- results_df %>%
  group_by(across(all_of(grp_fold))) %>%
  summarise(
    mean_r2 = mean(r2, na.rm = TRUE),
    mean_rmse = mean(rmse, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(across(all_of(grp_method))) %>%
  summarise(
    r2 = mean(mean_r2, na.rm = TRUE),
    sd_r2 = sd(mean_r2, na.rm = TRUE),
    rmse = mean(mean_rmse, na.rm = TRUE),
    sd_rmse = sd(mean_rmse, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(method_label = vapply(method, method_to_display, character(1)))

p_method <- summary_method %>%
  tidyr::pivot_longer(cols = c(r2, rmse), names_to = "metric", values_to = "value") %>%
  ggplot(aes(x = reorder(method_label, -value), y = value, fill = metric)) +
  geom_col(position = "dodge") +
  labs(x = "CV method", y = "Mean (across folds and models)", title = if (show_titles) "Ensemble performance by CV method" else NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
if ("predictor_set" %in% names(summary_method)) {
  p_method <- p_method + facet_wrap(~predictor_set, ncol = 1)
}
ggsave(file.path(out_dir, "cv_ensemble_by_method.png"), p_method, width = 10, height = 5, dpi = dpi)

## ================================ MODEL SELECTION FOR GAP-FILLING ================================
# Goal: find a model with good performance at predicting carbon density at points *near* training data,
# validated by train-test splits. Most relevant CV strategies:
#   - random_split: random 80/20 (generalisation; test points are "near" in the sense of the full sample)
#   - distance_buffer_1km: test points >1 km from train (strict "fill gaps between samples" test)
# Rank models by mean R² (and mean RMSE) on these strategies; prefer distance-buffered for gap-filling.
gap_cv_methods <- c("random_split", paste0("distance_buffer_", buffer_km, "km"))
summary_gap <- summary_by_method_model %>%
  filter(
    method %in% gap_cv_methods,
    if ("predictor_set" %in% names(summary_by_method_model)) predictor_set == "all" else TRUE
  )
if (nrow(summary_gap) > 0) {
  # Rank by mean R² (descending) then by mean RMSE (ascending) within each CV method
  summary_gap <- summary_gap %>%
    group_by(method) %>%
    mutate(
      rank_r2 = row_number(desc(mean_r2)),
      rank_rmse = row_number(mean_rmse)
    ) %>%
    ungroup()
  # Overall rank: average rank across the two gap-filling CV methods (lower = better)
  rank_overall <- summary_gap %>%
    group_by(model) %>%
    summarise(
      mean_r2_avg = mean(mean_r2, na.rm = TRUE),
      mean_rmse_avg = mean(mean_rmse, na.rm = TRUE),
      rank_r2_avg = mean(rank_r2, na.rm = TRUE),
      rank_rmse_avg = mean(rank_rmse, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(rank_overall = (rank_r2_avg + rank_rmse_avg) / 2) %>%
    arrange(rank_overall)
  write.csv(rank_overall, file.path(out_dir, "cv_best_models_gap_filling.csv"), row.names = FALSE)
  write.csv(summary_gap, file.path(out_dir, "cv_gap_filling_by_method_model.csv"), row.names = FALSE)
  cat("\n========================================\n")
  cat("MODEL SELECTION FOR GAP-FILLING (predict near training data)\n")
  cat("========================================\n\n")
  cat("CV strategies used: ", paste(gap_cv_methods, collapse = ", "), "\n", sep = "")
  cat("(random_split = generalisation; distance_buffer = test points >1 km from train)\n\n")
  cat("Model ranking (by average rank of R² and RMSE across these CV methods):\n")
  print(rank_overall %>% select(model, mean_r2_avg, mean_rmse_avg, rank_overall))
  cat("\nTop 3 models for gap-filling:\n")
  top3 <- head(rank_overall, 3)
  for (i in seq_len(nrow(top3))) {
    cat("  ", i, ". ", top3$model[i],
      " (mean R² = ", round(top3$mean_r2_avg[i], 4),
      ", mean RMSE = ", round(top3$mean_rmse_avg[i], 6), ")\n",
      sep = ""
    )
  }
  cat("\nSaved: cv_best_models_gap_filling.csv, cv_gap_filling_by_method_model.csv\n")
} else {
  cat("\nNo results for gap-filling CV methods; skipping model ranking.\n")
}

## ================================ FINAL GRAPHS ================================
max_block_m <- max(block_sizes_m, na.rm = TRUE)

# 1) Effect of excluding very large autocorrelation variables (all vs without_high_range)
if ("predictor_set" %in% names(results_df) && "without_high_range" %in% results_df$predictor_set) {
  sum_exclude <- summary_by_method_model %>%
    group_by(model, predictor_set) %>%
    summarise(
      mean_r2 = mean(mean_r2, na.rm = TRUE),
      sd_r2 = sqrt(mean(sd_r2^2, na.rm = TRUE)),
      mean_rmse = mean(mean_rmse, na.rm = TRUE),
      sd_rmse = sqrt(mean(sd_rmse^2, na.rm = TRUE)),
      .groups = "drop"
    )
  p_exclude_r2 <- ggplot(sum_exclude, aes(x = predictor_set, y = mean_r2, fill = predictor_set)) +
    geom_col(position = "dodge", show.legend = FALSE) +
    geom_errorbar(aes(ymin = mean_r2 - sd_r2, ymax = mean_r2 + sd_r2), width = 0.2, linewidth = 0.5) +
    facet_wrap(~model, ncol = 3) +
    labs(x = "Predictor set", y = "Mean R² (± SD)", title = if (show_titles) "Effect of excluding high-autocorrelation variables" else NULL) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  p_exclude_rmse <- ggplot(sum_exclude, aes(x = predictor_set, y = mean_rmse, fill = predictor_set)) +
    geom_col(position = "dodge", show.legend = FALSE) +
    geom_errorbar(aes(ymin = pmax(0, mean_rmse - sd_rmse), ymax = mean_rmse + sd_rmse), width = 0.2, linewidth = 0.5) +
    facet_wrap(~model, ncol = 3) +
    labs(x = "Predictor set", y = "Mean RMSE (± SD)", title = if (show_titles) "Effect of excluding high-autocorrelation variables" else NULL) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  p_exclude_combined <- p_exclude_r2 + p_exclude_rmse + plot_layout(ncol = 1)
  ggsave(file.path(out_dir, "cv_effect_excluding_high_autocorrelation.png"), p_exclude_combined, width = 10, height = 8, dpi = dpi)
  cat("Saved", file.path(out_dir, "cv_effect_excluding_high_autocorrelation.png"), "\n")
}

# 2) Performance over block sizes (1, 5, 25, 50, 100 km)
spatial_results <- results_df %>%
  filter(grepl("^spatial_block_", method), is.na(block_size_m) == FALSE)
if ("predictor_set" %in% names(spatial_results)) {
  spatial_results <- spatial_results %>% filter(predictor_set == "all")
}
if (nrow(spatial_results) > 0 && "block_size_m" %in% names(spatial_results)) {
  spatial_results <- spatial_results %>%
    mutate(block_size_km = round(block_size_m / 1000))
  sum_block <- spatial_results %>%
    group_by(block_size_km, model) %>%
    summarise(mean_r2 = mean(r2, na.rm = TRUE), mean_rmse = mean(rmse, na.rm = TRUE), .groups = "drop")
  km_levels <- sort(unique(sum_block$block_size_km))
  sum_block <- sum_block %>%
    mutate(block_km = factor(block_size_km, levels = km_levels))
  p_block_r2 <- ggplot(sum_block, aes(x = block_km, y = mean_r2, colour = model, group = model)) +
    geom_line(linewidth = 1) +
    geom_point(size = 3) +
    scale_colour_manual(values = MODEL_COLOURS, na.value = "gray75") +
    labs(x = "Block size (km)", y = "Mean R²") +
    scale_x_discrete(labels = paste0(km_levels, " km")) +
    theme_minimal() +
    theme(legend.position = "bottom", plot.title = element_text(size = 10, hjust = 0.5))
  p_block_rmse <- ggplot(sum_block, aes(x = block_km, y = mean_rmse, colour = model, group = model)) +
    geom_line(linewidth = 1) +
    geom_point(size = 3) +
    scale_colour_manual(values = MODEL_COLOURS, na.value = "gray75") +
    labs(x = "Block size (km)", y = "Mean RMSE") +
    scale_x_discrete(labels = paste0(km_levels, " km")) +
    theme_minimal() +
    theme(legend.position = "bottom", plot.title = element_text(size = 10, hjust = 0.5))
  p_block_combined <- p_block_r2 + p_block_rmse + plot_layout(ncol = 1) + plot_annotation(title = if (show_titles) "Performance over spatial block sizes" else NULL)
  ggsave(file.path(out_dir, "cv_performance_over_block_sizes.png"), p_block_combined, width = 9, height = 7, dpi = dpi)
  cat("Saved", file.path(out_dir, "cv_performance_over_block_sizes.png"), "\n")
}

# 2b) Model selection for spatial prediction with uncertainty: highlight models that provide prediction intervals
# Compare key CV methods (random, 1 km, 5 km blocks); GAM and GPR are the only models that supply prediction intervals
key_methods <- c("random_split", "spatial_block_1000m", "spatial_block_5000m")
summary_unc <- summary_by_method_model %>%
  filter(method %in% key_methods)
if ("predictor_set" %in% names(summary_unc)) summary_unc <- summary_unc %>% filter(predictor_set == "all")
if (nrow(summary_unc) > 0) {
  summary_unc <- summary_unc %>%
    mutate(method_label = vapply(method, method_to_display, character(1)))
  method_order <- vapply(key_methods, method_to_display, character(1))
  summary_unc$method_label <- factor(summary_unc$method_label, levels = method_order)
  p_unc_simple <- ggplot(summary_unc, aes(x = method_label, y = mean_r2, fill = model)) +
    geom_col(position = position_dodge(0.85), width = 0.75) +
    geom_errorbar(aes(ymin = pmax(0, mean_r2 - sd_r2), ymax = pmin(1, mean_r2 + sd_r2)),
      position = position_dodge(0.85), width = 0.2, linewidth = 0.4
    ) +
    scale_fill_manual(values = MODEL_COLOURS, na.value = "gray75", name = "Model") +
    labs(
      x = "CV method",
      y = "Mean R² (± SD)",
      title = if (show_titles) "Model selection for spatial prediction with uncertainty" else NULL,
      subtitle = if (show_titles) "Kernel-based methods and GPR provide prediction intervals; others give point predictions only. Spatial block CV tests transfer to nearby, unsampled areas." else NULL
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "bottom", plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(size = 8))
  ggsave(file.path(out_dir, "cv_uncertainty_capable_models.png"), p_unc_simple, width = 8, height = 5, dpi = dpi)
  cat("Saved", file.path(out_dir, "cv_uncertainty_capable_models.png"), "\n")
}

# # 2c) All vs pruned predictors: compare R² by model and predictor set for key CV methods
# predictor_set_labels <- c(all = "All predictors", pruned = "Pruned (GAM)", pruned_xgboost = "Pruned (XGBoost)")
# if ("predictor_set" %in% names(summary_by_method_model) && length(unique(summary_by_method_model$predictor_set)) >= 2) {
#   key_methods_ps <- c("random_split", "spatial_block_1000m", "spatial_block_5000m")
#   summary_ps <- summary_by_method_model %>% filter(method %in% key_methods_ps)
#   if (nrow(summary_ps) > 0) {
#     summary_ps <- summary_ps %>%
#       mutate(
#         method_label = vapply(method, method_to_display, character(1)),
#         predictor_set_label = ifelse(predictor_set %in% names(predictor_set_labels),
#           predictor_set_labels[predictor_set], paste0("Pruned (", predictor_set, ")")
#         )
#       )
#     method_order_ps <- vapply(key_methods_ps, method_to_display, character(1))
#     summary_ps$method_label <- factor(summary_ps$method_label, levels = method_order_ps)
#     n_ps_facets <- length(unique(summary_ps$predictor_set_label))
#     p_all_vs_pruned <- ggplot(summary_ps, aes(x = method_label, y = mean_r2, fill = model)) +
#       geom_col(position = position_dodge(0.85), width = 0.75) +
#       geom_errorbar(aes(ymin = pmax(0, mean_r2 - sd_r2), ymax = pmin(1, mean_r2 + sd_r2)),
#         position = position_dodge(0.85), width = 0.2, linewidth = 0.4
#       ) +
#       scale_fill_manual(values = MODEL_COLOURS, na.value = "gray75", name = "Model") +
#       facet_wrap(~predictor_set_label, ncol = min(3, n_ps_facets)) +
#       labs(
#         x = "CV method",
#         y = "Mean R² (± SD)",
#         title = "All vs pruned predictor sets (GAM and XGBoost)",
#         subtitle = "Same folds; pruned sets reduce overfitting risk."
#       ) +
#       theme_minimal() +
#       theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "bottom", plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(size = 8))
#     ggsave(file.path(out_dir, "cv_all_vs_pruned_predictors.png"), p_all_vs_pruned, width = 4 + 3 * n_ps_facets, height = 5.5, dpi = 150)
#     cat("Saved", file.path(out_dir, "cv_all_vs_pruned_predictors.png"), "\n")
#   }
# }

# 3) Illustration: spatial blocks at 25%, 50%, 75%, 100% superposed on scatter of data points (one fold)
# Build list of (block_size_m, blocks_sf, folds) for each spatial block size; use stored blocks or re-run cv_spatial
block_illustration_data <- list()
for (idx in seq_along(cv_strategies)) {
  s <- cv_strategies[[idx]]
  if (is.na(s$block_size_m)) next
  blocks_sf <- s$blocks
  if (is.null(blocks_sf)) {
    cv_out <- tryCatch(
      {
        blockCV::cv_spatial(
          x = core_sf, k = n_folds, size = s$block_size_m, selection = "random",
          iteration = 50, progress = FALSE, biomod2 = FALSE, hexagon = FALSE, plot = FALSE
        )
      },
      error = function(e) NULL
    )
    if (!is.null(cv_out) && !is.null(cv_out$blocks)) {
      blocks_sf <- cv_out$blocks
      fold_list <- cv_out$folds_ids
    } else {
      next
    }
  } else {
    fold_list <- s$folds
  }
  block_illustration_data[[length(block_illustration_data) + 1]] <- list(
    block_size_m = s$block_size_m,
    blocks_sf = blocks_sf,
    fold_list = fold_list
  )
}
if (length(block_illustration_data) > 0) {
  highlight_fold <- 1L
  plot_list <- list()
  for (idx in seq_along(block_illustration_data)) {
    bl <- block_illustration_data[[idx]]
    blocks_sf <- bl$blocks_sf
    fold_list <- bl$fold_list
    fold_vec <- integer(n_cores)
    for (k in seq_along(fold_list)) fold_vec[fold_list[[k]]] <- k
    points_with_fold <- core_sf
    points_with_fold$fold <- fold_vec
    points_with_fold$is_test <- points_with_fold$fold == highlight_fold
    block_km <- round(bl$block_size_m / 1000)
    title <- paste0("Spatial block ", block_km, " km")
    fold_col <- if ("folds" %in% names(blocks_sf)) "folds" else if ("fold" %in% names(blocks_sf)) "fold" else setdiff(names(blocks_sf), "geometry")[1]
    if (is.na(fold_col) || length(fold_col) == 0) fold_col <- "folds"
    p_block <- ggplot() +
      geom_sf(data = blocks_sf, aes(fill = as.factor(.data[[fold_col]])), alpha = 0.2, colour = "gray50", linewidth = 0.4) +
      geom_sf(data = points_with_fold, aes(color = is_test), size = 1.5, alpha = 0.95, stroke = 0) +
      scale_fill_viridis_d(option = "E", name = "Fold", na.value = NA) +
      scale_color_manual(values = c("FALSE" = "gray20", "TRUE" = "red"), labels = c("Train", "Test (fold 1)"), name = NULL) +
      labs(title = if (show_titles) title else NULL, x = "Longitude", y = "Latitude") +
      theme_minimal() +
      theme(legend.position = "right", plot.title = element_text(size = 11, hjust = 0.5))
    plot_list[[length(plot_list) + 1]] <- p_block
  }
  n_panels <- length(plot_list)
  p_blocks_illustration <- patchwork::wrap_plots(plot_list, ncol = 2, byrow = TRUE) +
    plot_annotation(
      title = if (show_titles) "Spatial blocks (1, 5, 25, 50, 100 km) superposed on data points (fold 1 = test)" else NULL,
      theme = theme(plot.title = element_text(hjust = 0.5, size = 12))
    )
  ggsave(file.path(out_dir, "cv_spatial_blocks_and_points.png"), p_blocks_illustration,
    width = 10, height = max(5, 2.5 * ceiling(n_panels / 2)), dpi = dpi
  )
  cat("Saved", file.path(out_dir, "cv_spatial_blocks_and_points.png"), "\n")
}

## ================================ NOTES ================================
notes_file <- file.path(out_dir, "pipeline_notes.txt")
sink(notes_file)
cat("CV pipeline: random split vs spatial block CV\n")
cat("============================================\n\n")
cat(autocor_note, "\n")
cat("Block sizes (m):", paste(block_sizes_m, collapse = ", "), "\n")
cat("CV strategies: Random split + spatial blocks over specified block sizes (m).\n")
cat("Predictors: from raster config (build_covariate_config_from_dir), present in all_extracted_new.rds.\n")
cat("Response: median_carbon_density (from median_carbon_density_100cm, one value per core).\n")
cat("\nPlots: cv_metrics_by_method_model.png, cv_observed_vs_predicted.png\n")
cat("cv_performance_over_block_sizes.png,\n")
cat("cv_uncertainty_capable_models.png, cv_env_distributions_by_fold.png, cv_spatial_blocks_and_points.png.\n")
cat("\nAutocorrelation, env_cluster, and distance-buffered CV are commented out.\n")
sink(NULL)
cat("Notes written to", notes_file, "\n")

cat("\n========================================\n")
cat("PIPELINE COMPLETE\n")
cat("========================================\n")
