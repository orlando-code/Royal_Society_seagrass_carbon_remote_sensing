# Cross-Validation Pipeline: Random split vs spatial block CV
#
# Tests cross-fold random and spatial block CV; plots comparisons for each.
# Uses predefined spatial block sizes (1, 5, 25, 50, 100 km). No autocorrelation
# or distance-buffered CV (those remain commented out).

## ================================ SETUP ================================
# TODO: ensure specific included for all models
setwd(here::here())
set.seed(42)

source("modelling/R/helpers.R")
source("modelling/R/assign_region_from_latlon.R")
source("modelling/R/plot_config.R")
source("modelling/R/extract_covariates_from_rasters.R")
load_packages(c(
  "here", "mgcv", "tidyverse", "ggplot2", "randomForest", "blockCV",
  "corrplot", "patchwork", "scales", "GauPro", "xgboost", "sf"
))

## ================================ CONFIG ================================
n_folds <- get0("n_folds", envir = .GlobalEnv, ifnotfound = 5)
if (is.null(n_folds)) stop("n_folds must be set in the global environment.")
cv_type   <- get0("cv_type",   envir = .GlobalEnv, ifnotfound = "spatial")
cv_blocksize <- get0("cv_blocksize", envir = .GlobalEnv, ifnotfound = 5000L)
dpi         <- get0("dpi",         envir = .GlobalEnv, ifnotfound = 150)
show_titles <- get0("show_titles", envir = .GlobalEnv, ifnotfound = TRUE)
cat("NFOLDS:", n_folds, " | CV type:", cv_type, "\n")
use_hyperparameter_tuning <- FALSE
nested_tuning <- FALSE
quick_cv_check <- get0("quick_cv_check", envir = .GlobalEnv, ifnotfound = FALSE)
model_list <- get0("model_list", envir = .GlobalEnv, ifnotfound = c("GPR", "GAM", "XGB"))
post_tuning_validation <- isTRUE(get0("post_tuning_validation", envir = .GlobalEnv, ifnotfound = FALSE))
exclude_regions <- get0("exclude_regions", envir = .GlobalEnv, ifnotfound = character(0))


# Block sizes for spatial diagnostics (run in both Step 2 and Step 4).
# If `cv_blocksize_scan` is set, use it; otherwise fall back to a single
# `cv_blocksize` only when `cv_type == "spatial"`.
cv_blocksize_scan <- get0("cv_blocksize_scan", envir = .GlobalEnv, ifnotfound = NULL)
use_spatial_diags <- !is.null(cv_blocksize_scan) &&
  length(cv_blocksize_scan) > 0L &&
  is.numeric(cv_blocksize_scan)
default_block_sizes_m <- if (isTRUE(use_spatial_diags)) {
  as.integer(cv_blocksize_scan)
} else if (identical(cv_type, "spatial")) {
  as.integer(cv_blocksize)
} else {
  integer(0)
}
buffer_km <- 1
buffer_m  <- buffer_km * 1000

# n_folds set from run_paper.R globals (get0 above); do not override
# n_folds <- 5

out_dir <- file.path(get0("cv_output_dir", envir = .GlobalEnv, ifnotfound = "output"), "cv_pipeline")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

hyperparams_by_model <- list()
if (post_tuning_validation) {
  config_dir <- out_dir
  load_best_config <- function(model_name) {
    base <- file.path(config_dir, "best_config")
    paths <- switch(model_name,
      XGB = c(paste0(base, "_xgb.rds"), file.path(config_dir, "xgb_best_config.rds")),
      GAM = c(paste0(base, "_gam.rds"), file.path(config_dir, "gam_best_config.rds")),
      GPR = c(paste0(base, "_gpr.rds"), file.path(config_dir, "gpr_best_config.rds")),
      NULL
    )
    if (is.null(paths)) return(NULL)
    for (p in paths) if (file.exists(p)) return(readRDS(p))
    NULL
  }
  for (m in model_list) {
    cfg <- load_best_config(m)
    if (is.null(cfg)) next
    if (m == "GPR") {
      hyperparams_by_model[[m]] <- list(
        kernel = cfg$kernel,
        nug.min = cfg$nug.min,
        nug.max = cfg$nug.max,
        nug.est = TRUE
      )
    } else if (m == "XGB") {
      hyperparams_by_model[[m]] <- list(
        nrounds = cfg$nrounds,
        max_depth = cfg$max_depth,
        learning_rate = cfg$learning_rate %||% 0.1,
        subsample = cfg$subsample %||% 0.8,
        colsample_bytree = cfg$colsample_bytree %||% 0.8,
        min_child_weight = cfg$min_child_weight %||% 1L
      )
    } else if (m == "GAM") {
      hyperparams_by_model[[m]] <- list(
        k_covariate = cfg$k_covariate %||% 6L
      )
    }
  }
  cat("Post-tuning validation: using best hyperparameters for models:",
      paste(names(hyperparams_by_model), collapse = ", "), "\n")
}

## ================================ DATA ================================
cat("\n========================================\n")
cat("LOADING CORE-LEVEL DATA\n")
cat("========================================\n\n")

dat <- read_rds("data/all_extracted_new.rds")
# Exclude necessary regions
if (length(exclude_regions) > 0L) {
  if (!"region" %in% names(dat)) dat <- assign_region_from_latlon(dat)
  dat <- dat[is.na(dat$region) | !dat$region %in% exclude_regions, ]
  cat("Excluded region(s):", paste(exclude_regions, collapse = ", "), "\n")
}
target_var <- "median_carbon_density_100cm"
# Predictors: raster-derived covariates that exist in dat (updated predictor set)
predictor_vars <- raster_covariates[raster_covariates %in% colnames(dat)]
if (length(predictor_vars) == 0) {
  stop("No predictor_vars found in dat. Check that all_extracted_new.rds contains columns matching raster config.")
}

# select columns in predictor_vars; ensure median_carbon_density for run_cv
extra_cols <- c()
if ("seagrass_species" %in% names(dat)) extra_cols <- c(extra_cols, "seagrass_species")
if ("region" %in% names(dat)) extra_cols <- c(extra_cols, "region")
complete_dat <- dat %>%
  dplyr::select(longitude, latitude, target_var, dplyr::all_of(predictor_vars), dplyr::all_of(extra_cols)) %>%
  dplyr::filter(complete.cases(.))
complete_dat$median_carbon_density <- complete_dat[[target_var]]
predictor_vars <- predictor_vars[predictor_vars %in% colnames(complete_dat)]
log_response <- isTRUE(get0("log_transform_target", envir = .GlobalEnv, ifnotfound = TRUE))

n_raw <- nrow(complete_dat)
n_unique_locs <- nrow(unique(complete_dat[, c("longitude", "latitude")]))
cat("Cores (complete cases):", n_raw, " | Unique locations:", n_unique_locs, "\n")

n_cores <- nrow(complete_dat)
cat("Working dataset:", n_cores, "rows,", length(predictor_vars), "predictors.\n")

# Load model-specific predictor sets: prefer SHAP when use_shap_per_model, else permutation
predictor_vars_by_model <- list()
use_shap_per_model <- isTRUE(get0("use_shap_per_model", envir = .GlobalEnv, ifnotfound = FALSE))
cov_dir   <- file.path(get0("cv_output_dir", envir = .GlobalEnv, ifnotfound = "output"), "covariate_selection")
shap_file <- file.path(cov_dir, "pruned_model_variables_shap.csv")
perm_file <- file.path(cov_dir, "pruned_model_variables_perm.csv")
src_file  <- if (use_shap_per_model && file.exists(shap_file)) shap_file else if (file.exists(perm_file)) perm_file else shap_file
if (file.exists(src_file)) {
  df <- read.csv(src_file, stringsAsFactors = FALSE)
  if (all(c("model", "variable") %in% names(df))) {
    for (m in unique(df$model)) {
      vars <- intersect(df$variable[df$model == m], colnames(complete_dat))
      if (length(vars) >= 2) predictor_vars_by_model[[m]] <- vars
    }
    if (length(predictor_vars_by_model) > 0)
      cat("Predictors (per model from ", basename(src_file), "): ", paste(sprintf("%s=%d", names(predictor_vars_by_model), lengths(predictor_vars_by_model)), collapse = ", "), "\n", sep = "")
  }
}

block_sizes_m <- default_block_sizes_m
autocor_note <- if (length(block_sizes_m) > 0)
  sprintf("Spatial block CV: %s km.", paste(block_sizes_m %/% 1000, collapse = ", ")) else "Random CV only (no spatial blocks)."
cat("CV type:", cv_type, " | Block sizes (m):", if (length(block_sizes_m)) paste(block_sizes_m, collapse = ", ") else "none", "\n")
cat(autocor_note, "\n\n")


## ================================ CV STRATEGIES ================================
cat("========================================\n")
cat("BUILDING CV STRATEGIES\n")
cat("========================================\n\n")

core_sf <- sf::st_as_sf(
  complete_dat,
  coords = c("longitude", "latitude"),
  crs = 4326,
  remove = FALSE
)
light_core_sf <- sf::st_as_sf(
  complete_dat %>% dplyr::select(-dplyr::all_of(predictor_vars)),
  coords = c("longitude", "latitude"),
  crs = 4326,
  remove = FALSE
) # create a version without all the variables to save memory when calculating spatial folds
cat("Core-level data:", nrow(complete_dat), "cores,", ncol(complete_dat), "variables.\n")

# CV strategies driven by cv_type (from run_paper.R)
# 1. True random: naive i.i.d. split (duplicate locations can leak between train/test)
true_random_folds <- sample(rep(seq_len(n_folds), length.out = n_cores))

# 2. Location-grouped random: all rows at the same (lon, lat) go to the same fold
loc_id <- as.integer(factor(paste(complete_dat$longitude, complete_dat$latitude)))
unique_loc_ids <- unique(loc_id)
loc_fold_assign <- sample(rep(seq_len(n_folds), length.out = length(unique_loc_ids)))
loc_grouped_folds <- loc_fold_assign[match(loc_id, unique_loc_ids)]

# 3. Pixel-grouped random: all rows sharing identical covariate vectors go to the
#    same fold (strictly more conservative than location-grouped, since distinct
#    locations that map to same raster pixel will also be grouped together)
pixel_info <- make_pixel_grouped_folds(complete_dat, predictor_vars, n_folds, seed = 42L)
pixel_grouped_folds <- pixel_info$fold_indices

cat("  True random folds:", n_cores, "rows ->", n_folds, "folds.\n")
cat("  Location-grouped random folds:", length(unique_loc_ids), "unique locations ->", n_folds, "folds.\n")
cat("  Pixel-grouped random folds:", pixel_info$n_groups, "unique covariate vectors ->", n_folds, "folds.\n")

cv_strategies <- list(
  list(method = "random_split", folds = true_random_folds, block_size_m = NA),
  list(method = "location_grouped_random", folds = loc_grouped_folds, block_size_m = NA),
  list(method = "pixel_grouped_random", folds = pixel_grouped_folds, block_size_m = NA)
)
if (length(block_sizes_m) > 0) {
  cat("\tCreating spatial folds (spatial block diagnostics)\n")
  for (size_m in block_sizes_m) {
    cat("  Spatial blocks ", size_m, " m ...\n", sep = "")
    cat("N.B. the progress bar is poor – it generally stays at 0, then leaps to 100\n", sep="")
    fold_info <- get_cached_spatial_folds(
      dat = complete_dat,
      block_size = size_m,
      n_folds = n_folds,
      cache_tag = "cv_pipeline_spatial",
      exclude_regions = exclude_regions,
      progress = TRUE
    )
    fi <- fold_info$fold_indices
    n_actual_folds <- if (is.vector(fi) && !is.list(fi)) max(fi, na.rm = TRUE) else length(fi)
    if (n_actual_folds < 2) {
      cat("  spatial block size", size_m, "m: too few folds (", n_actual_folds, "), skipping.\n")
    } else {
      cat("  Added spatial_block_", size_m, "m (", n_actual_folds, " folds).\n", sep = "")
      cv_strategies[[length(cv_strategies) + 1]] <- list(
        method = paste0("spatial_block_", size_m, "m"),
        folds = fi,
        block_size_m = size_m,
        blocks = NULL
      )
    }
  }

  # Region-stratified spatial CV: build spatial blocks WITHIN each region so every fold
  # has representation from all regions (avoids leave-one-region-out artifacts).
  if ("region" %in% names(complete_dat) && length(unique(complete_dat$region)) >= 2L) {
    strat_size_m <- cv_blocksize
    cat("\n  Region-stratified spatial blocks (", strat_size_m, " m) ...\n", sep = "")
    regions <- unique(complete_dat$region)
    strat_folds <- integer(n_cores)
    strat_ok <- TRUE
    for (reg in regions) {
      reg_idx <- which(complete_dat$region == reg)
      if (length(reg_idx) < n_folds) {
        strat_folds[reg_idx] <- sample(rep(seq_len(n_folds), length.out = length(reg_idx)))
        next
      }
      reg_fi <- tryCatch({
        fi_reg <- get_cached_spatial_folds(
          dat = complete_dat[reg_idx, , drop = FALSE],
          block_size = strat_size_m,
          n_folds = n_folds,
          cache_tag = paste0("cv_stratified_", gsub(" ", "_", tolower(reg))),
          exclude_regions = character(0),
          progress = FALSE
        )
        folds_to_vector(fi_reg$fold_indices, length(reg_idx))
      }, error = function(e) {
        sample(rep(seq_len(n_folds), length.out = length(reg_idx)))
      })
      strat_folds[reg_idx] <- reg_fi
    }
    if (all(strat_folds > 0L) && max(strat_folds) >= 2L) {
      cv_strategies[[length(cv_strategies) + 1]] <- list(
        method = paste0("region_stratified_", strat_size_m, "m"),
        folds = strat_folds,
        block_size_m = strat_size_m,
        blocks = NULL
      )
      cat("  Added region_stratified_", strat_size_m, "m (",
          max(strat_folds), " folds, balanced across ", length(regions), " regions).\n", sep = "")
      for (reg in regions) {
        reg_idx <- which(complete_dat$region == reg)
        cat("    ", reg, ": ", paste(table(strat_folds[reg_idx]), collapse = "/"), " per fold\n", sep = "")
      }
    }
  }
} else {
  cat("\tUsing random folds only (cv_type = ", cv_type, ")\n", sep = "")
}

if (isTRUE(quick_cv_check)) {
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
fold_vectors <- setNames(
  lapply(cv_strategies, function(s) folds_to_vector(s$folds, n_cores)),
  sapply(cv_strategies, function(s) s$method)
)
# Long-format data: method, fold, and all predictor values
fold_env_list <- list()
for (m in names(fold_vectors)) {
  fold_env_list[[m]] <- complete_dat %>%
    dplyr::mutate(method = m, fold = fold_vectors[[m]]) %>%
    dplyr::select(method, fold, dplyr::all_of(predictor_vars))
}
fold_env_df <- dplyr::bind_rows(fold_env_list)
fold_env_long <- fold_env_df %>%
  tidyr::pivot_longer(cols = dplyr::all_of(predictor_vars), names_to = "variable", values_to = "value") %>%
  mutate(method_label = vapply(method, method_to_display, character(1)))


## ================================ RUN CV FOR EACH STRATEGY ================================
cat("========================================\n")
cat("RUNNING CROSS-VALIDATION\n")
cat("========================================\n\n")

# Use model-specific pruned variables when available; else all predictors for all models
predictor_sets <- list()
if (length(predictor_vars_by_model) > 0) {
  predictor_sets$pruned_per_model <- predictor_vars_by_model
  cat("Using pruned variables per model from covariate selection.\n")
} else {
  predictor_sets$all <- predictor_vars
  cat("No per-model pruned sets found; using all ", length(predictor_vars), " predictors.\n", sep = "")
}
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
    use_per_model <- is.list(pvars) && !is.null(names(pvars)) && all(vapply(pvars, is.character, logical(1)))
    res <- tryCatch(
      {
        if (use_per_model) {
          run_cv(
            cv_method_name = s$method,
            fold_indices = s$folds,
            core_data = complete_dat,
            predictor_vars = predictor_vars,
            predictor_vars_by_model = pvars,
            tune_hyperparams = use_hyperparameter_tuning,
            nested_tuning = nested_tuning,
            verbose = FALSE,
            buffer_m = buf_m,
            core_sf = core_sf,
            return_predictions = TRUE,
            models = model_list,
            hyperparams_by_model = if (post_tuning_validation) hyperparams_by_model else NULL,
            log_response = log_response
          )
        } else {
          run_cv(
            cv_method_name = s$method,
            fold_indices = s$folds,
            core_data = complete_dat,
            predictor_vars = pvars,
            tune_hyperparams = use_hyperparameter_tuning,
            nested_tuning = nested_tuning,
            verbose = FALSE,
            buffer_m = buf_m,
            core_sf = core_sf,
            return_predictions = TRUE,
            models = model_list,
            hyperparams_by_model = if (post_tuning_validation) hyperparams_by_model else NULL,
            log_response = log_response
          )
        }
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

has_ps <- "predictor_set" %in% names(results_df)
filter_pruned <- function(d) {
  if (!"predictor_set" %in% names(d)) return(d)
  d[d[["predictor_set"]] %in% c("all", "pruned", "pruned_per_model"), , drop = FALSE]
}
grp <- c("method", "model", if (has_ps) "predictor_set")

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


suffix <- if (post_tuning_validation) "post_tuning" else "pre_tuning"
for (pair in list(
  list(results_df, sprintf("%s_cv_results_detailed.csv", suffix)),
  list(summary_by_method_model, sprintf("%s_cv_results_summary.csv", suffix))
)) safe_write_csv(pair[[1]], file.path(out_dir, pair[[2]]))


# Use only colours for models that appear in the data (keep order for legend)
models_present <- unique(results_df$model)
MODEL_COLOURS <- MODEL_COLOURS[names(MODEL_COLOURS) %in% models_present]

# Summary plot: R² and RMSE by method and model (x-axis: random left, then spatial blocks ascending)
metric_labeller <- c(r2 = "R²", rmse = "RMSE", mae = "MAE", bias = "Bias")
method_order <- order_methods_by_size(unique(results_df$method))
method_label_order <- vapply(method_order, method_to_display, character(1))
p_metrics <- results_df %>%
  dplyr::select(method, model, r2, rmse, mae, bias) %>%
  mutate(
    method_label = vapply(method, method_to_display, character(1)),
    method_label = factor(method_label, levels = method_label_order)
  ) %>%
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
ggsave(file.path(out_dir, sprintf("%s_cv_metrics_by_method_model.png", suffix)), p_metrics, width = 12, height = 8, dpi = dpi)
cat("\nSaved", file.path(out_dir, sprintf("%s_cv_metrics_by_method_model.png", suffix)), "\n")

# Scatter: observed vs predicted for all methods (one panel per model), colour by fold
if (!is.null(predictions_df) && nrow(predictions_df) > 0) {
  pred_all <- filter_pruned(predictions_df)
  if (nrow(pred_all) > 0) {
    pred_all <- pred_all %>% mutate(method_label = vapply(method, method_to_display, character(1)))
    p_scatter <- ggplot(pred_all, aes(x = predicted, y = observed, colour = as.factor(fold))) +
      geom_point(alpha = 0.6, size = 1) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "red", linewidth = 0.5) +
      facet_wrap(~model, scales = "free", ncol = 3) +
      scale_colour_viridis_d(name = "Fold", option = "plasma") +
      labs(x = "Predicted", y = "Observed", 
           title = if (show_titles) "Observed vs predicted by model (all CV folds)" else NULL) +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5),
            legend.position = "right")
    ggsave(file.path(out_dir, sprintf("%s_cv_observed_vs_predicted.png", suffix)), p_scatter, width = 12, height = 10, dpi = dpi)
    cat("Saved", file.path(out_dir, sprintf("%s_cv_observed_vs_predicted.png", suffix)), "\n")
  }
}

# Bar summary: mean R² and mean RMSE by method (averaged across models)
grp_fold   <- c("method", "fold", if (has_ps) "predictor_set")
grp_method <- c("method", if (has_ps) "predictor_set")
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

# To control bar order such that non-spatial methods appear first,
# followed by spatial methods ordered from smallest to largest block size,
# define an explicit method order.

# Function to extract a numeric km value from method strings (eg "spatial_block_1km" => 1)
get_block_km <- function(x) {
  m <- regmatches(x, regexec("spatial_block_(\\d+)km", x))
  as.numeric(vapply(m, function(val) if (length(val) > 1) val[2] else NA_real_, numeric(1)))
}

# Find all method types
methods_all <- unique(summary_method$method)
spatial_methods <- grep("^spatial_block_\\d+km$", methods_all, value = TRUE)
non_spatial_methods <- setdiff(methods_all, spatial_methods)

# Order spatial methods by their km value
spatial_methods_ordered <- spatial_methods[order(get_block_km(spatial_methods))]
method_order <- c(non_spatial_methods, spatial_methods_ordered)

# Ensure this order is reflected in method_label via a factor
summary_method$method_label <- factor(
  vapply(summary_method$method, method_to_display, character(1)),
  levels = vapply(method_order, method_to_display, character(1))
)

p_method <- summary_method %>%
  tidyr::pivot_longer(cols = c(r2, rmse), names_to = "metric", values_to = "value") %>%
  ggplot(aes(x = method_label, y = value, fill = metric)) +
  geom_col(position = "dodge") +
  labs(x = "CV method", y = "Mean (across folds and models)", title = if (show_titles) "Ensemble performance by CV method" else NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
if (has_ps) p_method <- p_method + facet_wrap(~predictor_set, ncol = 1)
ggsave(file.path(out_dir, sprintf("%s_cv_ensemble_by_method.png", suffix)), p_method, width = 10, height = 5, dpi = dpi)

# ## ================================ MODEL SELECTION FOR GAP-FILLING ================================
# # Goal: find a model with good performance at predicting carbon density at points *near* training data,
# # validated by train-test splits. Most relevant CV strategies:
# #   - random_split: random 80/20 (generalisation; test points are "near" in the sense of the full sample)
# #   - distance_buffer_1km: test points >1 km from train (strict "fill gaps between samples" test)
# # Rank models by mean R² (and mean RMSE) on these strategies; prefer distance-buffered for gap-filling.
# gap_cv_methods <- c("random_split", paste0("distance_buffer_", buffer_km, "km"))
# summary_gap <- summary_by_method_model %>%
#   filter(method %in% gap_cv_methods) %>%
#   filter_pruned()
# if (nrow(summary_gap) > 0) {
#   # Rank by mean R² (descending) then by mean RMSE (ascending) within each CV method
#   summary_gap <- summary_gap %>%
#     group_by(method) %>%
#     mutate(
#       rank_r2 = row_number(desc(mean_r2)),
#       rank_rmse = row_number(mean_rmse)
#     ) %>%
#     ungroup()
#   # Overall rank: average rank across the two gap-filling CV methods (lower = better)
#   rank_overall <- summary_gap %>%
#     group_by(model) %>%
#     summarise(
#       mean_r2_avg = mean(mean_r2, na.rm = TRUE),
#       mean_rmse_avg = mean(mean_rmse, na.rm = TRUE),
#       rank_r2_avg = mean(rank_r2, na.rm = TRUE),
#       rank_rmse_avg = mean(rank_rmse, na.rm = TRUE),
#       .groups = "drop"
#     ) %>%
#     mutate(rank_overall = (rank_r2_avg + rank_rmse_avg) / 2) %>%
#     arrange(rank_overall)
#   write.csv(rank_overall, file.path(out_dir, "cv_best_models_gap_filling.csv"), row.names = FALSE)
#   write.csv(summary_gap, file.path(out_dir, "cv_gap_filling_by_method_model.csv"), row.names = FALSE)
#   cat("\n========================================\n")
#   cat("MODEL SELECTION FOR GAP-FILLING (predict near training data)\n")
#   cat("========================================\n\n")
#   cat("CV strategies used: ", paste(gap_cv_methods, collapse = ", "), "\n", sep = "")
#   cat("(random_split = generalisation; distance_buffer = test points >1 km from train)\n\n")
#   cat("Model ranking (by average rank of R² and RMSE across these CV methods):\n")
#   print(rank_overall %>% select(model, mean_r2_avg, mean_rmse_avg, rank_overall))
#   cat("\nTop 3 models for gap-filling:\n")
#   top3 <- head(rank_overall, 3)
#   for (i in seq_len(nrow(top3))) {
#     cat("  ", i, ". ", top3$model[i],
#       " (mean R² = ", round(top3$mean_r2_avg[i], 4),
#       ", mean RMSE = ", round(top3$mean_rmse_avg[i], 6), ")\n",
#       sep = ""
#     )
#   }
#   cat("\nSaved: cv_best_models_gap_filling.csv, cv_gap_filling_by_method_model.csv\n")
# } else {
#   cat("\nNo results for gap-filling CV methods; skipping model ranking.\n")
# }

cat("\nTesting performance over block sizes\n")

# 2) Performance over block sizes
spatial_results <- results_df %>%
  filter(grepl("^spatial_block_", method), !is.na(block_size_m)) %>%
  filter_pruned()
if (nrow(spatial_results) > 0 && "block_size_m" %in% names(spatial_results)) {
  spatial_results <- spatial_results %>%
    mutate(block_size_km = round(block_size_m / 1000))
  sum_block <- spatial_results %>%
    group_by(block_size_km, model) %>%
    summarise(mean_r2 = mean(r2, na.rm = TRUE), mean_rmse = mean(rmse, na.rm = TRUE), .groups = "drop")
  km_levels <- sort(unique(sum_block$block_size_km))
  sum_block <- sum_block %>%
    mutate(block_km = factor(block_size_km, levels = km_levels))
  n_blocks <- length(km_levels)
  p_block_r2 <- ggplot(sum_block, aes(x = block_km, y = mean_r2, colour = model, group = model)) +
    geom_point(size = 3) +
    scale_colour_manual(values = MODEL_COLOURS, na.value = "gray75") +
    labs(x = "Block size (km)", y = "Mean R²") +
    scale_x_discrete(labels = paste0(km_levels, " km")) +
    theme_minimal() +
    theme(legend.position = "bottom", plot.title = element_text(size = 10, hjust = 0.5))
  if (n_blocks > 1) p_block_r2 <- p_block_r2 + geom_line(linewidth = 1)
  p_block_rmse <- ggplot(sum_block, aes(x = block_km, y = mean_rmse, colour = model, group = model)) +
    geom_point(size = 3) +
    scale_colour_manual(values = MODEL_COLOURS, na.value = "gray75") +
    labs(x = "Block size (km)", y = "Mean RMSE") +
    scale_x_discrete(labels = paste0(km_levels, " km")) +
    theme_minimal() +
    theme(legend.position = "bottom", plot.title = element_text(size = 10, hjust = 0.5))
  if (n_blocks > 1) p_block_rmse <- p_block_rmse + geom_line(linewidth = 1)
  p_block_combined <- p_block_r2 + p_block_rmse + plot_layout(ncol = 1) + plot_annotation(title = if (show_titles) "Performance over spatial block sizes" else NULL)
  ggsave(file.path(out_dir, sprintf("%s_cv_performance_over_block_sizes.png", suffix)), p_block_combined, width = 9, height = 7, dpi = dpi)
  cat("Saved", file.path(out_dir, sprintf("%s_cv_performance_over_block_sizes.png", suffix)), "\n")
}


# 3) Illustration: spatial blocks at 25%, 50%, 75%, 100% superposed on scatter of data points (one fold)
# Build list of (block_size_m, blocks_sf, folds) for each spatial block size; use stored blocks or re-run cv_spatial
block_illustration_data <- list()
for (idx in seq_along(cv_strategies)) {
  s <- cv_strategies[[idx]]
  if (is.na(s$block_size_m)) next
  blocks_sf <- s$blocks
  if (is.null(blocks_sf)) {
    cv_opts <- list(x = core_sf, k = n_folds, size = s$block_size_m, selection = "random",
      iteration = 50, progress = TRUE, biomod2 = FALSE, hexagon = FALSE, plot = FALSE)
    cv_out <- tryCatch(do.call(blockCV::cv_spatial, cv_opts), error = function(e) NULL)
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
    fold_vec <- folds_to_vector(bl$fold_list, n_cores)
    points_with_fold <- core_sf
    points_with_fold$fold <- fold_vec
    points_with_fold$is_test <- points_with_fold$fold == highlight_fold
    blocks_sf <- bl$blocks_sf
    fold_col <- if ("folds" %in% names(blocks_sf)) "folds" else if ("fold" %in% names(blocks_sf)) "fold" else setdiff(names(blocks_sf), "geometry")[1]
    if (is.na(fold_col) || length(fold_col) == 0) fold_col <- "folds"
    block_km <- round(bl$block_size_m / 1000)
    p_block <- ggplot() +
      geom_sf(data = blocks_sf, aes(fill = as.factor(.data[[fold_col]])), alpha = 0.2, colour = "gray50", linewidth = 0.4) +
      geom_sf(data = points_with_fold, aes(color = is_test), size = 1.5, alpha = 0.95, stroke = 0) +
      scale_fill_viridis_d(option = "E", name = "Fold", na.value = NA) +
      scale_color_manual(values = c("FALSE" = "gray20", "TRUE" = "red"), labels = c("Train", "Test (fold 1)"), name = NULL) +
      labs(title = if (show_titles) paste0("Spatial block ", block_km, " km") else NULL, x = "Longitude", y = "Latitude") +
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
  ggsave(file.path(out_dir, sprintf("%s_cv_spatial_blocks_and_points.png", suffix)), p_blocks_illustration,
    width = 10, height = max(5, 2.5 * ceiling(n_panels / 2)), dpi = dpi
  )
  cat("Saved", file.path(out_dir, sprintf("%s_cv_spatial_blocks_and_points.png", suffix)), "\n")
}

## ================================ NOTES ================================
cat("\nCV pipeline complete.\n")
cat(autocor_note, "\n")
cat("Block sizes (m):", paste(block_sizes_m, collapse = ", "), "\n")
cat("Plots saved to:", out_dir, "\n")
