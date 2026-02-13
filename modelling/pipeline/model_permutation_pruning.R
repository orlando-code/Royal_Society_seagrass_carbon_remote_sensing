# Model-wise permutation-importance pruning at 1 km spatial CV
# 
# For each model used in cv_pipeline, this script:
# - Uses the 1 km spatial block folds (spatial_block_1000m)
# - Runs a baseline CV with the full predictor set
# - Computes permutation importance for each covariate (RMSE increase when permuted)
# - Keeps variables up to a cumulative fraction of total importance (default 0.95)
# - Writes a single CSV with results for all models:
#   modelling/cv_pipeline_output/cv_model_permutation_pruning_1km.csv
# 
# This is an "investigation" script: it does not change cv_pipeline behaviour directly,
# but cv_pipeline_plots.R and other scripts can use the cached CSV.
# Cache is config-aware: if the CSV and a matching _config.rds exist, the script
# loads from cache and skips the pipeline. Change any config (e.g. permutation_keep_frac,
# cor_threshold) or delete the CSV/config to force recompute.

rm(list = ls())
setwd(here::here())
set.seed(42)

source("modelling/R/helpers.R")
source("modelling/R/extract_covariates_from_rasters.R")
load_packages(c("here", "tidyverse", "ggplot2", "blockCV", "sf"))

target_var <- "median_carbon_density_100cm"
out_dir <- "figures/cv_pipeline_output"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Config: overridden by run_paper_figures.R via paper_pipeline_config.rds when present
paper_config <- NULL
if (file.exists(file.path(out_dir, "paper_pipeline_config.rds"))) {
  paper_config <- tryCatch(readRDS(file.path(out_dir, "paper_pipeline_config.rds")), error = function(e) NULL)
}
# - use_correlation_filter / correlation_filter_on: drop one of each highly correlated predictor pair before permutation
# - cor_threshold / correlation_filter_threshold: absolute correlation above which a pair is "high"
# - permutation_keep_frac / model_permutation_pruning_keep_frac: keep variables up to this cumulative fraction of importance
# - run_permutation / model_permutation_pruning_on: if FALSE, only apply correlation filter and output that set (no permutation importance)
use_correlation_filter <- if (!is.null(paper_config) && "correlation_filter_on" %in% names(paper_config)) paper_config$correlation_filter_on else TRUE
cor_threshold <- if (!is.null(paper_config) && "correlation_filter_threshold" %in% names(paper_config)) paper_config$correlation_filter_threshold else 0.8
permutation_keep_frac <- if (!is.null(paper_config) && "model_permutation_pruning_keep_frac" %in% names(paper_config)) paper_config$model_permutation_pruning_keep_frac else 0.95
run_permutation <- if (!is.null(paper_config) && "model_permutation_pruning_on" %in% names(paper_config)) paper_config$model_permutation_pruning_on else TRUE
n_permutations <- 1        # per variable per model (can be increased if needed)

block_size_m <- 1000
block_size_name <- paste0("spatial_block_", block_size_m, "m")

cache_csv <- file.path(out_dir, paste0("cv_model_permutation_pruning_", block_size_name, ".csv"))
cache_config_path <- sub("\\.csv$", "_config.rds", cache_csv)
current_config <- list(
  permutation_keep_frac = permutation_keep_frac,
  n_permutations = n_permutations,
  use_correlation_filter = use_correlation_filter,
  cor_threshold = cor_threshold,
  block_size_m = block_size_m,
  run_permutation = run_permutation
)
loaded_from_cache <- FALSE
# Only use cache when we did a full permutation run (saved and current both have run_permutation = TRUE)
if (run_permutation && file.exists(cache_csv) && file.exists(cache_config_path)) {
  saved_config <- tryCatch(readRDS(cache_config_path), error = function(e) NULL)
  if (!is.null(saved_config) && identical(current_config, saved_config)) {
    importance_all <- read.csv(cache_csv, stringsAsFactors = FALSE)
    loaded_from_cache <- TRUE
    cat("Loaded results from cache (config unchanged).\n")
    cat("  ", cache_csv, "\n\n", sep = "")
  }
}

if (!loaded_from_cache) {

cat("\n========================================\n")
cat("MODEL-WISE PERMUTATION PRUNING AT", block_size_m, "m\n")
cat("========================================\n\n")

# ---------------------------------------------------------------------------
# 1. Load core-level data (same as cv_pipeline.R)
# ---------------------------------------------------------------------------
dat <- read_rds("data/all_extracted_new.rds")
RASTER_DIR <- "data/env_rasters"
raster_covariates <- names(build_covariate_config_from_dir(RASTER_DIR))
predictor_vars <- raster_covariates[raster_covariates %in% colnames(dat)]
if (length(predictor_vars) == 0) {
  stop("No predictor_vars found in dat. Check that all_extracted_new.rds contains columns matching raster config.")
}

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
cat("Predictors (all):", length(predictor_vars), "\n\n")

# Optional: correlation-based filtering to remove redundant predictors before permutation
if (use_correlation_filter && length(predictor_vars) >= 2) {
  cat("Applying correlation-based predictor filter with threshold |r| >", cor_threshold, "...\n")
  # core_data_complete stores the response as 'median_carbon_density'
  dat_complete <- core_data_complete[, c("median_carbon_density", predictor_vars), drop = FALSE]
  dat_complete <- dat_complete[complete.cases(dat_complete), ]
  if (nrow(dat_complete) >= 3) {
    # Correlations with target (for deciding which variable in a pair to keep)
    cor_with_target <- sapply(predictor_vars, function(v) {
      if (v %in% names(dat_complete)) {
        stats::cor(dat_complete[[v]], dat_complete[["median_carbon_density"]], use = "complete.obs")
      } else {
        NA_real_
      }
    })
    cor_df <- data.frame(
      variable = predictor_vars,
      correlation = cor_with_target,
      abs_correlation = abs(cor_with_target),
      stringsAsFactors = FALSE
    )
    cor_df <- cor_df[!is.na(cor_df$correlation), ]
    vars_for_cor <- cor_df$variable

    if (length(vars_for_cor) >= 2) {
      covar_data <- dat_complete[, vars_for_cor, drop = FALSE]
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
        # Work only on upper triangle to avoid duplicates, and process strongest pairs first
        pair_df <- pair_df[pair_df$var1 < pair_df$var2, ]
        pair_df <- pair_df[order(abs(pair_df$correlation), decreasing = TRUE), ]

        for (i in seq_len(nrow(pair_df))) {
          v1 <- pair_df$var1[i]
          v2 <- pair_df$var2[i]
          if (v1 %in% vars_removed || v2 %in% vars_removed) next
          c1 <- cor_df$abs_correlation[cor_df$variable == v1]
          c2 <- cor_df$abs_correlation[cor_df$variable == v2]
          if (length(c1) == 0 || length(c2) == 0) next
          # If correlations are identical to 3 dp, skip automatic removal (keep both)
          if (round(c1, 3) == round(c2, 3)) next
          if (c1 < c2) {
            vars_removed <- c(vars_removed, v1)
          } else {
            vars_removed <- c(vars_removed, v2)
          }
        }
      }

      vars_removed <- unique(vars_removed)
      if (length(vars_removed) > 0) {
        cat("  Correlation filter removed", length(vars_removed), "predictors:\n")
        cat("   ", paste(vars_removed, collapse = ", "), "\n")
      } else {
        cat("  Correlation filter did not remove any predictors.\n")
      }
      predictor_vars <- setdiff(predictor_vars, vars_removed)
    } else {
      cat("  Not enough predictors with valid correlations; skipping correlation filter.\n")
    }
  } else {
    cat("  Not enough complete cases for correlation filter; skipping.\n")
  }
  cat("Predictors after correlation filter:", length(predictor_vars), "\n\n")
}

  if (run_permutation) {
# ---------------------------------------------------------------------------
# 2. Load N km spatial folds (from cv_pipeline cache)
# ---------------------------------------------------------------------------
spatial_folds_cache_path <- file.path(out_dir, "spatial_folds_cache.rds")
if (!file.exists(spatial_folds_cache_path)) {
  stop("spatial_folds_cache.rds not found in ", out_dir, ". Run cv_pipeline.R first.")
}

cached <- tryCatch(readRDS(spatial_folds_cache_path), error = function(e) NULL)
if (is.null(cached) || !"spatial_strategies" %in% names(cached)) {
  stop("spatial_folds_cache.rds does not contain spatial_strategies; re-run cv_pipeline.R.")
}

strategy_block_size <- NULL
for (s in cached$spatial_strategies) {
  if (s$method == block_size_name) {
    strategy_block_size <- s
    break
  }
}
if (is.null(strategy_block_size)) {
  stop(paste0("spatial strategy (", block_size_name, ") not found in spatial_folds_cache.rds."))
}

folds_ids <- strategy_block_size$folds
if (length(folds_ids) != nrow(core_data_complete)) {
  stop(paste0(block_size_name, " folds length (", length(folds_ids), ") != nrow(core_data_complete) (", nrow(core_data_complete), ")."))
}
cat("Using", length(unique(folds_ids)), "folds for ", block_size_name, " spatial CV.\n\n")

core_sf <- sf::st_as_sf(
  core_data_complete,
  coords = c("longitude", "latitude"),
  crs = 4326,
  remove = FALSE
)

# ---------------------------------------------------------------------------
# 3. Helper: run CV once and return mean RMSE per model
# ---------------------------------------------------------------------------
run_cv_rmse_by_model <- function(core_data_use, predictor_vars_use) {
  res <- run_cv(
    cv_method_name = block_size_name,
    fold_indices = folds_ids,
    core_data = core_data_use,
    predictor_vars = predictor_vars_use,
    tune_hyperparams = FALSE,
    nested_tuning = TRUE,
    verbose = FALSE,
    buffer_m = NA_real_,
    core_sf = core_sf,
    return_predictions = FALSE,
    skip_inla = TRUE
  )
  if (!("model" %in% names(res) && "rmse" %in% names(res))) {
    stop("run_cv did not return expected columns 'model' and 'rmse'.")
  }
  res %>%
    dplyr::group_by(model) %>%
    dplyr::summarise(baseline_rmse = mean(rmse, na.rm = TRUE), .groups = "drop")
}

# ---------------------------------------------------------------------------
# 4. Baseline RMSE for all models with full predictor set
# ---------------------------------------------------------------------------
cat("Computing baseline ", block_size_name, " spatial CV RMSE for all models...\n")
baseline_rmse_by_model <- run_cv_rmse_by_model(core_data_complete, predictor_vars)
print(baseline_rmse_by_model)
cat("\n")

# ---------------------------------------------------------------------------
# 5. Permutation importance by model and variable
# ---------------------------------------------------------------------------
all_importance <- list()
n_vars <- length(predictor_vars)
cat("Starting permutation importance for", n_vars, "predictors.\n\n")
for (i in seq_along(predictor_vars)) {
  v <- predictor_vars[i]
  remaining <- n_vars - i
  cat(
    "Permutation importance for variable ", v,
    " (", i, " of ", n_vars, "; remaining ", remaining, ")\n",
    sep = ""
  )
  rmse_inc_list <- list()
  
  for (perm in seq_len(n_permutations)) {
    core_perm <- core_data_complete
    core_perm[[v]] <- sample(core_perm[[v]])
    
    rmse_perm <- run_cv_rmse_by_model(core_perm, predictor_vars)
    
    joined <- dplyr::left_join(
      rmse_perm,
      baseline_rmse_by_model,
      by = "model",
      suffix = c("_perm", "_base")
    )
    joined$rmse_increase <- joined$baseline_rmse_perm - joined$baseline_rmse_base
    rmse_inc_list[[length(rmse_inc_list) + 1]] <- joined[, c("model", "rmse_increase")]
  }
  
  rmse_inc_df <- dplyr::bind_rows(rmse_inc_list)
  rmse_inc_summary <- rmse_inc_df %>%
    dplyr::group_by(model) %>%
    dplyr::summarise(rmse_increase = mean(rmse_increase, na.rm = TRUE), .groups = "drop")
  rmse_inc_summary$variable <- v
  all_importance[[length(all_importance) + 1]] <- rmse_inc_summary
}

importance_all <- dplyr::bind_rows(all_importance)
importance_all$rmse_increase[!is.finite(importance_all$rmse_increase)] <- 0

# ---------------------------------------------------------------------------
# 6. For each model, sort by importance, compute cumulative fraction, choose kept vars
# ---------------------------------------------------------------------------
importance_all <- importance_all %>%
  dplyr::group_by(model) %>%
  dplyr::arrange(dplyr::desc(rmse_increase), .by_group = TRUE) %>%
  dplyr::mutate(
    cum_rmse = cumsum(pmax(rmse_increase, 0)),
    total_rmse = max(cum_rmse, na.rm = TRUE),
    cum_frac = dplyr::if_else(total_rmse > 0, cum_rmse / total_rmse, 1),
    keep = cum_frac <= permutation_keep_frac
  ) %>%
  dplyr::ungroup()

importance_all$permutation_keep_frac <- permutation_keep_frac

  saveRDS(current_config, cache_config_path)
  write.csv(importance_all, cache_csv, row.names = FALSE)
  cat("Saved model-wise permutation pruning results and config to:\n  ", cache_csv, "\n  ", cache_config_path, "\n\n", sep = "")
  } else {
    # Correlation-only: no permutation; output same CSV format so run_paper_figures.R can copy to pruned CSVs
    model_names <- c("Random Forest", "XGBoost", "SVM", "Neural Network", "GPR", "Ordinary Kriging", "Universal Kriging (env drift)", "GAM", "INLA")
    importance_all <- expand.grid(model = model_names, variable = predictor_vars, stringsAsFactors = FALSE)
    importance_all$keep <- TRUE
    importance_all$rmse_increase <- 0
    importance_all$cum_rmse <- 0
    importance_all$total_rmse <- 0
    importance_all$cum_frac <- 1
    importance_all$permutation_keep_frac <- NA_real_
    write.csv(importance_all, cache_csv, row.names = FALSE)
    cat("Saved correlation-filtered variable list (permutation skipped).\n  ", cache_csv, "\n\n", sep = "")
  }
}

cat("Permutation-importance pruning summary (head):\n")
print(head(importance_all))
cat("\n")
if (loaded_from_cache) {
  cat("(Cache used; delete ", cache_csv, " and ", cache_config_path, " to force recompute.)\n\n", sep = "")
}
cat("Done.\n")

