## Robust hyperparameter tuning for `pixel_grouped_random`.
##
## Hyperparameter selection objective:
##   - For each candidate hyperparameter set, evaluate CV performance
##     on *multiple* pixel-grouped fold instantiations (`robust_fold_seed_list`).
##   - Select the candidate with the minimum mean RMSE averaged across seeds.
##
## This script expects robust per-model predictor sets to exist, typically from:
##   robust_shap_covariate_pruning.R
##
## Outputs:
##   output/<cv_regime_name>/cv_pipeline/robust_pixel_grouped_random_tuning/
##     - best_config_<model>_robust.rds
##     - robust_cv_metrics_<model>.csv
##
project_root <- here::here()
setwd(project_root)

source(file.path(project_root, "modelling/R/helpers.R"))
source(file.path(project_root, "modelling/R/assign_region_from_latlon.R"))
load_packages(c("here", "dplyr", "readr", "mgcv", "randomForest", "GauPro", "xgboost", "sf"))

cv_regime_name <- get0("cv_regime_name", envir = .GlobalEnv, ifnotfound = "pixel_grouped")
cv_type <- get0("cv_type", envir = .GlobalEnv, ifnotfound = "pixel_grouped_random")
stopifnot(identical(cv_type, "pixel_grouped_random"))
cv_type_hash <- "pixel_grouped"

target_var <- get0("target_var", envir = .GlobalEnv, ifnotfound = "median_carbon_density_100cm")
log_transform_target <- isTRUE(get0("log_transform_target", envir = .GlobalEnv, ifnotfound = TRUE))

exclude_regions <- get0("exclude_regions", envir = .GlobalEnv, ifnotfound = character(0))
n_folds <- as.integer(get0("n_folds", envir = .GlobalEnv, ifnotfound = 5L))
cv_blocksize <- get0("cv_blocksize", envir = .GlobalEnv, ifnotfound = 5000L) # unused for pixel_grouped_random

robust_fold_seed_list <- get0("robust_fold_seed_list", envir = .GlobalEnv, ifnotfound = c(42L, 43L))
robust_fold_seed_list <- as.integer(robust_fold_seed_list)

robust_dir <- file.path(project_root, "output", cv_regime_name, "cv_pipeline", "robust_pixel_grouped_random_tuning")
dir.create(robust_dir, recursive = TRUE, showWarnings = FALSE)

cat("Robust hyperparameter tuning (pixel_grouped_random)\n")
cat("  robust_fold_seed_list:", paste(robust_fold_seed_list, collapse = ", "), "\n")
cat("  n_folds:", n_folds, "\n")
cat("  writing outputs to:", robust_dir, "\n")

# ---------------------------------------------------------------------------
# Load data and build consistent complete-case frame for folds
# ---------------------------------------------------------------------------
dat <- readr::read_rds(file.path(project_root, "data/all_extracted_new.rds"))
if (length(exclude_regions) > 0L) {
  if (!"region" %in% names(dat)) dat <- assign_region_from_latlon(dat)
  dat <- dat[is.na(dat$region) | !dat$region %in% exclude_regions, ]
}

predictor_vars_full <- setdiff(
  colnames(dat),
  c("latitude", "longitude", "number_id_final_version", "seagrass_species", "region", target_var)
)
predictor_vars_full <- predictor_vars_full[predictor_vars_full %in% colnames(dat)]

species_region_cols <- intersect(c("seagrass_species", "region"), names(dat))

complete_dat <- dat %>%
  dplyr::select(
    longitude,
    latitude,
    dplyr::all_of(target_var),
    dplyr::all_of(predictor_vars_full),
    dplyr::all_of(species_region_cols)
  ) %>%
  dplyr::filter(complete.cases(.))

complete_dat$median_carbon_density <- complete_dat[[target_var]]
predictor_vars_full <- predictor_vars_full[predictor_vars_full %in% colnames(complete_dat)]

stopifnot("median_carbon_density" %in% names(complete_dat))
core_data <- as.data.frame(complete_dat)

cat("  Complete-case rows:", nrow(core_data), "\n")
cat("  Predictor vars for fold hashing:", length(predictor_vars_full), "\n")

# ---------------------------------------------------------------------------
# Load robust pruned predictor sets
# ---------------------------------------------------------------------------
seeds_str <- paste(robust_fold_seed_list, collapse = "-")
robust_pruned_importance_type <- get0("robust_pruned_importance_type", envir = .GlobalEnv, ifnotfound = "perm")
robust_pruned_importance_type <- match.arg(robust_pruned_importance_type, choices = c("perm", "shap"))
robust_pruned_csv_override <- get0("robust_pruned_csv_override", envir = .GlobalEnv, ifnotfound = NA_character_)

robust_cov_dir <- file.path(project_root, "output", cv_regime_name, "covariate_selection", "robust_pixel_grouped_random")
default_robust_pruned_csv <- if (identical(robust_pruned_importance_type, "shap")) {
  file.path(
    robust_cov_dir,
    paste0("pruned_model_variables_shap_robust_pixel_grouped_random_seeds_", seeds_str, ".csv")
  )
} else {
  file.path(
    robust_cov_dir,
    paste0("pruned_model_variables_perm_robust_pixel_grouped_random_seeds_", seeds_str, ".csv")
  )
}

robust_pruned_csv <- if (!is.na(robust_pruned_csv_override) && nzchar(robust_pruned_csv_override)) {
  robust_pruned_csv_override
} else {
  default_robust_pruned_csv
}

if (!file.exists(robust_pruned_csv)) {
  fallback_csv <- file.path(
    project_root, "output", cv_regime_name, "covariate_selection",
    if (identical(robust_pruned_importance_type, "shap")) "pruned_model_variables_shap.csv" else "pruned_model_variables_perm.csv"
  )
  if (!file.exists(fallback_csv)) stop("Missing fallback pruned vars CSV: ", fallback_csv)
  cat("\nWARNING: robust pruned vars not found:\n  ", robust_pruned_csv,
      "\nFalling back to non-robust pruned vars:\n  ", fallback_csv, "\n")
  robust_pruned_csv <- fallback_csv
}

pruned_df <- read.csv(robust_pruned_csv, stringsAsFactors = FALSE)
stopifnot(all(c("model", "variable") %in% names(pruned_df)))

predictor_vars_by_model <- list()
for (m in unique(pruned_df$model)) {
  vars <- intersect(pruned_df$variable[pruned_df$model == m], colnames(core_data))
  if (length(vars) >= 2) predictor_vars_by_model[[m]] <- vars
}
models <- intersect(names(predictor_vars_by_model), c("GPR", "GAM", "XGB"))
if (length(models) == 0L) stop("No models found after loading robust pruned variables.")
cat("  Robust predictor sets loaded for:", paste(models, collapse = ", "), "\n")

# ---------------------------------------------------------------------------
# Baseline configs (only used to set candidate ranges, not for objective)
# ---------------------------------------------------------------------------
config_dir <- file.path(project_root, "output", cv_regime_name, "cv_pipeline")
cfg_xgb_path <- file.path(config_dir, "best_config_xgb.rds")
cfg_gam_path <- file.path(config_dir, "best_config_gam.rds")
cfg_gpr_path <- file.path(config_dir, "best_config_gpr.rds")

cfg_xgb <- if (file.exists(cfg_xgb_path)) readRDS(cfg_xgb_path) else NULL
cfg_gam <- if (file.exists(cfg_gam_path)) readRDS(cfg_gam_path) else NULL
cfg_gpr <- if (file.exists(cfg_gpr_path)) readRDS(cfg_gpr_path) else NULL

# ---------------------------------------------------------------------------
# Fold builder per seed
# ---------------------------------------------------------------------------
get_fold_indices_for_seed <- function(seed) {
  pf <- make_cv_folds(
    dat = core_data,
    covariate_cols = predictor_vars_full,
    n_folds = n_folds,
    cv_type = cv_type_hash,
    cv_blocksize = cv_blocksize,
    exclude_regions = exclude_regions,
    cache_tag = paste0("robust_hp_", cv_type, "_seed_", seed),
    seed = seed
  )
  pf$fold_indices
}

# ---------------------------------------------------------------------------
# Robust objective helper: mean RMSE across seeds and folds
# ---------------------------------------------------------------------------
robust_mean_rmse <- function(model_name, hp_by_model) {
  rmses <- c()
  for (seed in robust_fold_seed_list) {
    fold_indices <- get_fold_indices_for_seed(seed)
    res <- run_cv(
      cv_method_name = paste0("robust_hp_", model_name, "_seed_", seed),
      fold_indices = fold_indices,
      core_data = core_data,
      predictor_vars = predictor_vars_full,
      predictor_vars_by_model = predictor_vars_by_model,
      hyperparams_by_model = setNames(list(hp_by_model), model_name),
      tune_hyperparams = FALSE,
      nested_tuning = FALSE,
      verbose = FALSE,
      return_predictions = FALSE,
      models = c(model_name),
      log_response = log_transform_target
    )
    if (!is.null(res) && nrow(res) > 0L) {
      rmses <- c(rmses, res$rmse)
    }
  }
  mean(rmses, na.rm = TRUE)
}

# ---------------------------------------------------------------------------
# Tune XGB robustly (candidate regularization only)
# ---------------------------------------------------------------------------
if ("XGB" %in% models) {
  cat("\n=== Robust tuning: XGB ===\n")
  if (is.null(cfg_xgb)) stop("Missing best_config_xgb.rds for XGB baseline.")

  nrounds_grid          <- c(50L, 100L)
  max_depth_grid        <- c(1L, 2L, 3L, 4L, 6L)
  learning_rate_grid    <- c(0.01, 0.05, 0.1)
  min_child_weight_grid <- c(1L, 5L, 10L, 20L)
  subsample_grid        <- c(0.7, 1.0)
  colsample_bytree_grid <- c(0.6, 1.0)
  min_split_loss_grid   <- c(0, 1, 5, 10)
  reg_reg_lambda_grid   <- c(1, 10, 50, 100)

  full_grid <- expand.grid(
    nrounds          = nrounds_grid,
    max_depth        = max_depth_grid,
    learning_rate    = learning_rate_grid,
    min_child_weight = min_child_weight_grid,
    subsample        = subsample_grid,
    colsample_bytree = colsample_bytree_grid,
    min_split_loss   = min_split_loss_grid,
    reg_reg_lambda   = reg_reg_lambda_grid,
    stringsAsFactors = FALSE
  )
  cat("  Full grid size:", nrow(full_grid), "— sampling up to 120 candidates\n")
  set.seed(42)
  max_candidates <- as.integer(get0("xgb_max_grid_candidates", envir = .GlobalEnv, ifnotfound = 120L))
  if (nrow(full_grid) > max_candidates) {
    sampled_idx <- sample.int(nrow(full_grid), max_candidates)
    baseline_row <- data.frame(
      nrounds = cfg_xgb$nrounds %||% 100L,
      max_depth = cfg_xgb$max_depth %||% 6L,
      learning_rate = cfg_xgb$learning_rate %||% 0.1,
      min_child_weight = cfg_xgb$min_child_weight %||% 1L,
      subsample = cfg_xgb$subsample %||% 0.8,
      colsample_bytree = cfg_xgb$colsample_bytree %||% 0.8,
      min_split_loss = 0,
      reg_reg_lambda = 1,
      stringsAsFactors = FALSE
    )
    grid <- rbind(full_grid[sampled_idx, , drop = FALSE], baseline_row)
    grid <- grid[!duplicated(grid), ]
  } else {
    grid <- full_grid
  }
  cat("  Evaluating", nrow(grid), "candidates\n")

  grid_rows <- list()
  for (i in seq_len(nrow(grid))) {
    row <- grid[i, , drop = FALSE]
    hp <- list(
      nrounds          = row$nrounds,
      max_depth        = row$max_depth,
      learning_rate    = row$learning_rate,
      min_child_weight = row$min_child_weight,
      subsample        = row$subsample,
      colsample_bytree = row$colsample_bytree,
      min_split_loss   = row$min_split_loss,
      reg_reg_lambda   = row$reg_reg_lambda
    )
    obj_rmse <- robust_mean_rmse("XGB", hp)
    row$robust_mean_rmse <- obj_rmse
    grid_rows[[i]] <- row
    if (i %% 10 == 0 || i == nrow(grid)) {
      cat(sprintf("  [%d/%d] nrounds=%d depth=%d lr=%.2f mcw=%d ss=%.1f csbt=%.1f gamma=%g lambda=%g -> RMSE=%.6f\n",
        i, nrow(grid), row$nrounds, row$max_depth, row$learning_rate,
        row$min_child_weight, row$subsample, row$colsample_bytree,
        row$min_split_loss, row$reg_reg_lambda, obj_rmse))
    }
  }
  xgb_grid_tbl <- dplyr::bind_rows(grid_rows)
  best_idx <- which.min(xgb_grid_tbl$robust_mean_rmse)
  best_row <- xgb_grid_tbl[best_idx, , drop = FALSE]
  cat("\n  Best XGB config (robust RMSE =", round(best_row$robust_mean_rmse, 6), "):\n")
  cat("    nrounds=", best_row$nrounds, " max_depth=", best_row$max_depth,
      " lr=", best_row$learning_rate, " mcw=", best_row$min_child_weight,
      " subsample=", best_row$subsample, " colsample_bytree=", best_row$colsample_bytree,
      " gamma=", best_row$min_split_loss, " lambda=", best_row$reg_reg_lambda, "\n", sep = "")

  xgb_config <- list(
    nrounds          = best_row$nrounds,
    max_depth        = best_row$max_depth,
    learning_rate    = best_row$learning_rate,
    subsample        = best_row$subsample,
    colsample_bytree = best_row$colsample_bytree,
    min_child_weight = best_row$min_child_weight,
    min_split_loss   = best_row$min_split_loss,
    reg_reg_lambda   = best_row$reg_reg_lambda,
    robust_cv_metrics = xgb_grid_tbl
  )
  saveRDS(xgb_config, file.path(robust_dir, "best_config_xgb_robust.rds"))
  write.csv(xgb_grid_tbl, file.path(robust_dir, "robust_cv_metrics_xgb.csv"), row.names = FALSE)
}

# ---------------------------------------------------------------------------
# Tune GAM robustly (k_covariate grid)
# ---------------------------------------------------------------------------
if ("GAM" %in% models) {
  cat("\n=== Robust tuning: GAM ===\n")
  gam_k_grid <- c(2L, 4L, 6L, 8L, 10L, 12L, 15L, 20L)
  rows <- list()
  for (k_cov in gam_k_grid) {
    hp <- list(k_covariate = k_cov)
    obj_rmse <- robust_mean_rmse("GAM", hp)
    rows[[length(rows) + 1L]] <- data.frame(k_covariate = k_cov, robust_mean_rmse = obj_rmse, stringsAsFactors = FALSE)
    cat("  k_covariate=", k_cov, " -> robust mean RMSE=", round(obj_rmse, 6), "\n", sep = "")
  }
  gam_tbl <- dplyr::bind_rows(rows)
  best_idx <- which.min(gam_tbl$robust_mean_rmse)
  best_k <- gam_tbl$k_covariate[best_idx]
  gam_config <- list(k_covariate = best_k, robust_cv_metrics = gam_tbl)
  saveRDS(gam_config, file.path(robust_dir, "best_config_gam_robust.rds"))
  write.csv(gam_tbl, file.path(robust_dir, "robust_cv_metrics_gam.csv"), row.names = FALSE)
}

# ---------------------------------------------------------------------------
# Tune GPR robustly (restricted to baseline kernel + nearby nuggets)
# ---------------------------------------------------------------------------
if ("GPR" %in% models) {
  cat("\n=== Robust tuning: GPR ===\n")
  if (is.null(cfg_gpr)) stop("Missing best_config_gpr.rds for GPR baseline.")

  kernel_candidates <- c(cfg_gpr$kernel %||% "matern52")
  nug_min_vals <- c(1e-8, 1e-6, 1e-4)
  nug_max_vals <- c(10, 50, 100)

  # pick the two closest (in log space) nug_min and nug_max values
  log_base_min <- log10(cfg_gpr$nug.min %||% 1e-6)
  log_base_max <- log10(cfg_gpr$nug.max %||% 50)
  nug_min_candidates <- nug_min_vals[order(abs(log10(nug_min_vals) - log_base_min))][1:2]
  nug_max_candidates <- nug_max_vals[order(abs(log10(nug_max_vals) - log_base_max))][1:2]

  rows <- list()
  for (kernel in kernel_candidates) {
    for (nug_min in nug_min_candidates) {
      for (nug_max in nug_max_candidates) {
        if (nug_min >= nug_max) next
        hp <- list(kernel = kernel, nug.min = nug_min, nug.max = nug_max, nug.est = TRUE)
        obj_rmse <- robust_mean_rmse("GPR", hp)
        rows[[length(rows) + 1L]] <- data.frame(
          kernel = kernel, nug.min = nug_min, nug.max = nug_max,
          robust_mean_rmse = obj_rmse, stringsAsFactors = FALSE
        )
        cat("  kernel=", kernel, " nug.min=", nug_min, " nug.max=", nug_max,
            " -> robust mean RMSE=", round(obj_rmse, 6), "\n", sep = "")
      }
    }
  }
  gpr_tbl <- dplyr::bind_rows(rows)
  if (nrow(gpr_tbl) == 0L) stop("No valid GPR hyperparameter candidates after restriction.")
  best_idx <- which.min(gpr_tbl$robust_mean_rmse)
  best_row <- gpr_tbl[best_idx, , drop = FALSE]
  gpr_config <- list(
    kernel = best_row$kernel,
    nug.min = best_row$`nug.min`,
    nug.max = best_row$`nug.max`,
    nug.est = TRUE,
    robust_cv_metrics = gpr_tbl
  )
  saveRDS(gpr_config, file.path(robust_dir, "best_config_gpr_robust.rds"))
  write.csv(gpr_tbl, file.path(robust_dir, "robust_cv_metrics_gpr.csv"), row.names = FALSE)
}

cat("\nRobust tuning complete.\n  Output dir:", robust_dir, "\n")

