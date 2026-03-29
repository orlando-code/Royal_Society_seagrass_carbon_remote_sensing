## Robust hyperparameter tuning for seed-resampled `pixel_grouped`.
##
## Hyperparameter selection objective:
##   - For each candidate hyperparameter set, evaluate CV performance
##     on *multiple* pixel-grouped fold instantiations (`robust_fold_seed_list`).
##   - Select the candidate with minimum robust RMSE score:
##       mean_rmse + robust_rmse_lambda * sd_rmse
##
## This script expects robust per-model predictor sets to exist, typically from:
##   robust_shap_covariate_pruning.R
##
## Outputs:
##   output/<cv_regime_name>/cv_pipeline/robust_pixel_grouped_tuning/
##     - best_config_<model>_robust.rds
##     - robust_cv_metrics_<model>.csv
##
project_root <- here::here()
setwd(project_root)

source(file.path(project_root, "modelling/R/helpers.R"))
source(file.path(project_root, "modelling/R/assign_region_from_latlon.R"))
source(file.path(project_root, "modelling/config/pipeline_config.R"))
load_packages(c("here", "dplyr", "readr", "mgcv", "randomForest", "GauPro", "xgboost", "sf"))

cfg <- get_pipeline_config()
apply_pipeline_defaults(
  cfg,
  c(
    "cv_regime_name", "cv_type", "target_var", "log_transform_target",
    "exclude_regions", "n_folds", "cv_blocksize", "robust_fold_seed_list",
    "robust_pruned_importance_type", "correlation_filter_threshold",
    "xgb_max_grid_candidates", "robust_pruned_csv_override",
    "include_seagrass_species", "robust_rmse_lambda"
  ),
  envir = .GlobalEnv
)

cv_regime_name <- get("cv_regime_name", envir = .GlobalEnv)
cv_type <- get("cv_type", envir = .GlobalEnv)
stopifnot(identical(cv_type, "pixel_grouped"))
cv_type_hash <- "pixel_grouped"

target_var <- get("target_var", envir = .GlobalEnv)
log_transform_target <- isTRUE(get("log_transform_target", envir = .GlobalEnv))

exclude_regions <- get("exclude_regions", envir = .GlobalEnv)
n_folds <- as.integer(get("n_folds", envir = .GlobalEnv))
cv_blocksize <- get("cv_blocksize", envir = .GlobalEnv) # unused for pixel_grouped

robust_fold_seed_list <- get("robust_fold_seed_list", envir = .GlobalEnv)
robust_fold_seed_list <- as.integer(robust_fold_seed_list)
include_seagrass_species <- isTRUE(get("include_seagrass_species", envir = .GlobalEnv))
robust_rmse_lambda <- max(0, as.numeric(get("robust_rmse_lambda", envir = .GlobalEnv)))

robust_dir_default <- file.path(project_root, "output", cv_regime_name, "cv_pipeline", "robust_pixel_grouped_tuning")
robust_dir_override <- if (exists("robust_tuning_dir_override", envir = .GlobalEnv, inherits = FALSE)) {
  get("robust_tuning_dir_override", envir = .GlobalEnv)
} else {
  NA_character_
}
robust_dir <- if (!is.na(robust_dir_override) && nzchar(robust_dir_override)) robust_dir_override else robust_dir_default
dir.create(robust_dir, recursive = TRUE, showWarnings = FALSE)

cat("Robust hyperparameter tuning (pixel_grouped)\n")
cat("  robust_fold_seed_list:", paste(robust_fold_seed_list, collapse = ", "), "\n")
cat("  include_seagrass_species:", include_seagrass_species, "\n")
cat("  n_folds:", n_folds, "\n")
cat("  robust_rmse_lambda:", robust_rmse_lambda, "\n")
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
  c(
    "latitude", "longitude", "number_id_final_version",
    "seagrass_species",
    "region", target_var
  )
)
predictor_vars_full <- predictor_vars_full[predictor_vars_full %in% colnames(dat)]

species_region_cols <- intersect(
  c(if (include_seagrass_species) "seagrass_species" else character(0), "region"),
  names(dat)
)

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
robust_pruned_importance_type <- get("robust_pruned_importance_type", envir = .GlobalEnv)
robust_pruned_importance_type <- match.arg(robust_pruned_importance_type, choices = c("perm", "shap"))
robust_pruned_csv_override <- get("robust_pruned_csv_override", envir = .GlobalEnv)

robust_cov_dir <- file.path(project_root, "output", cv_regime_name, "covariate_selection", "robust_pixel_grouped")
default_robust_pruned_csv <- if (identical(robust_pruned_importance_type, "shap")) {
  file.path(
    robust_cov_dir,
    paste0("pruned_model_variables_shap_robust_pixel_grouped_seeds_", seeds_str, ".csv")
  )
} else {
  file.path(
    robust_cov_dir,
    paste0("pruned_model_variables_perm_robust_pixel_grouped_seeds_", seeds_str, ".csv")
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
  if (file.exists(fallback_csv)) {
    cat("\nWARNING: robust pruned vars not found:\n  ", robust_pruned_csv,
        "\nFalling back to non-robust pruned vars:\n  ", fallback_csv, "\n")
    robust_pruned_csv <- fallback_csv
  } else {
    cat("\nWARNING: no pruned-vars CSV found; falling back to correlation-filtered predictors for all models.\n")
  }
}

predictor_vars_by_model <- list()
if (file.exists(robust_pruned_csv)) {
  pruned_df <- read.csv(robust_pruned_csv, stringsAsFactors = FALSE)
  stopifnot(all(c("model", "variable") %in% names(pruned_df)))
  for (m in unique(pruned_df$model)) {
    vars <- intersect(pruned_df$variable[pruned_df$model == m], colnames(core_data))
    if (length(vars) >= 2) predictor_vars_by_model[[m]] <- vars
  }
} else {
  corr_vars <- prune_by_correlation(
    data = core_data,
    predictor_vars = predictor_vars_full,
    target_var = "median_carbon_density",
    cor_threshold = as.numeric(get("correlation_filter_threshold", envir = .GlobalEnv))
  )
  if (isTRUE(include_seagrass_species) && "seagrass_species" %in% names(core_data)) {
    corr_vars <- unique(c(corr_vars, "seagrass_species"))
  }
  for (m in c("GPR", "GAM", "XGB")) {
    vars <- intersect(corr_vars, colnames(core_data))
    if (length(vars) >= 2) predictor_vars_by_model[[m]] <- vars
  }
}
models <- intersect(names(predictor_vars_by_model), c("GPR", "GAM", "XGB"))
if (length(models) == 0L) stop("No models found after loading robust pruned variables.")
cat("  Robust predictor sets loaded for:", paste(models, collapse = ", "), "\n")

# ---------------------------------------------------------------------------
# Baseline configs (optional; used to enrich candidate ranges)
# ---------------------------------------------------------------------------
config_dir <- file.path(project_root, "output", cv_regime_name, "cv_pipeline")
cfg_xgb_path <- file.path(config_dir, "best_config_xgb.rds")
cfg_gam_path <- file.path(config_dir, "best_config_gam.rds")
cfg_gpr_path <- file.path(config_dir, "best_config_gpr.rds")

cfg_xgb <- if (file.exists(cfg_xgb_path)) readRDS(cfg_xgb_path) else list()
cfg_gam <- if (file.exists(cfg_gam_path)) readRDS(cfg_gam_path) else list()
cfg_gpr <- if (file.exists(cfg_gpr_path)) readRDS(cfg_gpr_path) else list()

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
# Robust objective helper: mean/sd RMSE across seeds and folds
# ---------------------------------------------------------------------------
robust_rmse_stats <- function(model_name, hp_by_model) {
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
  robust_mean_rmse <- mean(rmses, na.rm = TRUE)
  robust_sd_rmse <- stats::sd(rmses, na.rm = TRUE)
  if (!is.finite(robust_sd_rmse)) robust_sd_rmse <- 0
  robust_rmse_score <- robust_mean_rmse + robust_rmse_lambda * robust_sd_rmse
  list(
    robust_mean_rmse = robust_mean_rmse,
    robust_sd_rmse = robust_sd_rmse,
    robust_n_rmse = sum(is.finite(rmses)),
    robust_rmse_lambda = robust_rmse_lambda,
    robust_rmse_score = robust_rmse_score
  )
}

# ---------------------------------------------------------------------------
# Deterministic, representative candidate selection for large XGB grids
# ---------------------------------------------------------------------------
select_balanced_rows <- function(grid, n_select, key_cols) {
  n_select <- as.integer(n_select)
  if (n_select <= 0L) return(grid[FALSE, , drop = FALSE])
  if (nrow(grid) <= n_select) return(grid)

  ord_cols <- unique(c(key_cols, names(grid)))
  ord <- do.call(order, grid[ord_cols])
  g <- grid[ord, , drop = FALSE]

  counts <- lapply(key_cols, function(col) {
    lv <- sort(unique(g[[col]]))
    list(
      target = ceiling(n_select / length(lv)),
      counts = stats::setNames(integer(length(lv)), as.character(lv))
    )
  })
  names(counts) <- key_cols

  selected <- rep(FALSE, nrow(g))
  row_score <- function(i) {
    s <- 0
    for (col in key_cols) {
      key <- as.character(g[[col]][i])
      deficit <- counts[[col]]$target - counts[[col]]$counts[[key]]
      s <- s + max(deficit, 0)
    }
    s
  }

  for (k in seq_len(n_select)) {
    cand <- which(!selected)
    scores <- vapply(cand, row_score, numeric(1))
    best <- cand[which.max(scores)] # deterministic tie-break by order in cand
    selected[best] <- TRUE
    for (col in key_cols) {
      key <- as.character(g[[col]][best])
      counts[[col]]$counts[[key]] <- counts[[col]]$counts[[key]] + 1L
    }
  }
  g[selected, , drop = FALSE]
}

build_representative_xgb_grid <- function(full_grid, max_candidates, baseline_row, key_cols) {
  max_candidates <- as.integer(max_candidates)
  if (max_candidates <= 1L) return(baseline_row[1, , drop = FALSE])
  if (nrow(full_grid) <= max_candidates) {
    out <- unique(rbind(baseline_row, full_grid))
    return(out)
  }

  n_core <- max_candidates - 1L # keep one slot for baseline anchor
  core <- select_balanced_rows(full_grid, n_select = n_core, key_cols = key_cols)
  out <- unique(rbind(baseline_row, core))

  # Deterministic top-up in case baseline duplicated a selected row.
  if (nrow(out) < max_candidates) {
    ord <- do.call(order, full_grid[unique(c(key_cols, names(full_grid)))])
    full_ord <- full_grid[ord, , drop = FALSE]
    sig <- function(df) apply(df, 1, paste, collapse = "|")
    missing <- full_ord[!sig(full_ord) %in% sig(out), , drop = FALSE]
    need <- min(max_candidates - nrow(out), nrow(missing))
    if (need > 0L) out <- rbind(out, missing[seq_len(need), , drop = FALSE])
  }
  out
}

# ---------------------------------------------------------------------------
# Tune XGB robustly (candidate regularization only)
# ---------------------------------------------------------------------------
if ("XGB" %in% models) {
  cat("\n=== Robust tuning: XGB ===\n")

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
  cat("  Full grid size:", nrow(full_grid), "— selecting representative candidates (deterministic)\n")
  max_candidates <- as.integer(get("xgb_max_grid_candidates", envir = .GlobalEnv))
  baseline_row <- data.frame(
    nrounds = cfg_xgb$nrounds %||% 100L,
    max_depth = cfg_xgb$max_depth %||% 6L,
    learning_rate = cfg_xgb$learning_rate %||% 0.1,
    min_child_weight = cfg_xgb$min_child_weight %||% 1L,
    subsample = cfg_xgb$subsample %||% 0.8,
    colsample_bytree = cfg_xgb$colsample_bytree %||% 0.8,
    min_split_loss = cfg_xgb$min_split_loss %||% 0,
    reg_reg_lambda = cfg_xgb$reg_reg_lambda %||% 1,
    stringsAsFactors = FALSE
  )
  key_cols <- c("max_depth", "learning_rate", "min_child_weight", "reg_reg_lambda")
  grid <- build_representative_xgb_grid(
    full_grid = full_grid,
    max_candidates = max_candidates,
    baseline_row = baseline_row,
    key_cols = key_cols
  )
  for (kc in key_cols) {
    cat("    coverage ", kc, ": ", length(unique(grid[[kc]])), "/", length(unique(full_grid[[kc]])), " levels\n", sep = "")
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
    obj <- robust_rmse_stats("XGB", hp)
    row$robust_mean_rmse <- obj$robust_mean_rmse
    row$robust_sd_rmse <- obj$robust_sd_rmse
    row$robust_n_rmse <- obj$robust_n_rmse
    row$robust_rmse_lambda <- obj$robust_rmse_lambda
    row$robust_rmse_score <- obj$robust_rmse_score
    grid_rows[[i]] <- row
    if (i %% 10 == 0 || i == nrow(grid)) {
      cat(sprintf("  [%d/%d] nrounds=%d depth=%d lr=%.2f mcw=%d ss=%.1f csbt=%.1f gamma=%g lambda=%g -> meanRMSE=%.6f sdRMSE=%.6f score=%.6f\n",
        i, nrow(grid), row$nrounds, row$max_depth, row$learning_rate,
        row$min_child_weight, row$subsample, row$colsample_bytree,
        row$min_split_loss, row$reg_reg_lambda,
        obj$robust_mean_rmse, obj$robust_sd_rmse, obj$robust_rmse_score))
    }
  }
  xgb_grid_tbl <- dplyr::bind_rows(grid_rows)
  best_idx <- which.min(xgb_grid_tbl$robust_rmse_score)
  best_row <- xgb_grid_tbl[best_idx, , drop = FALSE]
  cat("\n  Best XGB config (robust score =", round(best_row$robust_rmse_score, 6),
      "; mean RMSE =", round(best_row$robust_mean_rmse, 6),
      "; sd RMSE =", round(best_row$robust_sd_rmse, 6), "):\n")
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
    obj <- robust_rmse_stats("GAM", hp)
    rows[[length(rows) + 1L]] <- data.frame(
      k_covariate = k_cov,
      robust_mean_rmse = obj$robust_mean_rmse,
      robust_sd_rmse = obj$robust_sd_rmse,
      robust_n_rmse = obj$robust_n_rmse,
      robust_rmse_lambda = obj$robust_rmse_lambda,
      robust_rmse_score = obj$robust_rmse_score,
      stringsAsFactors = FALSE
    )
    cat(
      "  k_covariate=", k_cov,
      " -> meanRMSE=", round(obj$robust_mean_rmse, 6),
      " sdRMSE=", round(obj$robust_sd_rmse, 6),
      " score=", round(obj$robust_rmse_score, 6), "\n",
      sep = ""
    )
  }
  gam_tbl <- dplyr::bind_rows(rows)
  best_idx <- which.min(gam_tbl$robust_rmse_score)
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
  kernel_candidates <- unique(c(cfg_gpr$kernel %||% "matern52", "matern52", "matern32", "gaussian"))
  nug_min_vals <- c(1e-8, 1e-6, 1e-4)
  nug_max_vals <- c(10, 50, 100)

  # If baseline exists, prefer nearby values; otherwise evaluate full defaults.
  if (!is.null(cfg_gpr$nug.min) && !is.null(cfg_gpr$nug.max)) {
    log_base_min <- log10(cfg_gpr$nug.min %||% 1e-6)
    log_base_max <- log10(cfg_gpr$nug.max %||% 50)
    nug_min_candidates <- nug_min_vals[order(abs(log10(nug_min_vals) - log_base_min))][1:2]
    nug_max_candidates <- nug_max_vals[order(abs(log10(nug_max_vals) - log_base_max))][1:2]
  } else {
    nug_min_candidates <- nug_min_vals
    nug_max_candidates <- nug_max_vals
  }

  rows <- list()
  for (kernel in kernel_candidates) {
    for (nug_min in nug_min_candidates) {
      for (nug_max in nug_max_candidates) {
        if (nug_min >= nug_max) next
        hp <- list(kernel = kernel, nug.min = nug_min, nug.max = nug_max, nug.est = TRUE)
        obj <- robust_rmse_stats("GPR", hp)
        rows[[length(rows) + 1L]] <- data.frame(
          kernel = kernel, nug.min = nug_min, nug.max = nug_max,
          robust_mean_rmse = obj$robust_mean_rmse,
          robust_sd_rmse = obj$robust_sd_rmse,
          robust_n_rmse = obj$robust_n_rmse,
          robust_rmse_lambda = obj$robust_rmse_lambda,
          robust_rmse_score = obj$robust_rmse_score,
          stringsAsFactors = FALSE
        )
        cat("  kernel=", kernel, " nug.min=", nug_min, " nug.max=", nug_max,
            " -> meanRMSE=", round(obj$robust_mean_rmse, 6),
            " sdRMSE=", round(obj$robust_sd_rmse, 6),
            " score=", round(obj$robust_rmse_score, 6), "\n", sep = "")
      }
    }
  }
  gpr_tbl <- dplyr::bind_rows(rows)
  if (nrow(gpr_tbl) == 0L) stop("No valid GPR hyperparameter candidates after restriction.")
  best_idx <- which.min(gpr_tbl$robust_rmse_score)
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

