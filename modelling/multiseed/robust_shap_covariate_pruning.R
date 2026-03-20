## Robust SHAP-based covariate selection for `pixel_grouped_random`.
##
## This replaces permutation-based importance with SHAP (iml::Shapley).
## To make it robust to `pixel_grouped_random` split instantiations,
## we compute SHAP on fold-specific training sets for multiple `fold_seed`s
## and average the absolute SHAP values across folds and seeds.
##
## Output:
##   output/<cv_regime_name>/covariate_selection/robust_pixel_grouped_random/
##     pruned_model_variables_shap_robust_pixel_grouped_random_seeds_<...>.csv
##
project_root <- here::here()
setwd(project_root)

source(file.path(project_root, "modelling/R/helpers.R"))
source(file.path(project_root, "modelling/R/assign_region_from_latlon.R"))
load_packages(c("dplyr", "readr", "ggplot2", "iml", "sf", "here"))

cv_regime_name <- get0("cv_regime_name", envir = .GlobalEnv, ifnotfound = "pixel_grouped")
cv_type_hash <- "pixel_grouped"
cv_type_label <- get0("cv_type_label", envir = .GlobalEnv, ifnotfound = "pixel_grouped_random")

target_var <- get0("target_var", envir = .GlobalEnv, ifnotfound = "median_carbon_density_100cm")
log_transform_target <- isTRUE(get0("log_transform_target", envir = .GlobalEnv, ifnotfound = TRUE))

exclude_regions <- get0("exclude_regions", envir = .GlobalEnv, ifnotfound = character(0))

n_folds <- as.integer(get0("n_folds", envir = .GlobalEnv, ifnotfound = 5L))
cv_blocksize <- get0("cv_blocksize", envir = .GlobalEnv, ifnotfound = 5000L)

robust_fold_seed_list <- get0("robust_fold_seed_list", envir = .GlobalEnv, ifnotfound = c(42L, 43L))
robust_fold_seed_list <- as.integer(robust_fold_seed_list)

use_correlation_filter <- isTRUE(get0("use_correlation_filter", envir = .GlobalEnv, ifnotfound = TRUE))
correlation_filter_threshold <- get0("correlation_filter_threshold", envir = .GlobalEnv, ifnotfound = 0.8)

# SHAP controls: defaults tuned to reduce runtime.
shap_n_points <- as.integer(get0("shap_n_points", envir = .GlobalEnv, ifnotfound = 20L))
shap_folds_per_seed <- as.integer(get0("shap_folds_per_seed", envir = .GlobalEnv, ifnotfound = 2L))
shap_max_gpr_train <- as.integer(get0("shap_max_gpr_train", envir = .GlobalEnv, ifnotfound = 200L))

# Selection controls: reuse the existing tuning/pruning defaults
permutation_coverage <- get0("permutation_coverage", envir = .GlobalEnv, ifnotfound = 0.99)
permutation_max_vars <- as.integer(get0("permutation_max_vars", envir = .GlobalEnv, ifnotfound = 15L))

model_list <- get0("model_list", envir = .GlobalEnv, ifnotfound = c("GPR", "GAM", "XGB"))
model_list <- intersect(model_list, c("GPR", "GAM", "XGB"))
if (length(model_list) == 0L) stop("No supported models for SHAP pruning.")

cat("Robust SHAP covariate pruning\n")
cat("  cv_regime_name:", cv_regime_name, "\n")
cat("  cv_type label:", cv_type_label, "\n")
cat("  cv_type hash:", cv_type_hash, "\n")
cat("  robust_fold_seed_list:", paste(robust_fold_seed_list, collapse = ", "), "\n")
cat("  n_folds:", n_folds, " cv_blocksize:", cv_blocksize, "\n")
cat("  models:", paste(model_list, collapse = ", "), "\n")
cat("  SHAP: n_points=", shap_n_points, " folds_per_seed=", shap_folds_per_seed,
    " max_gpr_train=", shap_max_gpr_train, "\n")
cat("  select: max_vars=", permutation_max_vars, " coverage=", permutation_coverage, "\n")
cat("  use_correlation_filter:", use_correlation_filter, " (threshold:", correlation_filter_threshold, ")\n")

dat <- readr::read_rds(file.path(project_root, "data/all_extracted_new.rds"))
if (length(exclude_regions) > 0L) {
  if (!"region" %in% names(dat)) dat <- assign_region_from_latlon(dat)
  dat <- dat[is.na(dat$region) | !dat$region %in% exclude_regions, ]
}

# Candidate predictors (environmental) - species is treated as an unprunable factor.
predictor_vars_full <- setdiff(
  colnames(dat),
  c("latitude", "longitude", "number_id_final_version", "seagrass_species", "region", target_var)
)
predictor_vars_full <- predictor_vars_full[predictor_vars_full %in% colnames(dat)]

cor_predictor_vars <- predictor_vars_full
if (isTRUE(use_correlation_filter)) {
  cat("\nPruning predictors by correlation (threshold:", correlation_filter_threshold, ") ...\n")
  cor_predictor_vars <- prune_by_correlation(
    data = dat,
    predictor_vars = cor_predictor_vars,
    target_var = target_var,
    cor_threshold = correlation_filter_threshold
  )
  cor_predictor_vars <- cor_predictor_vars[cor_predictor_vars %in% colnames(dat)]
  cat("  After correlation filter:", length(cor_predictor_vars), "predictors.\n")
}

importance_predictor_vars <- cor_predictor_vars
if ("seagrass_species" %in% names(dat)) importance_predictor_vars <- c(importance_predictor_vars, "seagrass_species")
importance_predictor_vars <- unique(importance_predictor_vars)

species_region_cols <- intersect(c("seagrass_species", "region"), names(dat))

# Build complete-case frame for consistent folds and SHAP inputs.
complete_dat <- dat %>%
  dplyr::select(
    longitude,
    latitude,
    dplyr::all_of(target_var),
    dplyr::all_of(unique(c(predictor_vars_full, cor_predictor_vars, importance_predictor_vars))),
    dplyr::all_of(species_region_cols)
  ) %>%
  dplyr::filter(complete.cases(.))

if (nrow(complete_dat) < 20L) stop("Too few complete cases for SHAP pruning.")
complete_dat$median_carbon_density <- complete_dat[[target_var]]

predictor_vars_full <- predictor_vars_full[predictor_vars_full %in% names(complete_dat)]
cor_predictor_vars <- cor_predictor_vars[cor_predictor_vars %in% names(complete_dat)]
importance_predictor_vars <- importance_predictor_vars[importance_predictor_vars %in% names(complete_dat)]

cat("\n  Complete-case rows:", nrow(complete_dat), "\n")
cat("  SHAP predictor vars (incl. species factor if present):", length(importance_predictor_vars), "\n")

# ---------------------------------------------------------------------------
# Load robust hyperparameters (so SHAP uses the same tuned settings)
# ---------------------------------------------------------------------------
robust_tuning_dir <- file.path(project_root, "output", cv_regime_name, "cv_pipeline", "robust_pixel_grouped_random_tuning")
load_best_config <- function(model_name) {
  p <- file.path(robust_tuning_dir, paste0("best_config_", tolower(model_name), "_robust.rds"))
  if (file.exists(p)) return(readRDS(p))
  # fallback: non-robust baseline
  p2 <- file.path(project_root, "output", cv_regime_name, "cv_pipeline", paste0("best_config_", tolower(model_name), ".rds"))
  if (file.exists(p2)) return(readRDS(p2))
  NULL
}

hyperparams_by_model <- list()
for (m in model_list) {
  cfg <- load_best_config(m)
  if (is.null(cfg)) stop("Missing best config for ", m, " in ", robust_tuning_dir, " (or baseline cv_pipeline).")
  if (m == "XGB") {
    hyperparams_by_model[[m]] <- list(
      nrounds = cfg$nrounds,
      max_depth = cfg$max_depth,
      learning_rate = cfg$learning_rate,
      subsample = cfg$subsample,
      colsample_bytree = cfg$colsample_bytree,
      min_child_weight = cfg$min_child_weight,
      min_split_loss = cfg$min_split_loss %||% 0,
      reg_reg_lambda = cfg$reg_reg_lambda %||% 1
    )
  } else if (m == "GAM") {
    hyperparams_by_model[[m]] <- list(k_covariate = cfg$k_covariate %||% 6L)
  } else if (m == "GPR") {
    hyperparams_by_model[[m]] <- list(
      kernel = cfg$kernel %||% "matern52",
      nug.min = cfg$nug.min %||% 1e-8,
      nug.max = cfg$nug.max %||% 50,
      nug.est = TRUE
    )
  }
}

# ---------------------------------------------------------------------------
# Compute SHAP importance on fold-specific training sets (robust across seeds)
# ---------------------------------------------------------------------------
robust_cov_dir <- file.path(project_root, "output", cv_regime_name, "covariate_selection", "robust_pixel_grouped_random")
dir.create(robust_cov_dir, recursive = TRUE, showWarnings = FALSE)

seeds_str <- paste(robust_fold_seed_list, collapse = "-")
out_csv <- file.path(
  robust_cov_dir,
  paste0(
    "pruned_model_variables_shap_robust_pixel_grouped_random_seeds_",
    seeds_str,
    ".csv"
  )
)

per_seed_model_fold_imp <- list() # key: model_seed_fold

for (seed in robust_fold_seed_list) {
  cat("\n--- fold_seed =", seed, " ---\n")
  pf <- make_cv_folds(
    dat = complete_dat,
    covariate_cols = predictor_vars_full,
    n_folds = n_folds,
    cv_type = cv_type_hash,
    cv_blocksize = cv_blocksize,
    exclude_regions = exclude_regions,
    cache_tag = paste0("robust_shap_pruning_", cv_type_hash, "_seed_", seed),
    seed = seed
  )
  fold_indices <- pf$fold_indices
  nf <- max(fold_indices)
  folds_to_use <- seq_len(min(nf, shap_folds_per_seed))

  for (model_name in model_list) {
    hp <- hyperparams_by_model[[model_name]]
    pvars <- importance_predictor_vars
    # Ensure SHAP sees only available predictor columns
    pvars <- pvars[pvars %in% names(complete_dat)]
    if (length(pvars) < 2L) next

    for (fold_k in folds_to_use) {
      train_df <- complete_dat[fold_indices != fold_k, , drop = FALSE]
      if (nrow(train_df) < 30L) next
      # Make SHAP's internal sampling deterministic per seed/fold/model
      set.seed(as.integer(seed + fold_k))

      imp <- tryCatch(
        compute_shap_importance(
          core_data = train_df,
          target_var = target_var,
          predictor_vars = pvars,
          model_name = model_name,
          n_points = shap_n_points,
          max_gpr_train = shap_max_gpr_train,
          log_response = log_transform_target,
          hyperparams = hp
        ),
        error = function(e) {
          cat("  SHAP FAILED:", model_name, " seed=", seed, " fold=", fold_k, " -> ", conditionMessage(e), "\n")
          NULL
        }
      )
      if (is.null(imp) || nrow(imp) == 0L) next
      imp$model <- model_name
      imp$fold_seed <- seed
      imp$fold <- fold_k
      per_seed_model_fold_imp[[paste(model_name, seed, fold_k, sep = "_")]] <- imp
    }
  }
}

if (length(per_seed_model_fold_imp) == 0L) stop("No SHAP importances computed; check SHAP settings/models.")

imp_all <- dplyr::bind_rows(per_seed_model_fold_imp)

imp_avg <- imp_all %>%
  group_by(model, variable) %>%
  summarise(
    shap_importance_mean = mean(shap_importance, na.rm = TRUE),
    shap_importance_sd = sd(shap_importance, na.rm = TRUE),
    .groups = "drop"
  )

# Select top env vars per model, preserving the species factor (not pruned).
pruned_rows <- list()
for (model_name in model_list) {
  dfm <- imp_avg %>% filter(model == model_name) %>% select(variable, shap_importance_mean, shap_importance_sd)
  names(dfm)[names(dfm) == "shap_importance_mean"] <- "value"
  v <- select_top_env_then_species(
    df = dfm %>% rename(shap_importance = value),
    value_col = "shap_importance",
    max_vars = permutation_max_vars,
    coverage = permutation_coverage
  )
  if (length(v) >= 2L) {
    pruned_rows[[model_name]] <- data.frame(model = model_name, variable = v, stringsAsFactors = FALSE)
  }
}

pruned_df <- dplyr::bind_rows(pruned_rows)
if (nrow(pruned_df) == 0L) stop("SHAP-based pruned set came out empty.")

write.csv(pruned_df, out_csv, row.names = FALSE)
cat("\nWrote SHAP robust pruned variables to:\n", out_csv, "\n", sep = "")

