## Evaluate robustly-selected covariates + hyperparameters for `pixel_grouped_random`.
##
## Uses:
##   - Robust covariate selection output from
##       output/<cv_regime_name>/covariate_selection/robust_pixel_grouped_random/
##       pruned_model_variables_perm_robust_pixel_grouped_random_seeds_<...>.csv
##   - Robust tuning output from
##       output/<cv_regime_name>/cv_pipeline/robust_pixel_grouped_random_tuning/
##       best_config_<model>_robust.rds
##
## Then evaluates across `eval_fold_seed_list` by re-instantiating the
## pixel-grouped fold assignments with different seeds.
##
project_root <- here::here()
setwd(project_root)

source(file.path(project_root, "modelling/R/helpers.R"))
source(file.path(project_root, "modelling/R/assign_region_from_latlon.R"))
load_packages(c("here", "dplyr", "readr", "patchwork", "mgcv", "GauPro", "xgboost", "sf", "randomForest"))

cv_regime_name <- get0("cv_regime_name", envir = .GlobalEnv, ifnotfound = "pixel_grouped")
cv_type <- get0("cv_type", envir = .GlobalEnv, ifnotfound = "pixel_grouped_random")
stopifnot(identical(cv_type, "pixel_grouped_random"))
cv_type_hash <- "pixel_grouped"

target_var <- get0("target_var", envir = .GlobalEnv, ifnotfound = "median_carbon_density_100cm")
log_response <- isTRUE(get0("log_transform_target", envir = .GlobalEnv, ifnotfound = TRUE))

exclude_regions <- get0("exclude_regions", envir = .GlobalEnv, ifnotfound = character(0))
n_folds <- as.integer(get0("n_folds", envir = .GlobalEnv, ifnotfound = 5L))
cv_blocksize <- get0("cv_blocksize", envir = .GlobalEnv, ifnotfound = 5000L)

robust_fold_seed_list <- as.integer(get0("robust_fold_seed_list", envir = .GlobalEnv, ifnotfound = c(42L, 43L)))
eval_fold_seed_list <- as.integer(get0("eval_fold_seed_list", envir = .GlobalEnv, ifnotfound = 42L:46L))

cat("Robust evaluation for pixel_grouped_random\n")
cat("  robust_fold_seed_list:", paste(robust_fold_seed_list, collapse = ", "), "\n")
cat("  eval_fold_seed_list:", paste(eval_fold_seed_list, collapse = ", "), "\n")

seeds_str <- paste(robust_fold_seed_list, collapse = "-")

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
predictor_vars_full <- predictor_vars_full[predictor_vars_full %in% names(complete_dat)]

core_data <- as.data.frame(complete_dat)
cat("  Complete-case rows:", nrow(core_data), "\n")

# Load robust pruned variables
robust_cov_dir <- file.path(project_root, "output", cv_regime_name, "covariate_selection", "robust_pixel_grouped_random")
robust_pruned_csv_override <- get0("robust_pruned_csv_override", envir = .GlobalEnv, ifnotfound = NA_character_)
robust_pruned_importance_type <- get0("robust_pruned_importance_type", envir = .GlobalEnv, ifnotfound = "perm")
robust_pruned_importance_type <- match.arg(robust_pruned_importance_type, choices = c("perm", "shap"))

default_robust_pruned_csv <- if (identical(robust_pruned_importance_type, "shap")) {
  file.path(
    robust_cov_dir,
    paste0(
      "pruned_model_variables_shap_robust_pixel_grouped_random_seeds_",
      seeds_str,
      ".csv"
    )
  )
} else {
  file.path(
    robust_cov_dir,
    paste0(
      "pruned_model_variables_perm_robust_pixel_grouped_random_seeds_",
      seeds_str,
      ".csv"
    )
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
      "\nFalling back to:\n  ", fallback_csv, "\n")
  robust_pruned_csv <- fallback_csv
}

pruned_df <- read.csv(robust_pruned_csv, stringsAsFactors = FALSE)
stopifnot(all(c("model", "variable") %in% names(pruned_df)))

predictor_vars_by_model <- list()
for (m in unique(pruned_df$model)) {
  vars <- intersect(pruned_df$variable[pruned_df$model == m], colnames(core_data))
  if (length(vars) >= 2L) predictor_vars_by_model[[m]] <- vars
}
models <- intersect(names(predictor_vars_by_model), c("GPR", "GAM", "XGB"))
eval_models <- get0("eval_models", envir = .GlobalEnv, ifnotfound = models)
eval_models <- intersect(eval_models, models)
if (length(eval_models) == 0L) stop("No usable robust pruned predictor sets found for requested eval_models.")
models <- eval_models

cat("  Models to evaluate:", paste(models, collapse = ", "), "\n")

# Load robust hyperparameter configs
robust_tuning_dir <- file.path(project_root, "output", cv_regime_name, "cv_pipeline", "robust_pixel_grouped_random_tuning")
load_best_config <- function(model_name) {
  p <- file.path(robust_tuning_dir, paste0("best_config_", tolower(model_name), "_robust.rds"))
  if (file.exists(p)) return(readRDS(p))
  # fallback naming (just in case)
  p2 <- file.path(robust_tuning_dir, paste0("best_config_", model_name, "_robust.rds"))
  if (file.exists(p2)) return(readRDS(p2))
  NULL
}

hyperparams_by_model <- list()
for (m in models) {
  cfg <- load_best_config(m)
  if (is.null(cfg)) stop("Missing robust tuning config for ", m, " in ", robust_tuning_dir)
  if (m == "XGB") {
    hyperparams_by_model[[m]] <- list(
      nrounds = cfg$nrounds,
      max_depth = cfg$max_depth,
      learning_rate = cfg$learning_rate,
      subsample = cfg$subsample,
      colsample_bytree = cfg$colsample_bytree,
      min_child_weight = cfg$min_child_weight,
      min_split_loss = cfg$min_split_loss,
      reg_reg_lambda = cfg$reg_reg_lambda
    )
  } else if (m == "GAM") {
    hyperparams_by_model[[m]] <- list(k_covariate = cfg$k_covariate %||% cfg$k_covariate)
  } else if (m == "GPR") {
    hyperparams_by_model[[m]] <- list(
      kernel = cfg$kernel,
      nug.min = cfg$nug.min,
      nug.max = cfg$nug.max,
      nug.est = TRUE
    )
  }
}

out_dir <- file.path(
  project_root,
  "output",
  cv_regime_name,
  "cv_pipeline",
  paste0(
    "robust_pixel_grouped_random_evaluation_robustSeeds_",
    seeds_str
  )
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

by_seed_detailed_list <- list()
by_seed_summary_list <- list()
pooled_r2_list       <- list()

for (seed in eval_fold_seed_list) {
  cat("\n--- eval fold_seed =", seed, " ---\n")
  fold_info <- make_cv_folds(
    dat = core_data,
    covariate_cols = predictor_vars_full,
    n_folds = n_folds,
    cv_type = cv_type_hash,
    cv_blocksize = cv_blocksize,
    exclude_regions = exclude_regions,
    cache_tag = paste0("robust_eval_", cv_type, "_seed_", seed, "_", seeds_str),
    seed = seed
  )
  fold_indices <- fold_info$fold_indices

  res <- run_cv(
    cv_method_name = paste0("pixel_grouped_random_eval_seed_", seed),
    fold_indices = fold_indices,
    core_data = core_data,
    predictor_vars = predictor_vars_full,
    predictor_vars_by_model = predictor_vars_by_model,
    hyperparams_by_model = hyperparams_by_model,
    tune_hyperparams = FALSE,
    nested_tuning = FALSE,
    verbose = FALSE,
    return_predictions = FALSE,
    models = models,
    log_response = log_response
  )
  if (is.null(res) || nrow(res) == 0L) next

  seed_pooled_r2 <- attr(res, "pooled_r2")
  if (!is.null(seed_pooled_r2) && nrow(seed_pooled_r2) > 0L) {
    seed_pooled_r2$fold_seed <- seed
    pooled_r2_list[[length(pooled_r2_list) + 1L]] <- seed_pooled_r2
  }

  by_seed_detailed <- res %>%
    dplyr::mutate(fold_seed = seed, robust_fold_seed_list = seeds_str, cv_regime = cv_regime_name)
  by_seed_detailed_list[[length(by_seed_detailed_list) + 1L]] <- by_seed_detailed

  by_seed_summary <- by_seed_detailed %>%
    dplyr::group_by(model) %>%
    dplyr::summarise(
      n_folds = dplyr::n(),
      mean_r2 = mean(r2, na.rm = TRUE),
      sd_r2 = sd(r2, na.rm = TRUE),
      mean_rmse = mean(rmse, na.rm = TRUE),
      sd_rmse = sd(rmse, na.rm = TRUE),
      mean_bias = mean(bias, na.rm = TRUE),
      sd_bias = sd(bias, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(fold_seed = seed, robust_fold_seed_list = seeds_str)

  if (!is.null(seed_pooled_r2) && nrow(seed_pooled_r2) > 0L) {
    by_seed_summary <- by_seed_summary %>%
      dplyr::left_join(
        seed_pooled_r2 %>% dplyr::select(model, pooled_r2, pooled_rmse),
        by = "model"
      )
  }
  by_seed_summary_list[[length(by_seed_summary_list) + 1L]] <- by_seed_summary
}

by_seed_detailed_all <- dplyr::bind_rows(by_seed_detailed_list)
by_seed_summary_all <- dplyr::bind_rows(by_seed_summary_list)
pooled_r2_all <- dplyr::bind_rows(pooled_r2_list)

write.csv(by_seed_detailed_all, file.path(out_dir, "by_seed_detailed.csv"), row.names = FALSE)
write.csv(by_seed_summary_all, file.path(out_dir, "by_seed_summary.csv"), row.names = FALSE)
if (nrow(pooled_r2_all) > 0L) {
  write.csv(pooled_r2_all, file.path(out_dir, "pooled_r2_by_seed.csv"), row.names = FALSE)
}

has_pooled <- "pooled_r2" %in% names(by_seed_summary_all)
if (!has_pooled) {
  by_seed_summary_all$pooled_r2   <- NA_real_
  by_seed_summary_all$pooled_rmse <- NA_real_
}

across_seeds_summary <- by_seed_summary_all %>%
  dplyr::group_by(model) %>%
  dplyr::summarise(
    mean_mean_rmse         = mean(mean_rmse, na.rm = TRUE),
    sd_mean_rmse           = sd(mean_rmse, na.rm = TRUE),
    mean_mean_r2           = mean(mean_r2, na.rm = TRUE),
    sd_mean_r2             = sd(mean_r2, na.rm = TRUE),
    mean_pooled_r2         = mean(pooled_r2, na.rm = TRUE),
    sd_pooled_r2           = sd(pooled_r2, na.rm = TRUE),
    mean_pooled_rmse       = mean(pooled_rmse, na.rm = TRUE),
    mean_bias              = mean(mean_bias, na.rm = TRUE),
    sd_bias                = sd(mean_bias, na.rm = TRUE),
    neg_r2_fraction        = mean(mean_r2 < 0, na.rm = TRUE),
    neg_pooled_r2_fraction = mean(pooled_r2 < 0, na.rm = TRUE),
    min_mean_r2            = min(mean_r2, na.rm = TRUE),
    min_pooled_r2          = min(pooled_r2, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(across_seeds_summary, file.path(out_dir, "across_seeds_summary.csv"), row.names = FALSE)

cat("\nWrote robust evaluation outputs to:\n", out_dir, "\n", sep = "")

