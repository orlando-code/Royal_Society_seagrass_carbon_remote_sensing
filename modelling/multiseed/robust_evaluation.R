## Evaluate robustly-selected covariates + hyperparameters for seeded `pixel_grouped`.
##
## Uses:
##   - Robust covariate selection output from
##       output/<cv_regime_name>/covariate_selection/robust_pixel_grouped/
##       pruned_model_variables_perm_robust_pixel_grouped_seeds_<...>.csv
##   - Robust tuning output from
##       output/<cv_regime_name>/cv_pipeline/robust_pixel_grouped_tuning/
##       best_config_<model>_robust.rds
##
## Then evaluates across `eval_fold_seed_list` by re-instantiating the
## pixel-grouped fold assignments with different seeds.
##
project_root <- here::here()
setwd(project_root)

source(file.path(project_root, "modelling/R/helpers.R"))
source(file.path(project_root, "modelling/R/assign_region_from_latlon.R"))
source(file.path(project_root, "modelling/config/pipeline_config.R"))
load_packages(c("here", "dplyr", "readr", "patchwork", "mgcv", "GauPro", "xgboost", "sf", "randomForest"))

cfg <- get_pipeline_config()
apply_pipeline_defaults(
  cfg,
  c(
    "cv_regime_name", "cv_type", "target_var", "log_transform_target",
    "exclude_regions", "n_folds", "cv_blocksize",
    "robust_fold_seed_list", "eval_fold_seed_list",
    "robust_pruned_importance_type", "correlation_filter_threshold",
    "robust_pruned_csv_override", "include_seagrass_species"
  ),
  envir = .GlobalEnv
)

cv_regime_name <- get("cv_regime_name", envir = .GlobalEnv)
cv_type <- get("cv_type", envir = .GlobalEnv)
stopifnot(identical(cv_type, "pixel_grouped"))
cv_type_hash <- "pixel_grouped"

target_var <- get("target_var", envir = .GlobalEnv)
log_response <- isTRUE(get("log_transform_target", envir = .GlobalEnv))

exclude_regions <- get("exclude_regions", envir = .GlobalEnv)
n_folds <- as.integer(get("n_folds", envir = .GlobalEnv))
cv_blocksize <- get("cv_blocksize", envir = .GlobalEnv)

robust_fold_seed_list <- as.integer(get("robust_fold_seed_list", envir = .GlobalEnv))
eval_fold_seed_list <- as.integer(get("eval_fold_seed_list", envir = .GlobalEnv))
include_seagrass_species <- isTRUE(get("include_seagrass_species", envir = .GlobalEnv))
cv_type_label <- get0("cv_type_label", envir = .GlobalEnv, ifnotfound = cv_type_hash)

cat("Robust evaluation for pixel_grouped\n")
cat("  robust_fold_seed_list:", paste(robust_fold_seed_list, collapse = ", "), "\n")
cat("  eval_fold_seed_list:", paste(eval_fold_seed_list, collapse = ", "), "\n")
cat("  include_seagrass_species:", include_seagrass_species, "\n")

seeds_str <- paste(robust_fold_seed_list, collapse = "-")

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
predictor_vars_full <- predictor_vars_full[predictor_vars_full %in% names(complete_dat)]

core_data <- as.data.frame(complete_dat)
cat("  Complete-case rows:", nrow(core_data), "\n")

# Load robust pruned variables
robust_cov_dir <- file.path(project_root, "output", cv_regime_name, "covariate_selection", "robust_pixel_grouped")
robust_pruned_csv_override <- get("robust_pruned_csv_override", envir = .GlobalEnv)
robust_pruned_importance_type <- get("robust_pruned_importance_type", envir = .GlobalEnv)
robust_pruned_importance_type <- match.arg(robust_pruned_importance_type, choices = c("perm", "shap"))

default_robust_pruned_csv <- if (identical(robust_pruned_importance_type, "shap")) {
  file.path(
    robust_cov_dir,
    paste0(
      "pruned_model_variables_shap_robust_pixel_grouped_seeds_",
      seeds_str,
      ".csv"
    )
  )
} else {
  file.path(
    robust_cov_dir,
    paste0(
      "pruned_model_variables_perm_robust_pixel_grouped_seeds_",
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
  if (file.exists(fallback_csv)) {
    cat("\nWARNING: robust pruned vars not found:\n  ", robust_pruned_csv,
        "\nFalling back to:\n  ", fallback_csv, "\n")
    robust_pruned_csv <- fallback_csv
  } else {
    cat("\nWARNING: no pruned-vars CSV found; falling back to correlation-only covariates.\n")
    robust_pruned_csv <- NA_character_
  }
}

predictor_vars_by_model <- list()
if (!is.na(robust_pruned_csv) && file.exists(robust_pruned_csv)) {
  pruned_df <- read.csv(robust_pruned_csv, stringsAsFactors = FALSE)
  stopifnot(all(c("model", "variable") %in% names(pruned_df)))
  for (m in unique(pruned_df$model)) {
    vars <- intersect(pruned_df$variable[pruned_df$model == m], colnames(core_data))
    if (length(vars) >= 2L) predictor_vars_by_model[[m]] <- vars
  }
}
# Final fallback to correlation-only covariates if pruned sets are unavailable.
if (length(predictor_vars_by_model) == 0L) {
  corr_vars <- prune_by_correlation(
    data = core_data,
    predictor_vars = predictor_vars_full,
    target_var = "median_carbon_density",
    cor_threshold = as.numeric(get("correlation_filter_threshold", envir = .GlobalEnv))
  )
  if (isTRUE(include_seagrass_species) && "seagrass_species" %in% names(core_data)) {
    corr_vars <- unique(c(corr_vars, "seagrass_species"))
  }
  for (m in c("GPR", "GAM", "XGB", "LR")) {
    vars <- intersect(corr_vars, colnames(core_data))
    if (length(vars) >= 2L) predictor_vars_by_model[[m]] <- vars
  }
}
models <- intersect(names(predictor_vars_by_model), c("GPR", "GAM", "XGB", "LR"))
# eval_models <- get("eval_models", envir = .GlobalEnv)
# eval_models <- intersect(eval_models, models)
# if (length(eval_models) == 0L) stop("No usable robust pruned predictor sets found for requested eval_models.")
# models <- eval_models

cat("  Models to evaluate:", paste(models, collapse = ", "), "\n")

# Load robust hyperparameter configs
robust_tuning_dir_default <- file.path(project_root, "output", cv_regime_name, "cv_pipeline", "robust_pixel_grouped_tuning")
robust_tuning_dir_override <- if (exists("robust_tuning_dir_override", envir = .GlobalEnv, inherits = FALSE)) {
  get("robust_tuning_dir_override", envir = .GlobalEnv)
} else {
  NA_character_
}
robust_tuning_dir <- if (!is.na(robust_tuning_dir_override) && nzchar(robust_tuning_dir_override)) {
  robust_tuning_dir_override
} else {
  robust_tuning_dir_default
}
hp_bundle <- build_hyperparams_by_model(
  models = models,
  config_dir = file.path(project_root, "output", cv_regime_name, "cv_pipeline"),
  robust_config_dir = robust_tuning_dir,
  prefer_robust = TRUE,
  include_baseline = FALSE,
  include_legacy = FALSE,
  include_case_variant_robust = TRUE,
  include_only_with_config = FALSE
)
hyperparams_by_model <- hp_bundle$hyperparams_by_model

out_dir <- file.path(
  if (!is.null(get0("run_output_dir", envir = .GlobalEnv, ifnotfound = NULL))) {
    get("run_output_dir", envir = .GlobalEnv)
  } else {
    file.path(
      project_root,
      "output",
      cv_regime_name,
      "cv_pipeline",
      build_seeded_run_folder_name(
        cv_type_label = cv_type_label,
        folder_type = "evaluation",
        repeat_seed_list = eval_fold_seed_list,
        robust_seed_list = robust_fold_seed_list,
        include_seed_values = TRUE
      )
    )
  }
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
assign("run_output_dir", out_dir, envir = .GlobalEnv)

by_seed_detailed_list <- list()
by_seed_summary_list <- list()
pooled_r2_list       <- list()
pooled_counts_list  <- list()

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
    cv_method_name = paste0("pixel_grouped_eval_seed_", seed),
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

  seed_pooled_counts <- attr(res, "pooled_counts_by_model_fold")
  if (!is.null(seed_pooled_counts) && nrow(seed_pooled_counts) > 0L) {
    seed_pooled_counts$fold_seed <- seed
    pooled_counts_list[[length(pooled_counts_list) + 1L]] <- seed_pooled_counts
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

pooled_counts_all <- dplyr::bind_rows(pooled_counts_list)
if (nrow(pooled_counts_all) > 0L) {
  write.csv(
    pooled_counts_all,
    file.path(out_dir, "n_predictions_by_model_fold_seed.csv"),
    row.names = FALSE
  )
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
    sd_pooled_rmse         = sd(pooled_rmse, na.rm = TRUE),
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

