## Robust SHAP-based covariate selection for seed-resampled `pixel_grouped`.
##
## To make it robust to seeded `pixel_grouped` split instantiations,
## we compute SHAP on fold-specific training sets for multiple `fold_seed`s
## and average the absolute SHAP values across folds and seeds.
##
## Output:
##   <run_output_dir>/covariate_selection/robust_pixel_grouped/ (main run), or
##   <subset_work>/covariates/ (tuning sweep subset run).
##     shap_importance_by_seed_fold_robust_pixel_grouped_seeds_<...>.csv
##     shap_importance_by_seed_robust_pixel_grouped_seeds_<...>.csv
##
if (!exists("seagrass_init_repo", mode = "function", inherits = TRUE)) source("modelling/R/init_repo.R")
project_root <- seagrass_init_repo(
  packages = c("here", "dplyr", "readr", "ggplot2", "iml", "sf"),
  source_files = c(
    "modelling/R/assign_region_from_latlon.R",
    "modelling/pipeline_config.R"
  ),
  include_helpers = TRUE,
  require_core_inputs = TRUE,
  check_renv = TRUE
)
cfg <- get_pipeline_config()
# Assign any absent global variables if they don't already exist in .GlobalEnv
apply_pipeline_defaults(
  cfg,
  c(
    "cv_regime_name", "cv_type_label", "target_var", "log_transform_target",
    "exclude_regions", "n_folds", "cv_blocksize", "robust_fold_seed_list",
    "use_correlation_filter", "correlation_filter_threshold",
    "shap_n_points", "shap_folds_per_seed", "shap_max_gpr_train",
    "permutation_coverage", "permutation_max_vars", "model_list",
    "shap_selection_policy", "shap_selection_coverage_grid", "shap_selection_max_vars_grid",
    "include_seagrass_species"
  ),
  envir = .GlobalEnv
)

cv_type_hash <- "pixel_grouped"
# cv_regime_name <- get("cv_regime_name", envir = .GlobalEnv)
# cv_type_label <- get("cv_type_label", envir = .GlobalEnv)

# target_var <- get("target_var", envir = .GlobalEnv)
# log_transform_target <- isTRUE(get("log_transform_target", envir = .GlobalEnv))

# exclude_regions <- get("exclude_regions", envir = .GlobalEnv)

# n_folds <- as.integer(get("n_folds", envir = .GlobalEnv))
# cv_blocksize <- get("cv_blocksize", envir = .GlobalEnv)

# robust_fold_seed_list <- get("robust_fold_seed_list", envir = .GlobalEnv)
# robust_fold_seed_list <- as.integer(robust_fold_seed_list)
# include_seagrass_species <- isTRUE(get("include_seagrass_species", envir = .GlobalEnv))

# use_correlation_filter <- isTRUE(get("use_correlation_filter", envir = .GlobalEnv))
# correlation_filter_threshold <- get("correlation_filter_threshold", envir = .GlobalEnv)

# # SHAP controls: defaults selected to reduce runtime.
# shap_n_points <- as.integer(get("shap_n_points", envir = .GlobalEnv))
# shap_folds_per_seed <- as.integer(get("shap_folds_per_seed", envir = .GlobalEnv))
# shap_max_gpr_train <- as.integer(get("shap_max_gpr_train", envir = .GlobalEnv))

# # Selection controls: reuse the existing tuning/pruning defaults
# permutation_coverage <- get("permutation_coverage", envir = .GlobalEnv)
# permutation_max_vars <- as.integer(get("permutation_max_vars", envir = .GlobalEnv))
# shap_selection_policy <- get("shap_selection_policy", envir = .GlobalEnv)
# shap_selection_policy <- match.arg(
#   shap_selection_policy,
#   choices = c("coverage_capped", "coverage_or_floor")
# )
# shap_selection_coverage_grid <- as.numeric(get("shap_selection_coverage_grid", envir = .GlobalEnv))
# shap_selection_max_vars_grid <- as.integer(get("shap_selection_max_vars_grid", envir = .GlobalEnv))

# model_list <- get("model_list", envir = .GlobalEnv)
model_list <- intersect(model_list, c("GPR", "GAM", "XGB", "LR"))
if (length(model_list) == 0L) stop("No supported models for SHAP pruning.")

cat("Robust SHAP covariate pruning set-up:\n")
cat("  cv_regime_name:", cv_regime_name, "\n")
cat("  cv_type label:", cv_type_label, "\n")
cat("  cv_type hash:", cv_type_hash, "\n")
cat("  robust_fold_seed_list:", paste(robust_fold_seed_list, collapse = ", "), "\n")
cat("  n_folds:", n_folds, " cv_blocksize:", cv_blocksize, "\n")
cat("  models:", paste(model_list, collapse = ", "), "\n")
cat("  SHAP: n_points=", shap_n_points, " folds_per_seed=", shap_folds_per_seed,
    " max_gpr_train=", shap_max_gpr_train, "\n")
cat("  select: max_vars=", permutation_max_vars, " coverage=", permutation_coverage, "\n")
cat("  select policy:", shap_selection_policy, "\n")
cat("  select sensitivity grid (coverage):", paste(shap_selection_coverage_grid, collapse = ", "), "\n")
cat("  select sensitivity grid (max_vars):", paste(shap_selection_max_vars_grid, collapse = ", "), "\n")
cat("  use_correlation_filter:", use_correlation_filter, " (threshold:", correlation_filter_threshold, ")\n")
cat("  include_seagrass_species:", include_seagrass_species, "\n")

dat <- readr::read_rds(file.path(project_root, "data/all_extracted_new.rds"))
if (length(exclude_regions) > 0L) {
  if (!"region" %in% names(dat)) dat <- assign_region_from_latlon(dat)
  dat <- dat[is.na(dat$region) | !dat$region %in% exclude_regions, ]
}

candidate_sets <- build_robust_candidate_predictor_sets(  # check this
  dat = dat,
  target_var = target_var,
  use_correlation_filter = use_correlation_filter,
  correlation_filter_threshold = correlation_filter_threshold,
  verbose = TRUE
)
predictor_vars_full <- candidate_sets$predictor_vars_full
cor_predictor_vars <- candidate_sets$candidate_predictor_vars

importance_predictor_vars <- cor_predictor_vars
if (isTRUE(include_seagrass_species) && "seagrass_species" %in% names(dat)) {
  importance_predictor_vars <- c(importance_predictor_vars, "seagrass_species")
}
importance_predictor_vars <- unique(importance_predictor_vars)

species_region_cols <- intersect(
  c(if (include_seagrass_species) "seagrass_species" else character(0), "region"),
  names(dat)
)

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

# Intersect between all predictor vars and the complete data frame
predictor_vars_full <- predictor_vars_full[predictor_vars_full %in% names(complete_dat)]
cor_predictor_vars <- cor_predictor_vars[cor_predictor_vars %in% names(complete_dat)]
importance_predictor_vars <- importance_predictor_vars[importance_predictor_vars %in% names(complete_dat)]

cat("\n  Complete-case rows:", nrow(complete_dat), "\n")
cat("  SHAP predictor vars:", length(importance_predictor_vars), "\n")

# ---------------------------------------------------------------------------
# Load robust hyperparameters (so that SHAP uses the same tuned settings)
# ---------------------------------------------------------------------------
if (is.na(run_output_dir) || !nzchar(as.character(run_output_dir)) || is.null(run_output_dir)) {
  stop("run_output_dir must be set in .GlobalEnv for robust SHAP pruning.")
}
subset_work_root <- if (basename(run_output_dir) == "evaluation") dirname(run_output_dir) else run_output_dir  # for run_tuning_seed_sweep.R
robust_tuning_dir <- if (basename(run_output_dir) == "evaluation") {
  file.path(subset_work_root, "tuning")
} else {
  file.path(run_output_dir, "cv_pipeline", "robust_pixel_grouped_tuning")
}
if (!dir.exists(robust_tuning_dir)) {
  stop("Required robust tuning directory not found: ", robust_tuning_dir)
}
hp_bundle <- build_hyperparams_by_model(
  models = model_list,
  config_dir = file.path(project_root, "output", cv_regime_name, "cv_pipeline"),
  robust_config_dir = robust_tuning_dir,
  prefer_robust = TRUE,
  include_baseline = TRUE,
  defaults_by_model = NULL,
  include_only_with_config = FALSE
)
hyperparams_by_model <- hp_bundle$hyperparams_by_model

# ---------------------------------------------------------------------------
# Compute SHAP importance on fold-specific training sets (robust across seeds)
# ---------------------------------------------------------------------------
robust_cov_dir <- if (basename(run_output_dir) == "evaluation") {
  file.path(subset_work_root, "covariates")
} else {
  file.path(run_output_dir, "covariate_selection", "robust_pixel_grouped")
}
dir.create(robust_cov_dir, recursive = TRUE, showWarnings = FALSE)

# Create output files signed by seed list
seeds_str <- paste(robust_fold_seed_list, collapse = "-")
out_csv <- file.path(
  robust_cov_dir,
  paste0(
    "pruned_model_variables_shap_robust_pixel_grouped_seeds_",
    seeds_str,
    ".csv"
  )
)

out_seed_fold_csv <- file.path(
  robust_cov_dir,
  paste0("shap_importance_by_seed_fold_robust_pixel_grouped_seeds_", seeds_str, ".csv")
)
out_seed_csv <- file.path(
  robust_cov_dir,
  paste0("shap_importance_by_seed_robust_pixel_grouped_seeds_", seeds_str, ".csv")
)
out_summary_csv <- file.path(
  robust_cov_dir,
  paste0("shap_importance_summary_robust_pixel_grouped_seeds_", seeds_str, ".csv")
)
out_selection_sensitivity_csv <- file.path(
  robust_cov_dir,
  paste0("shap_selection_sensitivity_robust_pixel_grouped_seeds_", seeds_str, ".csv")
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
      # Make SHAP's internal sampling deterministic per seed/fold/model for reproducibility
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

# Compute per-seed mean and SD of SHAP importance for each model
imp_by_seed <- imp_all %>%
  group_by(model, variable, fold_seed) %>%
  summarise(
    shap_importance_seed_mean = mean(shap_importance, na.rm = TRUE),
    shap_importance_seed_sd = sd(shap_importance, na.rm = TRUE),
    n_folds_used = dplyr::n(),
    .groups = "drop"
  )

# Compute per-model mean and SD of SHAP importance
imp_avg <- imp_all %>%
  group_by(model, variable) %>%
  summarise(
    shap_importance_mean = mean(shap_importance, na.rm = TRUE),
    shap_importance_sd = sd(shap_importance, na.rm = TRUE),
    shap_importance_cv = shap_importance_sd / pmax(shap_importance_mean, .Machine$double.eps),
    n_seed_fold = dplyr::n(),
    .groups = "drop"
  )

# Compute selection sensitivity grid (variables selected by 'coverage' and 'max_vars' policies)
selection_sensitivity_rows <- list()
for (model_name in model_list) {
  dfm <- imp_avg %>% filter(model == model_name) %>% select(variable, shap_importance_mean)
  if (nrow(dfm) == 0L) next
  names(dfm)[2] <- "shap_importance"
  for (covg in shap_selection_coverage_grid) {
    for (mv in shap_selection_max_vars_grid) {
      vars_sel <- select_top_env_then_species(
        df = dfm,
        value_col = "shap_importance",
        max_vars = mv,
        coverage = covg
      )
      selection_sensitivity_rows[[length(selection_sensitivity_rows) + 1L]] <- data.frame(
        model = model_name,
        coverage = covg,
        max_vars = mv,
        n_selected = length(vars_sel),
        selected_variables = paste(vars_sel, collapse = ";"),
        stringsAsFactors = FALSE
      )
    }
  }
}
selection_sensitivity_df <- dplyr::bind_rows(selection_sensitivity_rows)

# Select top environmental variables per model, preserving the species factor (which is not pruned)
select_vars_policy <- function(dfm, coverage, max_vars, policy) {
  if (nrow(dfm) == 0L) return(character(0))
  if (!identical(policy, "coverage_or_floor")) {
    return(select_top_env_then_species(
      df = dfm,
      value_col = "shap_importance",
      max_vars = max_vars,
      coverage = coverage
    ))
  }
  species_vars <- dfm$variable[is_species_var(dfm$variable)]
  env_df <- dfm[!is_species_var(dfm$variable), , drop = FALSE]
  if (nrow(env_df) == 0L) return(unique(species_vars))
  env_df <- env_df[order(-env_df$shap_importance), , drop = FALSE]
  cum <- cumsum(pmax(env_df$shap_importance, 0))
  tot <- max(cum, na.rm = TRUE)
  cov_frac <- if (tot > 0) cum / tot else seq_len(nrow(env_df)) / nrow(env_df)
  n_cov <- min(nrow(env_df), max(1L, sum(cov_frac <= coverage, na.rm = TRUE) + 1L))
  n_keep <- min(nrow(env_df), max(as.integer(max_vars), n_cov))
  unique(c(env_df$variable[seq_len(n_keep)], species_vars))
}

pruned_rows <- list()
for (model_name in model_list) {
  dfm <- imp_avg %>%
    filter(model == model_name) %>%
    select(variable, shap_importance_mean, shap_importance_sd) %>%
    rename(shap_importance = shap_importance_mean)
  v <- select_vars_policy(
    dfm = dfm,
    coverage = permutation_coverage,
    max_vars = permutation_max_vars,
    policy = shap_selection_policy
  )
  if (length(v) >= 2L) {
    pruned_rows[[model_name]] <- data.frame(model = model_name, variable = v, stringsAsFactors = FALSE)
  }
}

pruned_df <- dplyr::bind_rows(pruned_rows)
if (nrow(pruned_df) == 0L) stop("SHAP-based pruned set came out empty.")

# Write all the outputs
write.csv(imp_all, out_seed_fold_csv, row.names = FALSE)
write.csv(imp_by_seed, out_seed_csv, row.names = FALSE)
write.csv(imp_avg, out_summary_csv, row.names = FALSE)
write.csv(selection_sensitivity_df, out_selection_sensitivity_csv, row.names = FALSE)
write.csv(pruned_df, out_csv, row.names = FALSE)
cat("\nWrote SHAP robust pruned variables to:\n", out_csv, "\n", sep = "")
cat("Wrote SHAP seed/fold importances to:\n", out_seed_fold_csv, "\n", sep = "")
cat("Wrote SHAP per-seed importances to:\n", out_seed_csv, "\n", sep = "")
cat("Wrote SHAP summary importances to:\n", out_summary_csv, "\n", sep = "")
cat("Wrote SHAP selection sensitivity grid to:\n", out_selection_sensitivity_csv, "\n", sep = "")

