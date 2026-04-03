## Central pipeline configuration + helpers.
##
## Usage:
##   source("modelling/config/pipeline_config.R")
##   cfg <- get_pipeline_config()
##   apply_pipeline_defaults(cfg, c("target_var", "n_folds"))

get_pipeline_config <- function(overrides = list()) {
  cfg <- list(
    # Plotting / display
    dpi = 300,
    show_titles = FALSE,

    # Response
    target_var = "median_carbon_density_100cm",
    log_transform_target = TRUE,
    exclude_regions = c("Black Sea"),

    # Models / pruning defaults
    model_list = c("GPR", "GAM", "XGB", "LR"),
    include_seagrass_species = TRUE,
    use_correlation_filter = TRUE,
    correlation_filter_threshold = 0.8,
    permutation_max_vars = 15L,
    n_permutations = 1L,
    permutation_coverage = 0.99,
    use_shap_per_model = TRUE,
    shap_selection_policy = "coverage_capped",

    # CV defaults
    n_folds = 5L,
    cv_type = "pixel_grouped",
    cv_blocksize = 1000L,
    cv_regime_name = "pixel_grouped",
    cv_output_dir = file.path("output", "pixel_grouped"),
    cv_type_label = "pixel_grouped",

    # Robust seed policy (5 seeds promoted from tuning_seed_sweep using RMSE-robust
    robust_fold_seed_list = c(48L, 52L, 53L, 70L, 73L),
    # robust_fold_seed_list = c(45L, 59L, 63L, 69L, 77L),
    # Ten eval seeds, disjoint from tuning_seed_pool (42:81) and robust_fold_seed_list
    # (61 is in robust; former 52:61 overlapped).
    eval_fold_seed_list = 100L:121L,
    robust_pruned_importance_type = "shap",
    # Robust RMSE objective for hyperparameter and sweep-subset selection:
    # score = mean_rmse + robust_rmse_lambda * sd_rmse.
    robust_rmse_lambda = 0, # no correction fo sd of rmse

    # Robust SHAP controls
    shap_n_points = 100L,
    shap_folds_per_seed = 50L,
    shap_max_gpr_train = 400L,
    shap_selection_coverage_grid = c(0.90, 0.95, 0.99),
    shap_selection_max_vars_grid = c(10L, 15L, 20L),

    # Pipeline toggles
    do_warm_start = FALSE,
    do_shap_refined_tuning = TRUE,
    do_sensitivity = TRUE,
    do_tuning_seed_sweep = FALSE,
    do_tuning_seed_sweep_refined_tuning = FALSE,  # done tuning
    do_diagnostics = TRUE,
    do_fit_final_models = TRUE,
    do_supplement = TRUE,

    # Tuning sweep controls
    # tuning_seed_sweep_counts = c(2L, 5L, 10L, 15L, 20L),
    tuning_seed_sweep_counts = c(1L, 3L, 5L, 10L, 15L),
    tuning_seed_pool = 42L:81L,
    tuning_sweep_eval_seed_list = 100L:121L,
    tuning_seed_sampling = "random",
    # Subset registry (see tuning_seed_sweep.R): lowering repeats reuses existing runs; raising only tops up.
    tuning_seed_sweep_repeats = 3L,
    tuning_seed_sweep_random_seed = 42L,
    tuning_seed_sweep_skip_existing = TRUE,
    # If TRUE, do not trust cached eval artifacts during sweep runs; recompute
    # at least evaluation to avoid stale/overwritten summary reuse.
    tuning_seed_sweep_force_recompute = FALSE,
    tuning_seed_sweep_unique_subsets = TRUE,
    tuning_seed_sweep_parallel_jobs = 12L,

    # XGB candidate budget
    xgb_max_grid_candidates = 120L,

    # Sensitivity-suite defaults
    cv_types_to_check = c("pixel_grouped"),
    n_folds_list = c(1L, 3L, 5L, 10L),
    fold_seed_list = 42L:62L,
    include_loo = FALSE,
    max_loo_groups = 250L,
    seed_plateau_tolerance = 0.02,
    seed_plateau_min_n = 3L,
    seed_plateau_stable_steps = 2L,
    r2_sensitivity_models = c("GPR", "GAM", "XGB", "LR"),
    fold_sensitivity_models = c("GPR", "GAM", "XGB", "LR"),

    # Diagnostics / supplement map defaults
    n_lons = 1000L,
    n_lats = 1000L,

    # Optional overrides / flags used by downstream scripts
    robust_pruned_csv_override = NA_character_,
    perf_detailed_csv_override = NA_character_,
    use_robust_final_configs = FALSE
  )

  if (length(overrides) > 0L) {
    for (nm in names(overrides)) cfg[[nm]] <- overrides[[nm]]
  }

  # Keep derived defaults consistent.
  if (is.null(cfg$cv_output_dir) || !nzchar(cfg$cv_output_dir)) {
    cfg$cv_output_dir <- file.path("output", cfg$cv_regime_name)
  }
  cfg
}

apply_pipeline_defaults <- function(cfg, keys = names(cfg), envir = .GlobalEnv) {
  keys <- intersect(keys, names(cfg))
  for (nm in keys) {
    if (!exists(nm, envir = envir, inherits = FALSE)) {
      assign(nm, cfg[[nm]], envir = envir)
    }
  }
  invisible(NULL)
}

# Build a run-folder name that encodes split regime + seed-budget geometry.
# Example:
#   pixel_grouped_evaluation_62x3_seeds_48-66-77
build_seeded_run_folder_name <- function(cv_type_label, folder_type,
                                         repeat_seed_list, robust_seed_list,
                                         include_seed_values = TRUE) {
  n_repeats <- length(as.integer(repeat_seed_list))
  n_seeds <- length(as.integer(robust_seed_list))
  base <- paste0(cv_type_label, "_", folder_type, "_", n_repeats, "x", n_seeds, "_seeds")
  if (!isTRUE(include_seed_values)) return(base)
  seed_str <- paste(as.integer(robust_seed_list), collapse = "-")
  if (!nzchar(seed_str)) return(base)
  paste0(base, "_", seed_str)
}
