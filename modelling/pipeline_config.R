## Central pipeline configuration + helpers.
##
## Usage:
##   source("modelling/pipeline_config.R")
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

    # Seed registry (single source of truth for all seed vectors).
    # - active_robust_fold_seed_list: active default for robust selection/tuning
    # - paper_robust_fold_seed_list: frozen paper/report robust seeds (opt-in)
    # - eval_fold_seed_list: held-out evaluation seeds for final model evaluation
    # - tuning_seed_pool: candidate seed pool for tuning-seed subset sweep
    # - fold_seed_list: fold seeds for sensitivity checks
    seed_registry = list(
      active_robust_fold_seed_list = c(48L, 52L, 53L, 70L, 73L),
      paper_robust_fold_seed_list   = c(48L, 52L, 53L, 70L, 73L),
      eval_fold_seed_list           = 100L:121L,
      tuning_seed_pool              = 42L:81L,
      fold_seed_list                = 42L:62L
    ),
    # If TRUE, use seed_registry$paper_robust_fold_seed_list instead of default robust seeds.
    use_paper_seed_registry = FALSE,
    robust_pruned_importance_type = "shap",

    # Robust SHAP controls
    shap_n_points = 100L,
    shap_folds_per_seed = 50L,
    shap_max_gpr_train = 400L,
    shap_selection_coverage_grid = c(0.90, 0.95, 0.99),
    shap_selection_max_vars_grid = c(10L, 15L, 20L),

    # Pipeline toggles
    do_shap_refined_tuning = TRUE,
    do_sensitivity = TRUE,
    do_tuning_seed_sweep = FALSE,
    do_tuning_seed_sweep_refined_tuning = TRUE,  # done tuning
    do_diagnostics = TRUE,
    do_fit_final_models = TRUE,
    do_supplement = TRUE,

    # Full multiseed runs write to output/pixel_grouped_<robust_seed_subset>/.
    # This value is retained for backward compatibility but is not used for path construction.
    multiseed_run_output_id = NULL,
    # If TRUE, robust_fold_seed_list (and eval_fold_seed_list if present) are read from
    # output/tuning_seed_sweep_runs/chosen_seeds_latest.rds
    # after a tuning seed sweep (see tuning_seed_sweep.R).
    use_robust_seeds_from_tuning_sweep = FALSE,

    # Tuning sweep controls
    # tuning_seed_sweep_counts = c(2L, 5L, 10L, 15L, 20L),
    tuning_seed_sweep_counts = c(1L, 3L, 5L, 10L, 15L),
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
    # NULL -> sweep folder name uses timestamp; set a string to name the run directory.
    tuning_seed_sweep_run_id = NULL,

    # XGB candidate budget
    xgb_max_grid_candidates = 120L,

    # Sensitivity-suite defaults
    cv_types_to_check = c("pixel_grouped"),
    n_folds_list = c(1L, 3L, 5L, 10L),
    include_loo = FALSE,
    max_loo_groups = 250L,
    seed_plateau_tolerance = 0.02,
    seed_plateau_min_n = 3L,
    seed_plateau_stable_steps = 2L,

    # Diagnostics / supplement map defaults
    n_lons = 1000L,
    n_lats = 1000L,

    # Optional overrides / flags used by downstream scripts
    perf_detailed_csv_override = NA_character_,
    use_robust_final_configs = FALSE
  )

  if (length(overrides) > 0L) {
    for (nm in names(overrides)) cfg[[nm]] <- overrides[[nm]]
  }

  if (is.list(cfg$seed_registry) && length(cfg$seed_registry) > 0L) {
    sr <- cfg$seed_registry
    # Backward compatibility: allow top-level seed overrides in `overrides`,
    # but normalize them into `seed_registry` so there is one canonical source.
    if (!is.null(cfg$robust_fold_seed_list)) sr$active_robust_fold_seed_list <- as.integer(cfg$robust_fold_seed_list)
    if (!is.null(cfg$tuning_seed_pool)) sr$tuning_seed_pool <- as.integer(cfg$tuning_seed_pool)
    if (!is.null(cfg$eval_fold_seed_list)) sr$eval_fold_seed_list <- as.integer(cfg$eval_fold_seed_list)
    if (!is.null(cfg$fold_seed_list)) sr$fold_seed_list <- as.integer(cfg$fold_seed_list)

    robust_default <- if (!is.null(sr$active_robust_fold_seed_list)) {
      as.integer(sr$active_robust_fold_seed_list)
    } else {
      integer()
    }
    robust_paper <- if (!is.null(sr$paper_robust_fold_seed_list)) {
      as.integer(sr$paper_robust_fold_seed_list)
    } else {
      robust_default
    }
    eval_default <- if (!is.null(sr$eval_fold_seed_list)) {
      as.integer(sr$eval_fold_seed_list)
    } else {
      integer()
    }
    cfg$robust_fold_seed_list <- if (isTRUE(cfg$use_paper_seed_registry)) robust_paper else robust_default
    cfg$eval_fold_seed_list <- eval_default
    cfg$tuning_seed_pool <- if (!is.null(sr$tuning_seed_pool)) as.integer(sr$tuning_seed_pool) else integer()
    cfg$fold_seed_list <- if (!is.null(sr$fold_seed_list)) as.integer(sr$fold_seed_list) else integer()
    cfg$seed_registry <- sr
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

seagrass_strip_config_keys <- function(x, exclude_keys = character()) {
  if (!is.list(x) || length(exclude_keys) == 0L) return(x)
  out <- x
  drop_here <- intersect(names(out), exclude_keys)
  if (length(drop_here) > 0L) {
    out[drop_here] <- NULL
  }
  for (nm in names(out)) {
    if (is.list(out[[nm]])) {
      out[[nm]] <- seagrass_strip_config_keys(out[[nm]], exclude_keys = exclude_keys)
    }
  }
  out
}

seagrass_find_matching_configs <- function(current_cfg,
                                           search_root = "output",
                                           exclude_keys = character(),
                                           path_regex = NULL) {
  if (!dir.exists(search_root)) return(character())
  cfg_files <- list.files(
    path = search_root,
    pattern = "pipeline_config_effective\\.rds$",
    recursive = TRUE,
    full.names = TRUE
  )
  cfg_files <- cfg_files[grepl("/_run_metadata/pipeline_config_effective\\.rds$", cfg_files)]
  if (length(cfg_files) == 0L) return(character())
  if (!is.null(path_regex) && nzchar(as.character(path_regex))) {
    cfg_files <- cfg_files[grepl(path_regex, cfg_files, perl = TRUE)]
  }
  if (length(cfg_files) == 0L) return(character())

  normalized_current <- seagrass_strip_config_keys(current_cfg, exclude_keys = exclude_keys)
  matches <- character()
  for (fp in cfg_files) {
    old_cfg <- tryCatch(readRDS(fp), error = function(e) NULL)
    if (is.null(old_cfg) || !is.list(old_cfg)) next
    normalized_old <- seagrass_strip_config_keys(old_cfg, exclude_keys = exclude_keys)
    if (isTRUE(all.equal(normalized_old, normalized_current, check.attributes = FALSE))) {
      matches <- c(matches, fp)
    }
  }
  unique(matches)
}

seagrass_confirm_same_config <- function(current_cfg,
                                         search_root = "output",
                                         exclude_keys = character(),
                                         prompt_prefix = "run",
                                         path_regex = NULL) {
  matches <- seagrass_find_matching_configs(
    current_cfg = current_cfg,
    search_root = search_root,
    exclude_keys = exclude_keys,
    path_regex = path_regex
  )
  if (length(matches) == 0L) return(invisible(TRUE))

  if (!interactive()) {
    stop(
      "The current configuration is the same as ",
      paste(matches, collapse = ", "),
      ": it looks like the same code would run again. Re-run interactively to confirm yes/no.",
      call. = FALSE
    )
  }

  msg <- paste0(
    "The current configuration is the same as ",
    paste(matches, collapse = ", "),
    ": it looks like the same code would run again. Do you want to proceed? [yes/no]: "
  )
  repeat {
    ans <- tolower(trimws(readline(msg)))
    if (ans %in% c("y", "yes")) return(invisible(TRUE))
    if (ans %in% c("n", "no")) {
      stop("Stopped by user: duplicate ", prompt_prefix, " configuration.", call. = FALSE)
    }
    message("Please answer 'yes' or 'no'.")
  }
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
