## Train-test fraction diagnostic for seed-resampled `pixel_grouped`.
##
## For each `fold_seed` and each fold number, compute:
##   - n_test, n_train
##   - test_fraction (n_test / N)
##   - y_sd in the test fold (sd(observed))
##
## Then join to per-fold performance from robust evaluation output
## (or legacy repeated seed-check output, if provided) to see whether
## folds with small test fraction / low y variance are driving negative R².
##
if (!exists("seagrass_init_repo", mode = "function", inherits = TRUE)) {
  init_path <- file.path("modelling", "R", "init_repo.R")
  if (!file.exists(init_path)) {
    ff <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
    if (!length(ff)) stop("Run from repo root or with: Rscript /path/to/this_script.R", call. = FALSE)
    script_path <- normalizePath(sub("^--file=", "", ff[[1]]), winslash = "/", mustWork = FALSE)
    init_path <- normalizePath(file.path(dirname(script_path), "..", "R", "init_repo.R"), winslash = "/", mustWork = FALSE)
  }
  if (!file.exists(init_path)) stop("Missing bootstrap helper: modelling/R/init_repo.R", call. = FALSE)
  sys.source(init_path, envir = .GlobalEnv)
}
project_root <- seagrass_init_repo(
  include_helpers = FALSE,
  require_core_inputs = FALSE,
  check_renv = FALSE
)
project_root <- getwd()

source(file.path(project_root, "modelling/R/helpers.R"))
source(file.path(project_root, "modelling/R/assign_region_from_latlon.R"))
source(file.path(project_root, "modelling/pipeline_config.R"))
load_packages(c("here", "dplyr", "readr", "ggplot2", "patchwork"))

cfg <- get_pipeline_config()
apply_pipeline_defaults(
  cfg,
  c(
    "cv_regime_name", "target_var", "log_transform_target", "exclude_regions",
    "n_folds", "cv_type", "fold_seed_list", "cv_blocksize", "robust_fold_seed_list",
    "perf_detailed_csv_override", "include_seagrass_species"
  ),
  envir = .GlobalEnv
)

cv_regime_name <- get("cv_regime_name", envir = .GlobalEnv)
target_var <- get("target_var", envir = .GlobalEnv)
log_transform_target <- isTRUE(get("log_transform_target", envir = .GlobalEnv))

exclude_regions <- get("exclude_regions", envir = .GlobalEnv)
n_folds <- as.integer(get("n_folds", envir = .GlobalEnv))
cv_type <- get("cv_type", envir = .GlobalEnv)
cv_type_hash <- if (cv_type %in% c("pixel_grouped")) "pixel_grouped" else cv_type
if (!identical(cv_type_hash, "pixel_grouped")) {
  stop("This diagnostic currently supports only pixel-grouped hashing; got cv_type=", cv_type)
}

fold_seed_list <- get("fold_seed_list", envir = .GlobalEnv)
fold_seed_list <- as.integer(fold_seed_list)
robust_fold_seed_list <- as.integer(get("robust_fold_seed_list", envir = .GlobalEnv))
include_seagrass_species <- isTRUE(get("include_seagrass_species", envir = .GlobalEnv))
cv_type_label <- get0("cv_type_label", envir = .GlobalEnv, ifnotfound = cv_type_hash)
eval_seed_list <- if (length(fold_seed_list) > 0L) fold_seed_list else robust_fold_seed_list

run_output_dir <- if (exists("run_output_dir", envir = .GlobalEnv, inherits = FALSE)) {
  get("run_output_dir", envir = .GlobalEnv)
} else {
  file.path(
    project_root, "output", cv_regime_name, "cv_pipeline",
    build_seeded_run_folder_name(
      cv_type_label = cv_type_label,
      folder_type = "evaluation",
      repeat_seed_list = eval_seed_list,
      robust_seed_list = robust_fold_seed_list,
      include_seed_values = TRUE
    )
  )
}

cat("Train-test fraction diagnostic: pixel_grouped\n")
cat("  n_folds:", n_folds, "\n")
cat("  fold_seed_list:", paste(fold_seed_list, collapse = ", "), "\n")

# Data
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

n_total <- nrow(core_data)
stopifnot(n_total > 20L)
cat("  Complete-case rows:", n_total, "\n")

out_dir <- file.path(run_output_dir, "train_test_fraction_diagnostic")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fold_stats_rows <- list()
for (seed in fold_seed_list) {
  cat("  seed", seed, "...\n")
  fold_info <- make_cv_folds(
    dat = core_data,
    covariate_cols = predictor_vars_full,
    n_folds = n_folds,
    cv_type = cv_type_hash,
    cv_blocksize = get("cv_blocksize", envir = .GlobalEnv),
    exclude_regions = exclude_regions,
    cache_tag = paste0("ttfrac_", cv_type, "_seed_", seed),
    seed = seed
  )
  fold_indices <- fold_info$fold_indices

  for (fold in seq_len(max(fold_indices))) {
    test_idx <- which(fold_indices == fold)
    y <- core_data$median_carbon_density[test_idx]
    fold_stats_rows[[length(fold_stats_rows) + 1L]] <- data.frame(
      fold_seed = seed,
      fold = fold,
      n_test = length(test_idx),
      n_train = n_total - length(test_idx),
      test_fraction = length(test_idx) / n_total,
      y_mean = mean(y, na.rm = TRUE),
      y_sd = sd(y, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
}

fold_stats <- dplyr::bind_rows(fold_stats_rows)
write.csv(fold_stats, file.path(out_dir, "fold_stats_by_seed_by_fold.csv"), row.names = FALSE)

# Join to repeated CV per-fold performance (deployment check)
perf_csv_override <- get("perf_detailed_csv_override", envir = .GlobalEnv)
robust_eval_default <- file.path(run_output_dir, "by_seed_detailed.csv")
legacy_eval_default <- if (length(robust_fold_seed_list) > 0L) {
  file.path(
    project_root, "output", cv_regime_name, "cv_pipeline",
    paste0("robust_pixel_grouped_evaluation_robustSeeds_", paste(robust_fold_seed_list, collapse = "-")),
    "by_seed_detailed.csv"
  )
} else ""
legacy_default <- file.path(project_root, "output", cv_regime_name, "cv_pipeline", "repeated_pixel_grouped_seed_check", "by_seed_detailed.csv")
perf_csv_default <- if (nzchar(robust_eval_default) && file.exists(robust_eval_default)) {
  robust_eval_default
} else if (nzchar(legacy_eval_default) && file.exists(legacy_eval_default)) {
  legacy_eval_default
} else {
  legacy_default
}
perf_csv <- if (!is.na(perf_csv_override) && nzchar(perf_csv_override)) perf_csv_override else perf_csv_default
if (file.exists(perf_csv)) {
  perf <- readr::read_csv(perf_csv, show_col_types = FALSE)
  perf$fold_seed <- as.integer(perf$fold_seed)

  joined <- perf %>%
    left_join(fold_stats, by = c("fold_seed", "fold")) %>%
    filter(model %in% c("XGB", "GPR", "GAM"))

  summary_by_model <- joined %>%
    group_by(model) %>%
    summarise(
      n = n(),
      mean_test_fraction = mean(test_fraction, na.rm = TRUE),
      sd_test_fraction = sd(test_fraction, na.rm = TRUE),
      cor_test_fraction_r2 = cor(test_fraction, r2, use = "complete.obs"),
      cor_y_sd_r2 = cor(y_sd, r2, use = "complete.obs"),
      neg_r2_fraction = mean(r2 < 0, na.rm = TRUE),
      .groups = "drop"
    )

  write.csv(summary_by_model, file.path(out_dir, "summary_by_model.csv"), row.names = FALSE)

  cat("\nSummary by model (joined to repeated CV per-fold performance):\n")
  print(summary_by_model)

  p_r2_vs_y_sd <- ggplot(joined, aes(x = y_sd, y = r2, colour = model)) +
    geom_point(size = 2.5, alpha = 0.85) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
    theme_minimal(base_size = 12) +
    labs(title = "R2 vs test-fold y_sd (pixel_grouped)", x = "sd(observed) per fold", y = "R2 per fold")

  p_r2_vs_test_frac <- ggplot(joined, aes(x = test_fraction, y = r2, colour = model)) +
    geom_point(size = 2.5, alpha = 0.85) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
    theme_minimal(base_size = 12) +
    labs(title = "R2 vs test fold fraction (pixel_grouped)", x = "test_fraction", y = "R2 per fold")

  png_path <- file.path(out_dir, "r2_vs_y_sd_and_test_fraction.png")
  ggsave(png_path, patchwork::wrap_plots(p_r2_vs_y_sd, p_r2_vs_test_frac, ncol = 2),
         width = 12, height = 5, dpi = 200)

  cat("\nWrote:\n",
      file.path(out_dir, "fold_stats_by_seed_by_fold.csv"), "\n",
      file.path(out_dir, "summary_by_model.csv"), "\n",
      file.path(out_dir, "r2_vs_y_sd_and_test_fraction.png"), "\n",
      sep = "")
} else {
  cat("\nWARNING: No repeated check file found; skipping join to performance.\n  ", perf_csv, "\n")
  cat("\nWrote:\n",
      file.path(out_dir, "fold_stats_by_seed_by_fold.csv"), "\n",
      sep = "")
}

