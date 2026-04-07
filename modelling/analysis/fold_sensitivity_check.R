## Fold sensitivity check: evaluate how mean±sd performance changes
## when varying n_folds and re-instantiating grouped folds via different seeds.
##
## Scope (default):
##   - splitting regimes: random_split, location_grouped_random, pixel_grouped
##   - n_folds: 3, 5, 7
##   - seed_list: 42, 43, 44 (controls random assignment for each regime)
##
## Important:
##   This is a *diagnostic* check to justify fold count and quantify evaluation
##   uncertainty. It does NOT retune covariates/hyperparameters; it holds the
##   already selected model configs constant (from output/<cv_regime_name>).
##
## Outputs:
##   output/<cv_regime_name>/cv_pipeline/fold_sensitivity_check/
##     - fold_sensitivity_by_fold.csv
##     - fold_sensitivity_summary_by_seed.csv
##     - fold_sensitivity_summary_across_seeds.csv
##     - fold_sensitivity_rmse_plot.png
##     - fold_sensitivity_r2_plot.png

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
source(file.path(project_root, "modelling/pipeline_config.R"))
load_packages(c("here", "dplyr", "readr", "ggplot2", "tidyr"))

cfg <- get_pipeline_config()
apply_pipeline_defaults(
  cfg,
  c(
    "cv_regime_name", "target_var", "log_transform_target", "exclude_regions",
    "fold_seed_list", "n_folds_list", "cv_types_to_check",
    "use_shap_per_model", "cv_blocksize", "fold_sensitivity_models",
    "include_seagrass_species"
  ),
  envir = .GlobalEnv
)

# Guard against leaking loop variables (e.g. n_folds/cv_type) into .GlobalEnv
# when sourced from a long-lived R session.
had_global_n_folds <- exists("n_folds", envir = .GlobalEnv, inherits = FALSE)
old_global_n_folds <- if (had_global_n_folds) get("n_folds", envir = .GlobalEnv, inherits = FALSE) else NULL
had_global_cv_type <- exists("cv_type", envir = .GlobalEnv, inherits = FALSE)
old_global_cv_type <- if (had_global_cv_type) get("cv_type", envir = .GlobalEnv, inherits = FALSE) else NULL
on.exit({
  if (had_global_n_folds) {
    assign("n_folds", old_global_n_folds, envir = .GlobalEnv)
  } else if (exists("n_folds", envir = .GlobalEnv, inherits = FALSE)) {
    rm("n_folds", envir = .GlobalEnv)
  }

  if (had_global_cv_type) {
    assign("cv_type", old_global_cv_type, envir = .GlobalEnv)
  } else if (exists("cv_type", envir = .GlobalEnv, inherits = FALSE)) {
    rm("cv_type", envir = .GlobalEnv)
  }
}, add = TRUE)

cv_regime_name <- get("cv_regime_name", envir = .GlobalEnv)
target_var <- get("target_var", envir = .GlobalEnv)
log_response <- isTRUE(get("log_transform_target", envir = .GlobalEnv))

exclude_regions <- get("exclude_regions", envir = .GlobalEnv)

seed_list <- get("fold_seed_list", envir = .GlobalEnv)
n_folds_list <- get("n_folds_list", envir = .GlobalEnv)

cv_types_to_check <- get("cv_types_to_check", envir = .GlobalEnv)

# Only use models with tuned configs available.
use_shap_per_model <- isTRUE(get("use_shap_per_model", envir = .GlobalEnv))
include_seagrass_species <- isTRUE(get("include_seagrass_species", envir = .GlobalEnv))

models_default <- get("fold_sensitivity_models", envir = .GlobalEnv)

# CV parameter blocks
n_folds_block_for_spatial <- NA_integer_
cv_blocksize <- get("cv_blocksize", envir = .GlobalEnv)

out_base <- file.path(project_root, "output", cv_regime_name, "cv_pipeline", "fold_sensitivity_check")
dir.create(out_base, recursive = TRUE, showWarnings = FALSE)

cat("Fold sensitivity check\n")
cat("  cv_regime_name:", cv_regime_name, "\n")
cat("  target_var:", target_var, "\n")
cat("  log_response:", log_response, "\n")
cat("  exclude_regions:", ifelse(length(exclude_regions) == 0, "(none)", paste(exclude_regions, collapse = ", ")), "\n")
cat("  seed_list:", paste(seed_list, collapse = ", "), "\n")
cat("  n_folds_list:", paste(n_folds_list, collapse = ", "), "\n")
cat("  cv_types_to_check:", paste(cv_types_to_check, collapse = ", "), "\n")

# ---------------------------------------------------------------------------
# Load and prepare data
# ---------------------------------------------------------------------------
dat <- readr::read_rds(file.path(project_root, "data/all_extracted_new.rds"))

if (length(exclude_regions) > 0L) {
  if (!"region" %in% names(dat)) dat <- assign_region_from_latlon(dat)
  dat <- dat[is.na(dat$region) | !dat$region %in% exclude_regions, ]
}

coords_to_keep <- c("longitude", "latitude")
species_region_cols <- intersect(
  c(if (include_seagrass_species) "seagrass_species" else character(0), "region"),
  names(dat)
)

# Predictors used for pixel-group hashing: everything except coords and identifiers
predictor_vars_full <- setdiff(
  colnames(dat),
  c(
    "latitude", "longitude", "number_id_final_version",
    "seagrass_species",
    "region", target_var
  )
)
predictor_vars_full <- predictor_vars_full[predictor_vars_full %in% colnames(dat)]

complete_dat <- dat %>%
  dplyr::select(
    dplyr::all_of(coords_to_keep),
    dplyr::all_of(target_var),
    dplyr::all_of(predictor_vars_full),
    dplyr::all_of(species_region_cols)
  ) %>%
  dplyr::filter(complete.cases(.))

if (nrow(complete_dat) < 10L) stop("Too few complete-case rows after filtering.")

complete_dat$median_carbon_density <- complete_dat[[target_var]]
predictor_vars_full <- predictor_vars_full[predictor_vars_full %in% colnames(complete_dat)]

cat("  Complete-case rows:", nrow(complete_dat), "\n")
cat("  Full predictor count (hashing + run_cv predictor_vars):", length(predictor_vars_full), "\n")

core_data <- as.data.frame(complete_dat)

# ---------------------------------------------------------------------------
# Load tuned covariate sets and hyperparameters (held fixed during sensitivity check)
# ---------------------------------------------------------------------------
cov_dir <- file.path(project_root, "output", cv_regime_name, "covariate_selection")
config_dir <- file.path(project_root, "output", cv_regime_name, "cv_pipeline")

per_model_vars <- get_per_model_vars(cov_dir, colnames(core_data), use_shap_first = use_shap_per_model)

predictor_vars_by_model <- list()
for (m in models_default) {
  pvars <- load_model_vars(m, per_model_vars, use_shap_first = use_shap_per_model)
  if (!is.null(pvars) && length(pvars) >= 2L) predictor_vars_by_model[[m]] <- pvars
}
models <- intersect(names(predictor_vars_by_model), models_default)

hp_bundle <- build_hyperparams_by_model(
  models = models,
  config_dir = config_dir,
  robust_config_dir = NULL,
  prefer_robust = FALSE,
  include_baseline = TRUE,
  include_legacy = TRUE,
  include_only_with_config = TRUE
)
hyperparams_by_model <- hp_bundle$hyperparams_by_model

models <- intersect(models, names(hyperparams_by_model))
if (length(models) == 0L) stop("No models with both tuned covariate sets and tuned hyperparameters found.")

cat("  Using models:", paste(models, collapse = ", "), "\n")

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
compute_fold_y_stats <- function(fold_indices, y) {
  nf <- max(fold_indices, na.rm = TRUE)
  fold_ids <- seq_len(nf)
  idx_list <- lapply(fold_ids, function(k) which(fold_indices == k))
  n_test <- vapply(idx_list, length, integer(1))
  y_mean <- vapply(idx_list, function(ii) {
    if (length(ii) == 0L) return(NA_real_)
    mean(y[ii], na.rm = TRUE)
  }, numeric(1))
  y_sd <- vapply(idx_list, function(ii) {
    if (length(ii) < 2L) return(NA_real_)
    stats::sd(y[ii], na.rm = TRUE)
  }, numeric(1))
  tibble::tibble(fold = fold_ids, n_test = n_test, y_mean = y_mean, y_sd = y_sd)
}

# ---------------------------------------------------------------------------
# Main loop
# ---------------------------------------------------------------------------
by_fold_rows <- list()

for (n_folds in n_folds_list) {
  for (seed in seed_list) {
    for (cv_type in cv_types_to_check) {
      cv_fold_info <- make_cv_folds(
        dat = core_data,
        covariate_cols = predictor_vars_full,
        n_folds = n_folds,
        cv_type = cv_type,
        cv_blocksize = if (!is.na(n_folds_block_for_spatial)) n_folds_block_for_spatial else cv_blocksize,
        cache_tag = paste0("fold_sens_", cv_type, "_n", n_folds, "_seed", seed),
        exclude_regions = exclude_regions,
        seed = seed
      )

      fold_indices <- cv_fold_info$fold_indices
      method_name <- cv_fold_info$method_name

      fold_y_stats <- compute_fold_y_stats(fold_indices, core_data$median_carbon_density)

      res <- run_cv(
        cv_method_name = paste0(method_name, "_n", n_folds, "_seed", seed),
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
      res$fold_y_sd <- fold_y_stats$y_sd[match(res$fold, fold_y_stats$fold)]
      res$fold_y_mean <- fold_y_stats$y_mean[match(res$fold, fold_y_stats$fold)]
      res$n_folds <- n_folds
      res$fold_seed <- seed
      res$cv_type <- cv_type
      res$method_name <- method_name

      by_fold_rows[[length(by_fold_rows) + 1L]] <- res
    }
  }
}

by_fold <- dplyr::bind_rows(by_fold_rows)
if (nrow(by_fold) == 0L) stop("No results produced; check inputs and tuned configs.")

write.csv(by_fold, file.path(out_base, "fold_sensitivity_by_fold.csv"), row.names = FALSE)

summary_by_seed <- by_fold %>%
  group_by(method_name, cv_type, n_folds, fold_seed, model) %>%
  summarise(
    mean_r2 = mean(r2, na.rm = TRUE),
    sd_r2 = sd(r2, na.rm = TRUE),
    neg_r2_fraction = mean(r2 < 0, na.rm = TRUE),
    mean_rmse = mean(rmse, na.rm = TRUE),
    sd_rmse = sd(rmse, na.rm = TRUE),
    mean_fold_y_sd = mean(fold_y_sd, na.rm = TRUE),
    .groups = "drop"
  )
write.csv(summary_by_seed, file.path(out_base, "fold_sensitivity_summary_by_seed.csv"), row.names = FALSE)

summary_across_seeds <- summary_by_seed %>%
  group_by(method_name, cv_type, n_folds, model) %>%
  summarise(
    mean_mean_rmse = mean(mean_rmse, na.rm = TRUE),
    sd_mean_rmse = sd(mean_rmse, na.rm = TRUE),
    mean_mean_r2 = mean(mean_r2, na.rm = TRUE),
    sd_mean_r2 = sd(mean_r2, na.rm = TRUE),
    mean_neg_r2_fraction = mean(neg_r2_fraction, na.rm = TRUE),
    .groups = "drop"
  )
write.csv(summary_across_seeds, file.path(out_base, "fold_sensitivity_summary_across_seeds.csv"), row.names = FALSE)

# Plots
rmse_plot <- ggplot(summary_across_seeds, aes(x = n_folds, y = mean_mean_rmse, colour = model)) +
  geom_line() +
  geom_point(size = 2.2) +
  facet_wrap(~ method_name, scales = "free_y") +
  theme_minimal(base_size = 12) +
  labs(title = "Fold sensitivity (RMSE): mean across seeds", x = "n_folds", y = "Mean RMSE")
ggsave(file.path(out_base, "fold_sensitivity_rmse_plot.png"), rmse_plot, width = 11, height = 6, dpi = 200)

r2_plot <- ggplot(summary_across_seeds, aes(x = n_folds, y = mean_mean_r2, colour = model)) +
  geom_line() +
  geom_point(size = 2.2) +
  facet_wrap(~ method_name, scales = "free_y") +
  theme_minimal(base_size = 12) +
  labs(title = "Fold sensitivity (R2): mean across seeds", x = "n_folds", y = "Mean R2")
ggsave(file.path(out_base, "fold_sensitivity_r2_plot.png"), r2_plot, width = 11, height = 6, dpi = 200)

cat("\nWrote fold sensitivity outputs to:\n", out_base, "\n", sep = "")

