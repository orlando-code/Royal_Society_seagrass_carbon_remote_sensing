## Compare 5-fold CV performance of a fitted model vs a species-mean baseline.
##
## Goal:
##   Quantify whether environmental covariates improve predictive skill beyond
##   a simple baseline that predicts the training-fold mean for each
##   seagrass species.
##
## Outputs:
##   output/<cv_regime_name>/analysis/
##     - model_comparison_by_fold.csv
##     - model_comparison_summary.csv (mean across fold-wise metrics)
##     - model_comparison_pooled_by_seed.csv (RMSE/R² on all held-out points per seed)
##     - model_comparison_pooled_summary.csv
##     - model_comparison_paired_deltas.csv
##     - model_comparison_paired_pooled_deltas.csv
if (!exists("seagrass_init_repo", mode = "function", inherits = TRUE)) {
  init_path <- file.path("modelling", "R", "init_repo.R")
  if (!file.exists(init_path)) {
    ff <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
    if (!length(ff)) stop("Run from repo root or with: Rscript /path/to/model_comparison.R", call. = FALSE)
    script_path <- normalizePath(sub("^--file=", "", ff[[1]]), winslash = "/", mustWork = FALSE)
    init_path <- normalizePath(file.path(dirname(script_path), "..", "R", "init_repo.R"), winslash = "/", mustWork = FALSE)
  }
  if (!file.exists(init_path)) stop("Missing bootstrap helper: modelling/R/init_repo.R", call. = FALSE)
  sys.source(init_path, envir = .GlobalEnv)
}
project_root <- seagrass_init_repo(
  packages = c("here", "dplyr", "readr", "tidyr", "ggplot2", "patchwork"),
  source_files = c("modelling/pipeline_config.R", "modelling/plots/plot_config.R")
)
assign("plot_model_comparison_disable_autorun_on_source", TRUE, envir = .GlobalEnv)
source(file.path(project_root, "modelling/plots/plot_model_comparison.R"))
if (exists("plot_model_comparison_disable_autorun_on_source", envir = .GlobalEnv, inherits = FALSE)) {
  rm("plot_model_comparison_disable_autorun_on_source", envir = .GlobalEnv)
}

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("fit_gpr", "fit_rf", "fit_xgboost", "fit_gam"))
}

cfg <- get_pipeline_config()
apply_pipeline_defaults(
  cfg,
  c(
    "cv_regime_name", "cv_output_dir", "cv_type", "cv_blocksize", "run_output_dir",
    "target_var", "log_transform_target", "exclude_regions",
    "use_shap_per_model", "include_seagrass_species", "robust_fold_seed_list",
    "eval_fold_seed_list", "n_folds",
    "dpi", "model_list", "robust_pruned_importance_type"
  ), envir = .GlobalEnv)
# Optional strict alignment to a specific robust evaluation run directory.
# Example:
#   model_comparison_eval_run_dir <- "output/pixel_grouped/cv_pipeline/pixel_grouped_evaluation_22x5_seeds_45-59-63-69-77"
# When set, this script:
#   1) infers robust_fold_seed_list from the folder name,
#   2) prefers robust tuning configs for that exact seed set, and
#   3) writes outputs to a run-tagged subdirectory.
model_comparison_eval_run_dir <- get0("model_comparison_eval_run_dir", envir = .GlobalEnv, ifnotfound = NA_character_)
model_comparison_eval_run_dir <- as.character(model_comparison_eval_run_dir)
if (length(model_comparison_eval_run_dir) != 1L || !nzchar(model_comparison_eval_run_dir) || is.na(model_comparison_eval_run_dir)) {
  model_comparison_eval_run_dir <- NA_character_
}

parse_seed_list_from_run_dir <- function(run_dir) {
  b <- basename(normalizePath(run_dir, winslash = "/", mustWork = FALSE))
  # Strip multiseed timestamp suffix: ..._seeds_48-52-53_20260402_120000
  b <- sub("_\\d{8}_\\d{6}$", "", b)
  m <- regexec("seeds_([0-9]+(?:-[0-9]+)*)$", b)
  hit <- regmatches(b, m)[[1]]
  if (length(hit) < 2L) return(integer(0))
  as.integer(strsplit(hit[2], "-", fixed = TRUE)[[1]])
}

target_robust_seeds <- integer(0)
if (!is.na(model_comparison_eval_run_dir)) {
  if (!dir.exists(model_comparison_eval_run_dir)) {
    stop("model_comparison_eval_run_dir does not exist: ", model_comparison_eval_run_dir)
  }
  target_robust_seeds <- parse_seed_list_from_run_dir(model_comparison_eval_run_dir)
  if (length(target_robust_seeds) == 0L || any(!is.finite(target_robust_seeds))) {
    stop("Could not parse robust seed list from model_comparison_eval_run_dir basename: ", basename(model_comparison_eval_run_dir))
  }
  robust_fold_seed_list <- as.integer(target_robust_seeds)
}

n_folds_benchmark <- as.integer(get0("n_folds_benchmark", envir = .GlobalEnv, ifnotfound = n_folds))
fold_seed <- as.integer(get0("fold_seed", envir = .GlobalEnv, ifnotfound = 42L))
eval_seed_list <- as.integer(get0("eval_fold_seed_list", envir = .GlobalEnv, ifnotfound = fold_seed))
eval_seed_list <- unique(eval_seed_list[is.finite(eval_seed_list)])
if (length(eval_seed_list) == 0L) eval_seed_list <- fold_seed
cat("Model vs species-mean baseline benchmark\n")
cat("  models:", paste(model_list, collapse = ", "), "\n")
cat("  cv_type:", cv_type, "\n")
cat("  n_folds:", n_folds, "\n")
cat("  robust_fold_seed_list:", paste(robust_fold_seed_list, collapse = ", "), "\n")
cat("  eval_fold_seed_list:", paste(eval_seed_list, collapse = ", "), "\n")
if (!is.na(model_comparison_eval_run_dir)) {
  cat("  model_comparison_eval_run_dir:", model_comparison_eval_run_dir, "\n")
}

run_tag <- if (!is.na(model_comparison_eval_run_dir)) basename(normalizePath(model_comparison_eval_run_dir, winslash = "/", mustWork = FALSE)) else "current_config"
out_dir <- file.path(project_root, cv_output_dir, "analysis", run_tag)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

dat <- readr::read_rds(file.path(project_root, "data/all_extracted_new.rds"))
if (length(exclude_regions) > 0L) {
  if (!"region" %in% names(dat)) dat <- assign_region_from_latlon(dat)
  dat <- dat[is.na(dat$region) | !dat$region %in% exclude_regions, ]
}

if (!"seagrass_species" %in% names(dat)) {
  stop("Column 'seagrass_species' is required for species-mean baseline.")
}

predictor_vars_full <- setdiff(
  colnames(dat),
  c("latitude", "longitude", "number_id_final_version", "seagrass_species", "region", target_var)
)
predictor_vars_full <- predictor_vars_full[predictor_vars_full %in% colnames(dat)]

core_data <- dat %>%
  dplyr::select(
    longitude, latitude,
    dplyr::all_of(target_var),
    seagrass_species,
    dplyr::all_of(intersect("region", names(dat))),
    dplyr::all_of(predictor_vars_full)
  ) %>%
  dplyr::filter(complete.cases(.)) %>%
  as.data.frame()

if (nrow(core_data) < 10L) stop("Too few complete rows after filtering.")
core_data$median_carbon_density <- core_data[[target_var]]

active_run_output_dir <- if (!is.na(model_comparison_eval_run_dir)) {
  model_comparison_eval_run_dir
} else {
  get0("run_output_dir", envir = .GlobalEnv, ifnotfound = NA_character_)
}
if (is.na(active_run_output_dir) || !nzchar(as.character(active_run_output_dir))) {
  stop("run_output_dir must be set in .GlobalEnv (or set model_comparison_eval_run_dir) for model_comparison.R")
}
active_run_output_dir <- as.character(active_run_output_dir)
if (!dir.exists(active_run_output_dir)) {
  stop("Run output directory does not exist: ", active_run_output_dir)
}

config_dir <- file.path(project_root, cv_output_dir, "cv_pipeline")
robust_tuning_dir <- file.path(active_run_output_dir, "cv_pipeline", "robust_pixel_grouped_tuning")
if (!dir.exists(robust_tuning_dir)) {
  stop("Required robust tuning directory not found under run output dir: ", robust_tuning_dir)
}
robust_pruned_importance_type <- match.arg(
  get("robust_pruned_importance_type", envir = .GlobalEnv),
  choices = c("perm", "shap")
)
robust_vars <- load_run_scoped_robust_predictor_vars(
  run_output_dir = active_run_output_dir,
  robust_fold_seed_list = robust_fold_seed_list,
  robust_pruned_importance_type = robust_pruned_importance_type,
  valid_cols = colnames(core_data),
  model_list = model_list
)
cat("  robust_tuning_dir:", robust_tuning_dir, "\n")

predictor_vars_by_model <- list()
hyperparams_by_model <- list()
for (m in model_list) {
  pvars <- robust_vars$predictor_vars_by_model[[m]]
  if (is.null(pvars) || length(pvars) < 2L) next
  predictor_vars_by_model[[m]] <- pvars
  cfg_model <- load_best_model_config(
    model_name = m,
    config_dir = config_dir,
    robust_config_dir = robust_tuning_dir,
    prefer_robust = TRUE,
    include_baseline = FALSE,
  )
  hyperparams_by_model[[m]] <- model_hyperparams_from_config(m, cfg_model)
}
predictor_vars_by_model <- predictor_vars_by_model[intersect(names(predictor_vars_by_model), names(hyperparams_by_model))]
hyperparams_by_model <- hyperparams_by_model[intersect(names(hyperparams_by_model), names(predictor_vars_by_model))]
if (length(predictor_vars_by_model) == 0L) {
  stop("No models with both robust pruned covariates and robust tuning configs under run_output_dir.")
}

fit_model_once <- function(model_name, prep, hp = NULL) {
  switch(model_name,
    GPR = get("fit_gpr", mode = "function")(prep$train, prep$predictor_vars, test_data = prep$test, hyperparams = hp),
    RF  = get("fit_rf", mode = "function")(prep$train, prep$test, prep$predictor_vars, hyperparams = hp),
    XGB = get("fit_xgboost", mode = "function")(prep$train, prep$test, prep$predictor_vars, hyperparams = hp),
    GAM = get("fit_gam", mode = "function")(
      prep$train, prep$test, prep$predictor_vars,
      include_spatial = FALSE,
      k_covariate = if (is.list(hp) && !is.null(hp$k_covariate)) hp$k_covariate else 6L
    ),
    LR = {
      form <- stats::as.formula(
        paste0("median_carbon_density ~ ", paste(prep$predictor_vars, collapse = " + "))
      )
      fit <- stats::lm(form, data = prep$train)
      list(predictions = as.numeric(stats::predict(fit, newdata = prep$test)))
    },
    stop("Unsupported model_name: ", model_name)
  )
}

pooled_key <- function(seed, model_col, approach) {
  paste(as.integer(seed), as.character(model_col), as.character(approach), sep = "\t")
}

fold_rows <- list()
pooled_accum <- list()
for (seed in eval_seed_list) {
  cv_fold_info <- make_cv_folds(
    dat = core_data,
    covariate_cols = predictor_vars_full,
    n_folds = n_folds_benchmark,
    cv_type = cv_type,
    cv_blocksize = cv_blocksize,
    cache_tag = paste0("model_comparison_models_", cv_type, "_n", n_folds_benchmark, "_seed", as.integer(seed)),
    exclude_regions = exclude_regions,
    seed = as.integer(seed)
  )
  fold_indices <- cv_fold_info$fold_indices
  n_folds_seed <- max(fold_indices, na.rm = TRUE)

  for (m in model_list) {
    predictor_vars <- predictor_vars_by_model[[m]]
    hyperparams <- hyperparams_by_model[[m]]
    for (k in seq_len(n_folds_seed)) {
      train_raw <- as.data.frame(core_data[fold_indices != k, , drop = FALSE])
      test_raw <- as.data.frame(core_data[fold_indices == k, , drop = FALSE])
      if (nrow(train_raw) < 2L || nrow(test_raw) < 1L) next

      observed_orig <- test_raw$median_carbon_density
      train_fit <- train_raw
      test_fit <- test_raw
      if (log_transform_target) {
        train_fit <- transform_response(train_fit, "median_carbon_density", log = TRUE)
        test_fit <- transform_response(test_fit, "median_carbon_density", log = TRUE)
      }

      prep <- prepare_data_for_model(m, train_fit, test_fit, predictor_vars)
      fit <- tryCatch(fit_model_once(m, prep, hp = hyperparams), error = function(e) NULL)
      pred_model <- if (!is.null(fit)) fit$predictions else rep(NA_real_, nrow(test_raw))
      if (log_transform_target) pred_model <- inverse_response_transform(pred_model, log = TRUE)

      keep <- test_rows_with_factors_in_train(
        train_raw, test_raw, factor_vars = intersect(predictor_vars, c("seagrass_species", "region"))
      )
      obs_sub <- observed_orig[keep]
      pred_model_sub <- pred_model[keep]

      species_means <- tapply(train_raw$median_carbon_density, train_raw$seagrass_species, mean, na.rm = TRUE)
      global_mean <- mean(train_raw$median_carbon_density, na.rm = TRUE)
      pred_baseline <- as.numeric(species_means[as.character(test_raw$seagrass_species)])
      pred_baseline[!is.finite(pred_baseline)] <- global_mean
      pred_baseline_sub <- pred_baseline[keep]

      model_metrics <- calculate_metrics(obs_sub, pred_model_sub)
      baseline_metrics <- calculate_metrics(obs_sub, pred_baseline_sub)

      km <- pooled_key(seed, m, m)
      z_m <- pooled_accum[[km]]
      if (is.null(z_m)) z_m <- list(obs = numeric(), pred = numeric())
      z_m$obs <- c(z_m$obs, obs_sub)
      z_m$pred <- c(z_m$pred, pred_model_sub)
      pooled_accum[[km]] <- z_m

      kb <- pooled_key(seed, m, "species_mean_baseline")
      z_b <- pooled_accum[[kb]]
      if (is.null(z_b)) z_b <- list(obs = numeric(), pred = numeric())
      z_b$obs <- c(z_b$obs, obs_sub)
      z_b$pred <- c(z_b$pred, pred_baseline_sub)
      pooled_accum[[kb]] <- z_b

      fold_rows[[length(fold_rows) + 1L]] <- dplyr::bind_rows(
        data.frame(
          seed = as.integer(seed), model = m, fold = k, approach = m,
          n_test_raw = nrow(test_raw), n_eval = sum(keep),
          model_metrics, stringsAsFactors = FALSE
        ),
        data.frame(
          seed = as.integer(seed), model = m, fold = k, approach = "species_mean_baseline",
          n_test_raw = nrow(test_raw), n_eval = sum(keep),
          baseline_metrics, stringsAsFactors = FALSE
        )
      )
    }
  }
}

by_fold <- dplyr::bind_rows(fold_rows)
if (nrow(by_fold) == 0L) stop("No valid fold results produced.")

pooled_by_seed <- dplyr::bind_rows(lapply(names(pooled_accum), function(nm) {
  parts <- strsplit(nm, "\t", fixed = TRUE)[[1L]]
  seed_i <- as.integer(parts[[1L]])
  model_i <- parts[[2L]]
  approach_i <- parts[[3L]]
  z <- pooled_accum[[nm]]
  met <- calculate_metrics(z$obs, z$pred)
  dplyr::bind_cols(
    data.frame(
      seed = seed_i,
      model = model_i,
      approach = approach_i,
      stringsAsFactors = FALSE
    ),
    met
  )
}))
if (nrow(pooled_by_seed) == 0L) stop("No pooled results produced.")

pooled_summary_tbl <- pooled_by_seed %>%
  dplyr::group_by(model, approach) %>%
  dplyr::summarise(
    seeds = dplyr::n(),
    mean_r2 = mean(r2, na.rm = TRUE),
    sd_r2 = sd(r2, na.rm = TRUE),
    mean_rmse = mean(rmse, na.rm = TRUE),
    sd_rmse = sd(rmse, na.rm = TRUE),
    mean_mae = mean(mae, na.rm = TRUE),
    sd_mae = sd(mae, na.rm = TRUE),
    mean_bias = mean(bias, na.rm = TRUE),
    sd_bias = sd(bias, na.rm = TRUE),
    mean_n_eval = mean(n_eval, na.rm = TRUE),
    num_na = sum(is.na(r2)),
    .groups = "drop"
  )

summary_tbl <- by_fold %>%
  dplyr::group_by(model, approach) %>%
  dplyr::summarise(
    folds = dplyr::n(),
    mean_r2 = mean(r2, na.rm = TRUE),
    sd_r2 = sd(r2, na.rm = TRUE),
    mean_rmse = mean(rmse, na.rm = TRUE),
    sd_rmse = sd(rmse, na.rm = TRUE),
    mean_mae = mean(mae, na.rm = TRUE),
    sd_mae = sd(mae, na.rm = TRUE),
    mean_bias = mean(bias, na.rm = TRUE),
    sd_bias = sd(bias, na.rm = TRUE),
    mean_n_eval = mean(n_eval, na.rm = TRUE),
    num_na = sum(is.na(r2)),
    .groups = "drop"
  )

paired_model <- by_fold %>%
  dplyr::filter(approach != "species_mean_baseline") %>%
  dplyr::select(seed, model, fold, rmse_model = rmse, r2_model = r2, mae_model = mae, bias_model = bias)

paired_baseline <- by_fold %>%
  dplyr::filter(approach == "species_mean_baseline") %>%
  dplyr::select(seed, model, fold, rmse_baseline = rmse, r2_baseline = r2, mae_baseline = mae, bias_baseline = bias)

paired_delta <- paired_model %>%
  dplyr::left_join(paired_baseline, by = c("seed", "model", "fold")) %>%
  dplyr::mutate(
    delta_rmse = rmse_model - rmse_baseline,
    delta_r2 = r2_model - r2_baseline,
    delta_mae = mae_model - mae_baseline,
    delta_bias = bias_model - bias_baseline,
    percentage_improvement_r2 = (r2_model - r2_baseline) / r2_baseline * 100,
    percentage_improvement_rmse = (rmse_model - rmse_baseline) / rmse_baseline * 100,
    percentage_improvement_mae = (mae_model - mae_baseline) / mae_baseline * 100,
    percentage_improvement_bias = (bias_model - bias_baseline) / bias_baseline * 100
  )

paired_delta_summary <- paired_delta %>%
  dplyr::group_by(model) %>%
  dplyr::summarise(
    mean_delta_rmse = mean(delta_rmse, na.rm = TRUE),
    mean_delta_r2 = mean(delta_r2, na.rm = TRUE),
    mean_delta_mae = mean(delta_mae, na.rm = TRUE),
    mean_delta_bias = mean(delta_bias, na.rm = TRUE),
    mean_percentage_improvement_r2 = mean(percentage_improvement_r2, na.rm = TRUE),
    mean_percentage_improvement_rmse = mean(percentage_improvement_rmse, na.rm = TRUE),
    mean_percentage_improvement_mae = mean(percentage_improvement_mae, na.rm = TRUE),
    mean_percentage_improvement_bias = mean(percentage_improvement_bias, na.rm = TRUE),
    .groups = "drop"
  )

paired_pooled_model <- pooled_by_seed %>%
  dplyr::filter(.data$approach != "species_mean_baseline") %>%
  dplyr::select(seed, model, rmse_model = rmse, r2_model = r2, mae_model = mae, bias_model = bias)

paired_pooled_baseline <- pooled_by_seed %>%
  dplyr::filter(.data$approach == "species_mean_baseline") %>%
  dplyr::select(seed, model, rmse_baseline = rmse, r2_baseline = r2, mae_baseline = mae, bias_baseline = bias)

paired_pooled_delta <- paired_pooled_model %>%
  dplyr::left_join(paired_pooled_baseline, by = c("seed", "model")) %>%
  dplyr::mutate(
    delta_rmse = rmse_model - rmse_baseline,
    delta_r2 = r2_model - r2_baseline,
    delta_mae = mae_model - mae_baseline,
    delta_bias = bias_model - bias_baseline,
    percentage_improvement_r2 = (r2_model - r2_baseline) / r2_baseline * 100,
    percentage_improvement_rmse = (rmse_model - rmse_baseline) / rmse_baseline * 100,
    percentage_improvement_mae = (mae_model - mae_baseline) / mae_baseline * 100,
    percentage_improvement_bias = (bias_model - bias_baseline) / bias_baseline * 100
  )

paired_pooled_delta_summary <- paired_pooled_delta %>%
  dplyr::group_by(model) %>%
  dplyr::summarise(
    mean_delta_rmse = mean(delta_rmse, na.rm = TRUE),
    mean_delta_r2 = mean(delta_r2, na.rm = TRUE),
    mean_delta_mae = mean(delta_mae, na.rm = TRUE),
    mean_delta_bias = mean(delta_bias, na.rm = TRUE),
    mean_percentage_improvement_r2 = mean(percentage_improvement_r2, na.rm = TRUE),
    mean_percentage_improvement_rmse = mean(percentage_improvement_rmse, na.rm = TRUE),
    mean_percentage_improvement_mae = mean(percentage_improvement_mae, na.rm = TRUE),
    mean_percentage_improvement_bias = mean(percentage_improvement_bias, na.rm = TRUE),
    .groups = "drop"
  )

plot_model_comparison_outputs(
  summary_tbl = summary_tbl,
  pooled_summary_tbl = pooled_summary_tbl,
  by_fold = by_fold,
  model_list = model_list,
  out_dir = out_dir,
  dpi = dpi
)

readr::write_csv(by_fold, file.path(out_dir, "model_comparison_by_fold.csv"))
readr::write_csv(summary_tbl, file.path(out_dir, "model_comparison_summary.csv"))
readr::write_csv(pooled_by_seed, file.path(out_dir, "model_comparison_pooled_by_seed.csv"))
readr::write_csv(pooled_summary_tbl, file.path(out_dir, "model_comparison_pooled_summary.csv"))
readr::write_csv(paired_delta, file.path(out_dir, "model_comparison_paired_deltas.csv"))
readr::write_csv(paired_delta_summary, file.path(out_dir, "model_comparison_paired_deltas_summary.csv"))
readr::write_csv(paired_pooled_delta, file.path(out_dir, "model_comparison_paired_pooled_deltas.csv"))
readr::write_csv(paired_pooled_delta_summary, file.path(out_dir, "model_comparison_paired_pooled_deltas_summary.csv"))

cat("\nSaved outputs to:", out_dir, "\n")
cat("Mean RMSE / R2 by model/approach (fold-wise means over all fold×seed rows):\n")
print(summary_tbl[, c("model", "approach", "mean_rmse", "sd_rmse", "mean_r2", "sd_r2")])
cat("\nPooled metrics by model/approach (mean ± sd across evaluation seeds):\n")
print(pooled_summary_tbl[, c("model", "approach", "mean_rmse", "sd_rmse", "mean_r2", "sd_r2")])
cat("\nFold-wise deltas by model (model - species baseline):\n")
print(paired_delta_summary)
cat("\nPooled deltas by model (model - species baseline, per seed then averaged):\n")
print(paired_pooled_delta_summary)


