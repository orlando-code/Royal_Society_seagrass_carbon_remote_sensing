## Fast sensitivity suite for pooled vs mean-fold R2.
##
## What it does:
##   1) Varies fold count (k) and fold seed, with fixed tuned model settings.
##   2) Quantifies training-size effects using n_train_raw from each fold result.
##   3) Quantifies fold-composition sensitivity via SS_total / y-spread by fold.
##
## Design goal: fast diagnostics only.
##   - Reuses tuned covariate sets and tuned hyperparameters.
##   - No re-tuning, no re-selection.
##
## Outputs (under <run_output_dir>/sensitivity_suite/, where <run_output_dir> is
##   output/pixel_grouped_<robust_seed_subset>/ from run_multiseed_pixel_grouped.R):
##   - sensitivity_by_fold.csv
##   - sensitivity_pooled_by_seed.csv
##   - sensitivity_summary.csv
##   - sensitivity_fold_composition.csv
##   - sensitivity_train_size_effect.csv
##   - sensitivity_ss_total_r2_effect.csv
##   - sensitivity_seed_count_convergence.csv
##   - sensitivity_seed_count_plateau.csv
if (!exists("seagrass_init_repo", mode = "function", inherits = TRUE)) {
  init_path <- file.path("modelling", "R", "init_repo.R")
  if (!file.exists(init_path)) {
    ff <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
    if (!length(ff)) stop("Run from repo root or with: Rscript /path/to/sensitivity_suite.R", call. = FALSE)
    script_path <- normalizePath(sub("^--file=", "", ff[[1]]), winslash = "/", mustWork = FALSE)
    init_path <- normalizePath(file.path(dirname(script_path), "..", "R", "init_repo.R"), winslash = "/", mustWork = FALSE)
  }
  if (!file.exists(init_path)) stop("Missing bootstrap helper: modelling/R/init_repo.R", call. = FALSE)
  sys.source(init_path, envir = .GlobalEnv)
}
project_root <- seagrass_init_repo(
  packages = c("here", "dplyr", "readr", "tibble"),
  source_files = c("modelling/pipeline_config.R"),
  include_helpers = TRUE,
  require_core_inputs = TRUE,
  check_renv = TRUE
)

cfg <- get_pipeline_config()
apply_pipeline_defaults(
  cfg,
  c(
    "cv_regime_name", "target_var", "log_transform_target", "exclude_regions",
    "cv_type_label", "robust_fold_seed_list",
    "cv_types_to_check", "n_folds_list", "fold_seed_list", "include_loo", "max_loo_groups",
    "seed_plateau_tolerance", "seed_plateau_min_n", "seed_plateau_stable_steps",
    "use_shap_per_model", "model_list", "include_seagrass_species",
    "do_tuning_seed_sweep", "tuning_seed_sweep_counts", "tuning_seed_pool",
    "eval_fold_seed_list", "tuning_seed_sampling", "do_tuning_seed_sweep_refined_tuning",
    "robust_pruned_importance_type"
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

cv_types_to_check <- get("cv_types_to_check", envir = .GlobalEnv)
allowed_cv_types <- c("random", "location_grouped", "pixel_grouped", "spatial")
cv_types_to_check <- intersect(as.character(cv_types_to_check), allowed_cv_types)
if (length(cv_types_to_check) == 0L) {
  cv_types_to_check <- "pixel_grouped"
  warning("cv_types_to_check had no recognised entries; using 'pixel_grouped'.", call. = FALSE)
}
n_folds_list <- as.integer(get("n_folds_list", envir = .GlobalEnv))
fold_seed_list <- as.integer(get("fold_seed_list", envir = .GlobalEnv))
include_loo <- isTRUE(get("include_loo", envir = .GlobalEnv))
max_loo_groups <- as.integer(get("max_loo_groups", envir = .GlobalEnv))
seed_plateau_tolerance <- as.numeric(get("seed_plateau_tolerance", envir = .GlobalEnv))
seed_plateau_min_n <- as.integer(get("seed_plateau_min_n", envir = .GlobalEnv))
seed_plateau_stable_steps <- as.integer(get("seed_plateau_stable_steps", envir = .GlobalEnv))

use_shap_per_model <- isTRUE(get("use_shap_per_model", envir = .GlobalEnv))
models_default <- get("model_list", envir = .GlobalEnv)
include_seagrass_species <- isTRUE(get("include_seagrass_species", envir = .GlobalEnv))

# Optional expensive block:
# re-run robust tune/prune/eval for different numbers of tuning seeds.
do_tuning_seed_sweep <- isTRUE(get("do_tuning_seed_sweep", envir = .GlobalEnv))
tuning_seed_sweep_counts <- as.integer(get("tuning_seed_sweep_counts", envir = .GlobalEnv))
tuning_seed_pool <- as.integer(get("tuning_seed_pool", envir = .GlobalEnv))
eval_fold_seed_list <- as.integer(get("eval_fold_seed_list", envir = .GlobalEnv))
tuning_seed_sampling <- match.arg(
  get("tuning_seed_sampling", envir = .GlobalEnv),
  choices = c("prefix", "random")
)
do_tuning_seed_sweep_refined_tuning <- isTRUE(get("do_tuning_seed_sweep_refined_tuning", envir = .GlobalEnv))
robust_pruned_importance_type <- match.arg(
  get("robust_pruned_importance_type", envir = .GlobalEnv),
  choices = c("perm", "shap")
)

cv_type_label <- get("cv_type_label", envir = .GlobalEnv)
robust_fold_seed_list <- as.integer(get("robust_fold_seed_list", envir = .GlobalEnv))

resolve_run_output_dir <- function() {
  d <- get0("run_output_dir", envir = .GlobalEnv, ifnotfound = NA_character_)
  if (is.na(d) || !nzchar(as.character(d))) {
    stop("run_output_dir must be set in .GlobalEnv for sensitivity_suite.R")
  }
  as.character(d)
}

run_output_dir <- resolve_run_output_dir()
assign("run_output_dir", run_output_dir, envir = .GlobalEnv)
out_dir <- file.path(run_output_dir, "sensitivity_suite")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("R2 sensitivity suite (fast)\n")
cat("  cv_regime_name:", cv_regime_name, "\n")
cat("  target_var:", target_var, "\n")
cat("  log_response:", log_response, "\n")
cat("  cv_types_to_check:", paste(cv_types_to_check, collapse = ", "), "\n")
cat("  n_folds_list:", paste(n_folds_list, collapse = ", "), "\n")
cat("  fold_seed_list:", paste(fold_seed_list, collapse = ", "), "\n")
cat("  run_output_dir:", run_output_dir, "\n")
cat("  sensitivity output dir:", out_dir, "\n")
cat("  seed plateau tolerance:", seed_plateau_tolerance, "\n")

prepare_core_data <- function(dat, target_var, exclude_regions, include_seagrass_species = TRUE) {
  if (length(exclude_regions) > 0L) {
    if (!"region" %in% names(dat)) dat <- assign_region_from_latlon(dat)
    dat <- dat[is.na(dat$region) | !dat$region %in% exclude_regions, , drop = FALSE]
  }
  species_region_cols <- intersect(
    c(if (isTRUE(include_seagrass_species)) "seagrass_species" else character(0), "region"),
    names(dat)
  )
  predictor_vars_full <- setdiff(
    colnames(dat),
    c(
      "latitude", "longitude", "number_id_final_version",
      "seagrass_species",
      "region", target_var
    )
  )
  predictor_vars_full <- predictor_vars_full[predictor_vars_full %in% colnames(dat)]
  core_data <- dat %>%
    dplyr::select(
      dplyr::all_of(c("longitude", "latitude")),
      dplyr::all_of(target_var),
      dplyr::all_of(predictor_vars_full),
      dplyr::all_of(species_region_cols)
    ) %>%
    dplyr::filter(complete.cases(.)) %>%
    as.data.frame()
  core_data$median_carbon_density <- core_data[[target_var]]
  list(core_data = core_data, predictor_vars_full = predictor_vars_full)
}

load_fixed_model_setup <- function(project_root, run_output_dir, cv_regime_name, core_cols, models_default,
                                   use_shap_per_model, robust_fold_seed_list,
                                   robust_pruned_importance_type) {
  cfg_dir <- file.path(run_output_dir, "cv_pipeline")
  robust_tuning_dir <- file.path(cfg_dir, "robust_pixel_grouped_tuning")
  if (!dir.exists(robust_tuning_dir)) {
    stop("Required robust tuning directory not found under run_output_dir: ", robust_tuning_dir)
  }
  robust_vars <- load_run_scoped_robust_predictor_vars(
    run_output_dir = run_output_dir,
    robust_fold_seed_list = robust_fold_seed_list,
    robust_pruned_importance_type = robust_pruned_importance_type,
    valid_cols = core_cols,
    model_list = models_default
  )
  predictor_vars_by_model <- robust_vars$predictor_vars_by_model
  models <- intersect(models_default, names(predictor_vars_by_model))
  if (length(models) == 0L) {
    stop(
      "No models with valid pruned covariates found. ",
      "Expected run-scoped robust covariate file: ", robust_vars$robust_pruned_csv
    )
  }
  cat("  Using robust tuning config directory:\n    ", robust_tuning_dir, "\n", sep = "")
  hp_bundle <- build_hyperparams_by_model(
    models = models,
    config_dir = file.path(project_root, "output", cv_regime_name, "cv_pipeline"),
    robust_config_dir = robust_tuning_dir,
    prefer_robust = TRUE,
    include_baseline = TRUE,
    include_only_with_config = TRUE
  )
  hyperparams_by_model <- hp_bundle$hyperparams_by_model
  models <- intersect(models, names(hyperparams_by_model))
  if (length(models) == 0L) stop("No models with both vars + fixed tuned hyperparameters found.")

  list(models = models, predictor_vars_by_model = predictor_vars_by_model, hyperparams_by_model = hyperparams_by_model)
}

fold_composition_stats <- function(fold_indices, y) {
  tibble::tibble(row_id = seq_along(fold_indices), fold = fold_indices, observed = y) %>%
    dplyr::group_by(fold) %>%
    dplyr::summarise(
      n_test = dplyr::n(),
      y_mean = mean(observed, na.rm = TRUE),
      y_sd = stats::sd(observed, na.rm = TRUE),
      y_iqr = stats::IQR(observed, na.rm = TRUE),
      ss_total_fold = sum((observed - mean(observed, na.rm = TRUE))^2, na.rm = TRUE),
      sum_y_fold = sum(observed, na.rm = TRUE),
      sum_y2_fold = sum(observed^2, na.rm = TRUE),
      .groups = "drop"
    )
}

compute_pooled_ss <- function(df) {
  # df is seed-level fold metrics already grouped to one model/setting.
  n_total <- sum(df$n_eval, na.rm = TRUE)
  sum_y <- sum(df$sum_y, na.rm = TRUE)
  sum_y2 <- sum(df$sum_y2, na.rm = TRUE)
  ss_res_total <- sum(df$ss_residual, na.rm = TRUE)
  ss_tot_global <- sum_y2 - (sum_y^2 / n_total)
  pooled_r2_reconstructed <- if (is.finite(ss_tot_global) && ss_tot_global > 0) {
    1 - ss_res_total / ss_tot_global
  } else {
    NA_real_
  }
  tibble::tibble(
    n_eval_total = n_total,
    ss_residual_total = ss_res_total,
    ss_total_global = ss_tot_global,
    pooled_r2_reconstructed = pooled_r2_reconstructed
  )
}

safe_quantile <- function(x, probs) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(rep(NA_real_, length(probs)))
  as.numeric(stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE, type = 7))
}

compute_fold_environment_stats <- function(dat, fold_indices, predictor_vars_full) {
  pv <- intersect(
    predictor_vars_full,
    names(dat)[vapply(dat, is.numeric, logical(1))]
  )
  if (length(pv) == 0L) {
    return(tibble::tibble(
      fold = sort(unique(fold_indices)),
      env_n_predictors = 0L,
      env_mean_pairwise_dist = NA_real_,
      env_median_pairwise_dist = NA_real_,
      env_mean_abs_z = NA_real_,
      env_median_abs_z = NA_real_,
      env_mean_predictor_sd = NA_real_,
      env_mean_predictor_iqr = NA_real_
    ))
  }

  X <- as.matrix(dat[, pv, drop = FALSE])
  Xs <- scale(X)
  Xs[!is.finite(Xs)] <- NA_real_

  fold_ids <- sort(unique(fold_indices))
  out <- lapply(fold_ids, function(f) {
    idx <- which(fold_indices == f)
    if (length(idx) == 0L) {
      return(data.frame(
        fold = f,
        env_n_predictors = length(pv),
        env_mean_pairwise_dist = NA_real_,
        env_median_pairwise_dist = NA_real_,
        env_mean_abs_z = NA_real_,
        env_median_abs_z = NA_real_,
        env_mean_predictor_sd = NA_real_,
        env_mean_predictor_iqr = NA_real_
      ))
    }

    Xf <- Xs[idx, , drop = FALSE]
    row_mean_abs_z <- rowMeans(abs(Xf), na.rm = TRUE)
    env_mean_abs_z <- mean(row_mean_abs_z, na.rm = TRUE)
    env_median_abs_z <- stats::median(row_mean_abs_z, na.rm = TRUE)

    pred_sd <- apply(Xf, 2, stats::sd, na.rm = TRUE)
    pred_iqr <- apply(Xf, 2, stats::IQR, na.rm = TRUE)
    env_mean_predictor_sd <- mean(pred_sd, na.rm = TRUE)
    env_mean_predictor_iqr <- mean(pred_iqr, na.rm = TRUE)

    complete_rows <- stats::complete.cases(Xf)
    env_mean_pairwise_dist <- NA_real_
    env_median_pairwise_dist <- NA_real_
    if (sum(complete_rows) >= 2L) {
      d <- as.numeric(stats::dist(Xf[complete_rows, , drop = FALSE]))
      d <- d[is.finite(d)]
      if (length(d) > 0L) {
        env_mean_pairwise_dist <- mean(d, na.rm = TRUE)
        env_median_pairwise_dist <- stats::median(d, na.rm = TRUE)
      }
    }

    data.frame(
      fold = f,
      env_n_predictors = length(pv),
      env_mean_pairwise_dist = env_mean_pairwise_dist,
      env_median_pairwise_dist = env_median_pairwise_dist,
      env_mean_abs_z = env_mean_abs_z,
      env_median_abs_z = env_median_abs_z,
      env_mean_predictor_sd = env_mean_predictor_sd,
      env_mean_predictor_iqr = env_mean_predictor_iqr
    )
  })
  dplyr::bind_rows(out)
}

compute_seed_convergence <- function(vals, seeds, metric_name, tol = 0.02, min_n = 3L, stable_steps = 2L) {
  ord <- order(seeds)
  vals <- as.numeric(vals[ord])
  seeds <- as.integer(seeds[ord])
  keep <- is.finite(vals) & is.finite(seeds)
  vals <- vals[keep]
  seeds <- seeds[keep]
  n <- length(vals)
  if (n == 0L) {
    return(list(
      series = tibble::tibble(
        metric = character(), n_seeds = integer(), fold_seed_last = integer(),
        cumulative_mean = numeric(), cumulative_sd = numeric(),
        cumulative_se = numeric(), delta_from_prev = numeric()
      ),
      plateau = tibble::tibble(
        metric = metric_name, n_available_seeds = 0L, tolerance = tol,
        plateau_n_seeds = NA_integer_, plateau_fold_seed = NA_integer_,
        plateau_value = NA_real_, plateau_reached = FALSE
      )
    ))
  }

  cum_mean <- vapply(seq_len(n), function(i) mean(vals[seq_len(i)], na.rm = TRUE), numeric(1))
  cum_sd <- vapply(seq_len(n), function(i) stats::sd(vals[seq_len(i)], na.rm = TRUE), numeric(1))
  cum_se <- cum_sd / sqrt(seq_len(n))
  delta <- c(NA_real_, abs(diff(cum_mean)))

  plateau_idx <- NA_integer_
  if (n >= max(min_n, stable_steps + 1L)) {
    for (i in seq(max(min_n, stable_steps + 1L), n)) {
      idx <- (i - stable_steps + 1L):i
      if (all(is.finite(delta[idx])) && all(delta[idx] <= tol)) {
        plateau_idx <- i
        break
      }
    }
  }

  list(
    series = tibble::tibble(
      metric = metric_name,
      n_seeds = seq_len(n),
      fold_seed_last = seeds,
      cumulative_mean = cum_mean,
      cumulative_sd = cum_sd,
      cumulative_se = cum_se,
      delta_from_prev = delta
    ),
    plateau = tibble::tibble(
      metric = metric_name,
      n_available_seeds = n,
      tolerance = tol,
      plateau_n_seeds = ifelse(is.na(plateau_idx), NA_integer_, plateau_idx),
      plateau_fold_seed = ifelse(is.na(plateau_idx), NA_integer_, seeds[plateau_idx]),
      plateau_value = ifelse(is.na(plateau_idx), NA_real_, cum_mean[plateau_idx]),
      plateau_reached = !is.na(plateau_idx)
    )
  )
}

fit_lm_slope <- function(df, y_col, x_col) {
  ok <- is.finite(df[[y_col]]) & is.finite(df[[x_col]])
  d <- df[ok, , drop = FALSE]
  n <- nrow(d)
  if (n < 3L) {
    return(tibble::tibble(
      n_lm = n,
      slope = NA_real_,
      slope_se = NA_real_,
      slope_t = NA_real_,
      slope_p = NA_real_,
      slope_ci_low = NA_real_,
      slope_ci_high = NA_real_,
      intercept = NA_real_
    ))
  }
  fit <- stats::lm(stats::as.formula(paste0(y_col, " ~ ", x_col)), data = d)
  sm <- summary(fit)$coefficients
  row <- sm[x_col, , drop = FALSE]
  ci <- tryCatch(stats::confint(fit, parm = x_col, level = 0.95), error = function(e) c(NA_real_, NA_real_))
  tibble::tibble(
    n_lm = n,
    slope = unname(row[1, "Estimate"]),
    slope_se = unname(row[1, "Std. Error"]),
    slope_t = unname(row[1, "t value"]),
    slope_p = unname(row[1, "Pr(>|t|)"]),
    slope_ci_low = unname(ci[1]),
    slope_ci_high = unname(ci[2]),
    intercept = unname(sm["(Intercept)", "Estimate"])
  )
}

dat <- readr::read_rds(file.path(project_root, "data/all_extracted_new.rds"))
prep <- prepare_core_data(dat, target_var, exclude_regions, include_seagrass_species = include_seagrass_species)
core_data <- prep$core_data
predictor_vars_full <- prep$predictor_vars_full

if (nrow(core_data) < 10L) stop("Too few complete-case rows after filtering.")
cat("  complete-case rows:", nrow(core_data), "\n")

setup <- load_fixed_model_setup(
  project_root = project_root,
  run_output_dir = run_output_dir,
  cv_regime_name = cv_regime_name,
  core_cols = colnames(core_data),
  models_default = models_default,
  use_shap_per_model = use_shap_per_model,
  robust_fold_seed_list = robust_fold_seed_list,
  robust_pruned_importance_type = robust_pruned_importance_type
)
models <- setup$models
predictor_vars_by_model <- setup$predictor_vars_by_model
hyperparams_by_model <- setup$hyperparams_by_model
cat("  models:", paste(models, collapse = ", "), "\n")

if (include_loo) {
  group_estimate <- switch(
    cv_types_to_check[[1]],
    random = nrow(core_data),
    location_grouped = length(unique(paste(core_data$longitude, core_data$latitude))),
    pixel_grouped = nrow(unique(core_data[, predictor_vars_full, drop = FALSE])),
    nrow(core_data)
  )
  if (is.finite(group_estimate) && group_estimate <= max_loo_groups) {
    n_folds_list <- unique(c(n_folds_list, as.integer(group_estimate)))
    cat("  include_loo: TRUE (n_folds =", as.integer(group_estimate), ")\n")
  } else {
    cat("  include_loo skipped (estimated groups =", as.integer(group_estimate),
        "> max_loo_groups =", max_loo_groups, ")\n")
  }
}

fold_rows <- list()
pooled_rows <- list()
comp_rows <- list()
env_rows <- list()

# Loop indices must not reuse names `cv_type` / `n_folds` / `seed`: this file is
# typically source()d into .GlobalEnv, and `for (x in y)` assigns globals.
for (cv_type_loop in cv_types_to_check) {
  for (n_folds_k in n_folds_list) {
    for (seed_k in fold_seed_list) {
      cv_fold_info <- make_cv_folds(
        dat = core_data,
        covariate_cols = predictor_vars_full,
        n_folds = n_folds_k,
        cv_type = cv_type_loop,
        cache_tag = paste0("r2_sens_", cv_type_loop, "_k", n_folds_k, "_seed", seed_k),
        exclude_regions = exclude_regions,
        seed = seed_k
      )
      fold_indices <- cv_fold_info$fold_indices
      method_name <- cv_fold_info$method_name

      comp <- fold_composition_stats(fold_indices, core_data$median_carbon_density)
      comp$cv_type <- cv_type_loop
      comp$method_name <- method_name
      comp$n_folds <- n_folds_k
      comp$fold_seed <- seed_k
      comp_rows[[length(comp_rows) + 1L]] <- comp
      env_stats <- compute_fold_environment_stats(core_data, fold_indices, predictor_vars_full)
      env_stats$cv_type <- cv_type_loop
      env_stats$method_name <- method_name
      env_stats$n_folds <- n_folds_k
      env_stats$fold_seed <- seed_k
      env_rows[[length(env_rows) + 1L]] <- env_stats

      cv_name <- paste0(method_name, "_k", n_folds_k, "_seed", seed_k)
      cat("  running:", cv_name, "\n")
      res <- run_cv(
        cv_method_name = cv_name,
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

      res$cv_type <- cv_type_loop
      res$method_name <- method_name
      res$n_folds <- n_folds_k
      res$fold_seed <- seed_k
      fold_rows[[length(fold_rows) + 1L]] <- res

      pooled <- attr(res, "pooled_r2")
      if (!is.null(pooled) && nrow(pooled) > 0L) {
        pooled$cv_type <- cv_type_loop
        pooled$method_name <- method_name
        pooled$n_folds <- n_folds_k
        pooled$fold_seed <- seed_k
        pooled_rows[[length(pooled_rows) + 1L]] <- pooled
      }
    }
  }
}

by_fold <- dplyr::bind_rows(fold_rows)
by_pooled <- dplyr::bind_rows(pooled_rows)
by_comp <- dplyr::bind_rows(comp_rows)
by_env <- dplyr::bind_rows(env_rows)

if (nrow(by_fold) == 0L) stop("No CV results produced.")

by_fold <- by_fold %>%
  dplyr::left_join(
    by_comp %>% dplyr::select(cv_type, method_name, n_folds, fold_seed, fold, n_test, y_mean, y_sd, y_iqr, ss_total_fold),
    by = c("cv_type", "method_name", "n_folds", "fold_seed", "fold")
  ) %>%
  dplyr::left_join(
    by_env %>% dplyr::select(
      cv_type, method_name, n_folds, fold_seed, fold,
      env_n_predictors, env_mean_pairwise_dist, env_median_pairwise_dist,
      env_mean_abs_z, env_median_abs_z, env_mean_predictor_sd, env_mean_predictor_iqr
    ),
    by = c("cv_type", "method_name", "n_folds", "fold_seed", "fold")
  )

seed_level <- by_fold %>%
  dplyr::group_by(cv_type, method_name, n_folds, fold_seed, model) %>%
  dplyr::summarise(
    mean_fold_r2 = mean(r2, na.rm = TRUE),
    sd_fold_r2 = stats::sd(r2, na.rm = TRUE),
    mean_fold_rmse = mean(rmse, na.rm = TRUE),
    neg_fold_r2_fraction = mean(r2 < 0, na.rm = TRUE),
    mean_ss_total_fold = mean(ss_total, na.rm = TRUE),
    min_ss_total_fold = min(ss_total, na.rm = TRUE),
    max_ss_total_fold = max(ss_total, na.rm = TRUE),
    n_eval_sum = sum(n_eval, na.rm = TRUE),
    .groups = "drop"
  )

if (nrow(by_pooled) > 0L) {
  seed_level <- seed_level %>%
    dplyr::left_join(
      by_pooled %>%
        dplyr::rename(pooled_r2_attr = pooled_r2, pooled_rmse_attr = pooled_rmse, pooled_n_predictions = n_predictions),
      by = c("cv_type", "method_name", "n_folds", "fold_seed", "model")
    )
}

reconstructed <- by_fold %>%
  dplyr::group_by(cv_type, method_name, n_folds, fold_seed, model) %>%
  dplyr::group_modify(~ compute_pooled_ss(.x)) %>%
  dplyr::ungroup()

seed_level <- seed_level %>%
  dplyr::left_join(
    reconstructed,
    by = c("cv_type", "method_name", "n_folds", "fold_seed", "model")
  ) %>%
  dplyr::mutate(
    pooled_minus_meanfold_r2 = pooled_r2_reconstructed - mean_fold_r2,
    pooled_attr_minus_reconstructed = pooled_r2_attr - pooled_r2_reconstructed
  )

summary_across_seeds <- seed_level %>%
  dplyr::group_by(cv_type, method_name, n_folds, model) %>%
  dplyr::summarise(
    mean_mean_fold_r2 = mean(mean_fold_r2, na.rm = TRUE),
    sd_mean_fold_r2 = stats::sd(mean_fold_r2, na.rm = TRUE),
    mean_pooled_r2 = mean(pooled_r2_reconstructed, na.rm = TRUE),
    sd_pooled_r2 = stats::sd(pooled_r2_reconstructed, na.rm = TRUE),
    mean_delta_pooled_minus_meanfold = mean(pooled_minus_meanfold_r2, na.rm = TRUE),
    mean_mean_fold_rmse = mean(mean_fold_rmse, na.rm = TRUE),
    mean_pooled_rmse_attr = mean(pooled_rmse_attr, na.rm = TRUE),
    median_mean_fold_r2 = stats::median(mean_fold_r2, na.rm = TRUE),
    p01_mean_fold_r2 = safe_quantile(mean_fold_r2, 0.01),
    p05_mean_fold_r2 = safe_quantile(mean_fold_r2, 0.05),
    p10_mean_fold_r2 = safe_quantile(mean_fold_r2, 0.10),
    p90_mean_fold_r2 = safe_quantile(mean_fold_r2, 0.90),
    p95_mean_fold_r2 = safe_quantile(mean_fold_r2, 0.95),
    p99_mean_fold_r2 = safe_quantile(mean_fold_r2, 0.99),
    median_pooled_r2 = stats::median(pooled_r2_reconstructed, na.rm = TRUE),
    p01_pooled_r2 = safe_quantile(pooled_r2_reconstructed, 0.01),
    p05_pooled_r2 = safe_quantile(pooled_r2_reconstructed, 0.05),
    p10_pooled_r2 = safe_quantile(pooled_r2_reconstructed, 0.10),
    p90_pooled_r2 = safe_quantile(pooled_r2_reconstructed, 0.90),
    p95_pooled_r2 = safe_quantile(pooled_r2_reconstructed, 0.95),
    p99_pooled_r2 = safe_quantile(pooled_r2_reconstructed, 0.99),
    median_mean_fold_rmse = stats::median(mean_fold_rmse, na.rm = TRUE),
    p01_mean_fold_rmse = safe_quantile(mean_fold_rmse, 0.01),
    p05_mean_fold_rmse = safe_quantile(mean_fold_rmse, 0.05),
    p10_mean_fold_rmse = safe_quantile(mean_fold_rmse, 0.10),
    p90_mean_fold_rmse = safe_quantile(mean_fold_rmse, 0.90),
    p95_mean_fold_rmse = safe_quantile(mean_fold_rmse, 0.95),
    p99_mean_fold_rmse = safe_quantile(mean_fold_rmse, 0.99),
    median_pooled_rmse = stats::median(pooled_rmse_attr, na.rm = TRUE),
    p01_pooled_rmse = safe_quantile(pooled_rmse_attr, 0.01),
    p05_pooled_rmse = safe_quantile(pooled_rmse_attr, 0.05),
    p10_pooled_rmse = safe_quantile(pooled_rmse_attr, 0.10),
    p90_pooled_rmse = safe_quantile(pooled_rmse_attr, 0.90),
    p95_pooled_rmse = safe_quantile(pooled_rmse_attr, 0.95),
    p99_pooled_rmse = safe_quantile(pooled_rmse_attr, 0.99),
    mean_neg_fold_r2_fraction = mean(neg_fold_r2_fraction, na.rm = TRUE),
    .groups = "drop"
  )

train_size_effect <- by_fold %>%
  dplyr::group_by(cv_type, method_name, n_folds, model) %>%
  dplyr::summarise(
    n_rows = dplyr::n(),
    mean_n_train = mean(n_train_raw, na.rm = TRUE),
    min_n_train = min(n_train_raw, na.rm = TRUE),
    max_n_train = max(n_train_raw, na.rm = TRUE),
    cor_train_r2 = suppressWarnings(stats::cor(n_train_raw, r2, use = "complete.obs")),
    slope_r2_per_train_row = ifelse(
      sum(is.finite(n_train_raw) & is.finite(r2)) >= 3L,
      unname(stats::coef(stats::lm(r2 ~ n_train_raw))[["n_train_raw"]]),
      NA_real_
    ),
    .groups = "drop"
  )

ss_total_effect_base <- by_fold %>%
  dplyr::group_by(cv_type, method_name, n_folds, model) %>%
  dplyr::summarise(
    cor_ss_total_r2 = suppressWarnings(stats::cor(ss_total, r2, use = "complete.obs")),
    cor_y_sd_r2 = suppressWarnings(stats::cor(y_sd, r2, use = "complete.obs")),
    r2_sd = stats::sd(r2, na.rm = TRUE),
    ss_total_cv = stats::sd(ss_total, na.rm = TRUE) / mean(ss_total, na.rm = TRUE),
    .groups = "drop"
  )

ss_total_lm_effect <- by_fold %>%
  dplyr::group_by(cv_type, method_name, n_folds, model) %>%
  dplyr::group_modify(~ fit_lm_slope(.x, y_col = "r2", x_col = "ss_total")) %>%
  dplyr::ungroup() %>%
  dplyr::rename(
    lm_n = n_lm,
    lm_slope_r2_per_ss_total = slope,
    lm_slope_se = slope_se,
    lm_t = slope_t,
    lm_p_value = slope_p,
    lm_slope_ci_low = slope_ci_low,
    lm_slope_ci_high = slope_ci_high,
    lm_intercept = intercept
  )

ss_total_effect <- ss_total_effect_base %>%
  dplyr::left_join(
    ss_total_lm_effect,
    by = c("cv_type", "method_name", "n_folds", "model")
  )

fold_environment_effect <- by_fold %>%
  dplyr::group_by(cv_type, method_name, n_folds, model) %>%
  dplyr::summarise(
    cor_env_pairdist_r2 = suppressWarnings(stats::cor(env_mean_pairwise_dist, r2, use = "complete.obs")),
    cor_env_abs_z_r2 = suppressWarnings(stats::cor(env_mean_abs_z, r2, use = "complete.obs")),
    cor_env_predictor_sd_r2 = suppressWarnings(stats::cor(env_mean_predictor_sd, r2, use = "complete.obs")),
    cor_env_pairdist_ss_total = suppressWarnings(stats::cor(env_mean_pairwise_dist, ss_total, use = "complete.obs")),
    mean_env_pairdist = mean(env_mean_pairwise_dist, na.rm = TRUE),
    sd_env_pairdist = stats::sd(env_mean_pairwise_dist, na.rm = TRUE),
    mean_env_predictor_sd = mean(env_mean_predictor_sd, na.rm = TRUE),
    .groups = "drop"
  )

fold_environment_distribution <- by_fold %>%
  dplyr::group_by(cv_type, method_name, n_folds, model) %>%
  dplyr::summarise(
    median_env_pairdist = stats::median(env_mean_pairwise_dist, na.rm = TRUE),
    p01_env_pairdist = safe_quantile(env_mean_pairwise_dist, 0.01),
    p05_env_pairdist = safe_quantile(env_mean_pairwise_dist, 0.05),
    p10_env_pairdist = safe_quantile(env_mean_pairwise_dist, 0.10),
    p90_env_pairdist = safe_quantile(env_mean_pairwise_dist, 0.90),
    p95_env_pairdist = safe_quantile(env_mean_pairwise_dist, 0.95),
    p99_env_pairdist = safe_quantile(env_mean_pairwise_dist, 0.99),
    median_env_predictor_sd = stats::median(env_mean_predictor_sd, na.rm = TRUE),
    p01_env_predictor_sd = safe_quantile(env_mean_predictor_sd, 0.01),
    p05_env_predictor_sd = safe_quantile(env_mean_predictor_sd, 0.05),
    p10_env_predictor_sd = safe_quantile(env_mean_predictor_sd, 0.10),
    p90_env_predictor_sd = safe_quantile(env_mean_predictor_sd, 0.90),
    p95_env_predictor_sd = safe_quantile(env_mean_predictor_sd, 0.95),
    p99_env_predictor_sd = safe_quantile(env_mean_predictor_sd, 0.99),
    .groups = "drop"
  )

# Seed-count sensitivity: when do cumulative metrics stabilise?
seed_groups <- seed_level %>%
  dplyr::group_by(cv_type, method_name, n_folds, model) %>%
  dplyr::group_split()

seed_conv_rows <- list()
seed_plateau_rows <- list()
for (g in seed_groups) {
  keys <- g[1, c("cv_type", "method_name", "n_folds", "model"), drop = FALSE]

  pooled_conv <- compute_seed_convergence(
    vals = g$pooled_r2_reconstructed,
    seeds = g$fold_seed,
    metric_name = "pooled_r2",
    tol = seed_plateau_tolerance,
    min_n = seed_plateau_min_n,
    stable_steps = seed_plateau_stable_steps
  )
  meanfold_conv <- compute_seed_convergence(
    vals = g$mean_fold_r2,
    seeds = g$fold_seed,
    metric_name = "mean_fold_r2",
    tol = seed_plateau_tolerance,
    min_n = seed_plateau_min_n,
    stable_steps = seed_plateau_stable_steps
  )

  seed_conv_rows[[length(seed_conv_rows) + 1L]] <- dplyr::bind_rows(pooled_conv$series, meanfold_conv$series) %>%
    dplyr::mutate(
      cv_type = keys$cv_type[[1]],
      method_name = keys$method_name[[1]],
      n_folds = keys$n_folds[[1]],
      model = keys$model[[1]]
    )
  seed_plateau_rows[[length(seed_plateau_rows) + 1L]] <- dplyr::bind_rows(pooled_conv$plateau, meanfold_conv$plateau) %>%
    dplyr::mutate(
      cv_type = keys$cv_type[[1]],
      method_name = keys$method_name[[1]],
      n_folds = keys$n_folds[[1]],
      model = keys$model[[1]]
    )
}
seed_count_convergence <- dplyr::bind_rows(seed_conv_rows) %>%
  dplyr::select(cv_type, method_name, n_folds, model, metric, n_seeds, fold_seed_last,
                cumulative_mean, cumulative_sd, cumulative_se, delta_from_prev)
seed_count_plateau <- dplyr::bind_rows(seed_plateau_rows) %>%
  dplyr::select(cv_type, method_name, n_folds, model, metric, n_available_seeds, tolerance,
                plateau_n_seeds, plateau_fold_seed, plateau_value, plateau_reached)

write.csv(by_fold, file.path(out_dir, "sensitivity_by_fold.csv"), row.names = FALSE)
write.csv(by_pooled, file.path(out_dir, "sensitivity_pooled_by_seed.csv"), row.names = FALSE)
write.csv(summary_across_seeds, file.path(out_dir, "sensitivity_summary.csv"), row.names = FALSE)
write.csv(by_comp, file.path(out_dir, "sensitivity_fold_composition.csv"), row.names = FALSE)
write.csv(train_size_effect, file.path(out_dir, "sensitivity_train_size_effect.csv"), row.names = FALSE)
write.csv(ss_total_effect, file.path(out_dir, "sensitivity_ss_total_r2_effect.csv"), row.names = FALSE)
write.csv(by_env, file.path(out_dir, "sensitivity_fold_environment_by_fold.csv"), row.names = FALSE)
write.csv(fold_environment_effect, file.path(out_dir, "sensitivity_fold_environment_r2_effect.csv"), row.names = FALSE)
write.csv(fold_environment_distribution, file.path(out_dir, "sensitivity_fold_environment_distribution.csv"), row.names = FALSE)
write.csv(seed_count_convergence, file.path(out_dir, "sensitivity_seed_count_convergence.csv"), row.names = FALSE)
write.csv(seed_count_plateau, file.path(out_dir, "sensitivity_seed_count_plateau.csv"), row.names = FALSE)

cat("\nWrote sensitivity outputs to:\n", out_dir, "\n", sep = "")
cat("DONE\n")
