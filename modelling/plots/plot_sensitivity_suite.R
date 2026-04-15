## Standalone sensitivity plotting script.
##
## Sources at runtime:
## - modelling/R/init_repo.R (bootstrap; via seagrass_init_repo)
## - modelling/pipeline_config.R (run configuration, output directory resolution)
## - modelling/plots/plot_config.R (plot themes, model colours, linetypes)
##
## Reads CSVs from:
## - <run_output_dir>/sensitivity_suite/
## - latest available output/tuning_seed_sweep_runs/sweep_*/sensitivity_tuning_seed_sweep_summary.csv

if (!exists("seagrass_init_repo", mode = "function", inherits = TRUE)) source("modelling/R/init_repo.R")
project_root <- seagrass_init_repo(
  packages = c("dplyr", "ggplot2", "readr", "tidyr"),
  source_files = c("modelling/pipeline_config.R", "modelling/plots/plot_config.R"),
  include_helpers = TRUE,
  require_core_inputs = FALSE,
  check_renv = FALSE
)
library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)

plot_log_sp <- function(message_text) {
  message(sprintf("[plot:plot_sensitivity_suite] %s", message_text))
}

rel_path_sp <- function(path) {
  p <- normalizePath(path, winslash = "/", mustWork = FALSE)
  root <- normalizePath(project_root, winslash = "/", mustWork = FALSE)
  sub(paste0("^", root, "/?"), "", p)
}

filter_models_required <- function(df, model_list, df_label) {
  if (!"model" %in% names(df)) return(df)
  requested <- unique(as.character(model_list))
  requested <- requested[nzchar(requested)]
  if (length(requested) == 0L) {
    stop("model_list from current config is empty; nothing to plot.", call. = FALSE)
  }
  available <- unique(as.character(df$model))
  available <- available[nzchar(available)]
  missing_requested <- setdiff(requested, available)
  if (length(missing_requested) > 0L) {
    stop(
      "Requested model_list contains models not present in ", df_label, ".\n",
      "Requested: ", paste(requested, collapse = ", "), "\n",
      "Available: ", paste(available, collapse = ", "), "\n",
      "Missing: ", paste(missing_requested, collapse = ", "),
      call. = FALSE
    )
  }
  extra_available <- setdiff(available, requested)
  if (length(extra_available) > 0L) {
    plot_log_sp(paste0(df_label, " has extra models that will be excluded: ", paste(extra_available, collapse = ", ")))
  }
  df[df$model %in% requested, , drop = FALSE]
}

cfg <- get_pipeline_config()
apply_pipeline_defaults(
  cfg,
  c("cv_regime_name", "cv_type_label", "robust_fold_seed_list", "eval_fold_seed_list", "dpi", "show_titles", "model_list", "run_output_dir"),
  envir = .GlobalEnv
)

combine_env_r2_plots <- function(pairdist_plot, predictor_sd_plot, out_file,
                                 width = 16, height = 8, dpi = 300) {
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    warning("patchwork not installed; skipping combined environmental panel.")
    return(invisible(NULL))
  }
  combined <- (pairdist_plot + predictor_sd_plot) +
    patchwork::plot_layout(ncol = 1, guides = "collect") &
    theme(legend.position = "bottom")
  ggsave(filename = out_file, plot = combined, width = width, height = height, dpi = dpi)
  invisible(combined)
}

combine_fold_comp_plots <- function(ss_total_plot, y_sd_plot, out_file,
                                    width = 12, height = 11, dpi = 300) {
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    warning("patchwork not installed; skipping combined fold-composition panel.")
    return(invisible(NULL))
  }
  combined <- (ss_total_plot + y_sd_plot) +
    patchwork::plot_layout(ncol = 1, guides = "collect") &
    theme(legend.position = "bottom")
  ggsave(filename = out_file, plot = combined, width = width, height = height, dpi = dpi)
  invisible(combined)
}

cv_regime_name <- get("cv_regime_name", envir = .GlobalEnv)
cv_pipeline_dir <- file.path(project_root, "output", cv_regime_name, "cv_pipeline")

resolve_robust_eval_dir <- function() {
  robust_seed_str <- paste(as.integer(get("robust_fold_seed_list", envir = .GlobalEnv)), collapse = "-")
  preferred_seed_dir <- file.path(project_root, "output", paste0(get("cv_regime_name", envir = .GlobalEnv), "_", robust_seed_str))
  if (dir.exists(preferred_seed_dir)) {
    return(preferred_seed_dir)
  }

  run_output_dir <- get0("run_output_dir", envir = .GlobalEnv, ifnotfound = NA_character_)
  if (is.na(run_output_dir) || !nzchar(as.character(run_output_dir))) {
    candidates <- Sys.glob(file.path(project_root, "output", paste0(get("cv_regime_name", envir = .GlobalEnv), "_*")))
    candidates <- candidates[dir.exists(candidates)]
    if (length(candidates) == 0L) {
      stop("Could not resolve run_output_dir from current config.", call. = FALSE)
    }
    run_output_dir <- candidates[which.max(file.info(candidates)$mtime)]
    plot_log_sp(paste0("run_output_dir not set; falling back to most recent candidate: ", run_output_dir))
  }
  run_output_dir <- as.character(run_output_dir)
  if (!dir.exists(run_output_dir)) {
    stop("run_output_dir does not exist: ", run_output_dir)
  }
  run_output_dir
}

robust_eval_dir <- resolve_robust_eval_dir()
sens_dir <- file.path(robust_eval_dir, "sensitivity_suite")
dir.create(sens_dir, recursive = TRUE, showWarnings = FALSE)
plots_dir <- file.path(sens_dir, "plots")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

resolve_input_file <- function(filename) {
  p <- file.path(sens_dir, filename)
  if (file.exists(p)) return(p)
  NA_character_
}

required_inputs <- c(
  "sensitivity_summary.csv",
  "sensitivity_by_fold.csv",
  "sensitivity_train_size_effect.csv",
  "sensitivity_ss_total_r2_effect.csv",
  "sensitivity_seed_count_convergence.csv",
  "sensitivity_seed_count_plateau.csv",
  "sensitivity_fold_environment_by_fold.csv",
  "sensitivity_fold_environment_r2_effect.csv",
  "sensitivity_fold_environment_distribution.csv"
)
plot_log_sp("Sources: modelling/R/init_repo.R, modelling/pipeline_config.R, modelling/plots/plot_config.R")
plot_log_sp(paste0("Input sensitivity directory: ", rel_path_sp(sens_dir)))
plot_log_sp(paste0("Output plot directory: ", rel_path_sp(plots_dir)))
plot_log_sp(paste0("Requested model_list: ", paste(as.character(model_list), collapse = ", ")))
plot_log_sp(
  paste0(
    "Required CSVs: ",
    paste(file.path(rel_path_sp(sens_dir), required_inputs), collapse = ", "),
    " ; Optional sweep CSV: sensitivity_tuning_seed_sweep_summary.csv"
  )
)
input_paths <- setNames(vapply(required_inputs, resolve_input_file, character(1)), required_inputs)
path_present <- !is.na(input_paths) & nzchar(input_paths) & file.exists(input_paths)
missing_inputs <- names(input_paths)[!path_present]
has_full_sensitivity_inputs <- all(path_present)

tuning_seed_sweep_path <- {
  sweep_run_csvs <- Sys.glob(
    file.path(project_root, "output", "tuning_seed_sweep_runs", "sweep_*", "sensitivity_tuning_seed_sweep_summary.csv")
  )
  sweep_run_csvs <- sweep_run_csvs[file.exists(sweep_run_csvs)]
  cands <- c(
    file.path(sens_dir, "sensitivity_tuning_seed_sweep_summary.csv"),
    sweep_run_csvs
  )
  hit <- cands[file.exists(cands)]
  if (length(hit) == 0L) {
    NA_character_
  } else if (length(hit) == 1L) {
    hit[[1L]]
  } else {
    hit[which.max(file.info(hit)$mtime)]
  }
}

if (!has_full_sensitivity_inputs && (is.na(tuning_seed_sweep_path) || !file.exists(tuning_seed_sweep_path))) {
  stop(
    "Missing full sensitivity inputs and no sweep summary found.\n",
    "Missing:\n  ", paste(missing_inputs, collapse = "\n  "),
    "\nSearched in:\n  ", sens_dir
  )
}
if (!has_full_sensitivity_inputs) {
  plot_log_sp(paste0(
    "Missing full sensitivity-suite CSVs; generating sweep plots only.\nMissing:\n  ",
    paste(missing_inputs, collapse = "\n  ")
  ))
}

if (has_full_sensitivity_inputs) {
  summary_df <- readr::read_csv(input_paths[["sensitivity_summary.csv"]], show_col_types = FALSE)
  by_fold_df <- readr::read_csv(input_paths[["sensitivity_by_fold.csv"]], show_col_types = FALSE)
  train_effect_df <- readr::read_csv(input_paths[["sensitivity_train_size_effect.csv"]], show_col_types = FALSE)
  ss_effect_df <- readr::read_csv(input_paths[["sensitivity_ss_total_r2_effect.csv"]], show_col_types = FALSE)
  seed_conv_df <- readr::read_csv(input_paths[["sensitivity_seed_count_convergence.csv"]], show_col_types = FALSE)
  seed_plateau_df <- readr::read_csv(input_paths[["sensitivity_seed_count_plateau.csv"]], show_col_types = FALSE)
  env_by_fold_df <- readr::read_csv(input_paths[["sensitivity_fold_environment_by_fold.csv"]], show_col_types = FALSE)
  env_effect_df <- readr::read_csv(input_paths[["sensitivity_fold_environment_r2_effect.csv"]], show_col_types = FALSE)
  env_dist_df <- readr::read_csv(input_paths[["sensitivity_fold_environment_distribution.csv"]], show_col_types = FALSE)

  summary_df <- filter_models_required(summary_df, model_list, "sensitivity_summary.csv")
  by_fold_df <- filter_models_required(by_fold_df, model_list, "sensitivity_by_fold.csv")
  train_effect_df <- filter_models_required(train_effect_df, model_list, "sensitivity_train_size_effect.csv")
  ss_effect_df <- filter_models_required(ss_effect_df, model_list, "sensitivity_ss_total_r2_effect.csv")
  seed_conv_df <- filter_models_required(seed_conv_df, model_list, "sensitivity_seed_count_convergence.csv")
  seed_plateau_df <- filter_models_required(seed_plateau_df, model_list, "sensitivity_seed_count_plateau.csv")
  env_by_fold_df <- filter_models_required(env_by_fold_df, model_list, "sensitivity_fold_environment_by_fold.csv")
  env_effect_df <- filter_models_required(env_effect_df, model_list, "sensitivity_fold_environment_r2_effect.csv")
  env_dist_df <- filter_models_required(env_dist_df, model_list, "sensitivity_fold_environment_distribution.csv")
}
tuning_seed_sweep_df <- if (!is.na(tuning_seed_sweep_path) && file.exists(tuning_seed_sweep_path)) {
  z <- readr::read_csv(tuning_seed_sweep_path, show_col_types = FALSE)
  filter_models_required(z, model_list, basename(tuning_seed_sweep_path))
} else {
  NULL
}

if (has_full_sensitivity_inputs) {
# 1) Fold-count sensitivity: pooled vs mean-fold R2 by k
p_fold_count <- ggplot(
  summary_df,
  aes(x = n_folds, group = model, color = model)
) +
  geom_line(aes(y = mean_pooled_r2, linetype = "Pooled R2"), linewidth = 0.9) +
  geom_point(aes(y = mean_pooled_r2, shape = "Pooled R2"), size = 2.1) +
  geom_line(aes(y = mean_mean_fold_r2, linetype = "Mean-fold R2"), linewidth = 0.8, alpha = 0.8) +
  geom_point(aes(y = mean_mean_fold_r2, shape = "Mean-fold R2"), size = 1.8, alpha = 0.8) +
  facet_wrap(~ cv_type, scales = "free_y", labeller = labeller(cv_type = function(x) ifelse(x == "pixel_grouped", "", x))) +
  scale_color_manual(values = MODEL_COLOURS, name = "Model", drop = FALSE) +
  scale_linetype_manual(name = expression("R"^2*" type"), values = METRIC_LINESTYLES[c("Pooled R2", "Mean-fold R2")]) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey40") +
  scale_shape_manual(
    name = expression("R"^2*" type"),
    values = c("Pooled R2" = 16, "Mean-fold R2" = 17)
  ) +
  labs(
    title = if (show_titles) expression("Sensitivity of " * R^2 * " to fold count") else NULL,
    x = "Number of folds",
    y = expression("R"^2)
  )

ggsave(
  filename = file.path(plots_dir, "sensitivity_r2_fold_count.png"),
  plot = p_fold_count, width = 11, height = 6, dpi = dpi
)

# 2) Gap between pooled and mean-fold R2 vs k
p_gap <- ggplot(
  summary_df,
  aes(x = n_folds, y = mean_delta_pooled_minus_meanfold, color = model)
) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey40") +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.1) +
  facet_wrap(~ cv_type, scales = "free_y", labeller = labeller(cv_type = function(x) ifelse(x == "pixel_grouped", "", x))) +
  scale_color_manual(values = MODEL_COLOURS, name = "Model", drop = FALSE) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey40") +
  labs(
    title = if (show_titles) expression("Gap: pooled " * R^2 * " - mean-fold " * R^2) else NULL,
    x = "Number of folds",
    y = expression("Change in " * R^2 * " from pooled to mean-fold " * R^2)
  )

ggsave(
  filename = file.path(plots_dir, "sensitivity_r2_gap_by_fold_count.png"),
  plot = p_gap, width = 11, height = 6, dpi = dpi
)

# 3) Training size vs fold R2 (raw fold points + linear trend)
p_train_size <- ggplot(
  by_fold_df,
  aes(x = n_train_raw, y = r2, color = model)
) +
  geom_point(alpha = 0.55, size = 1.6) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.9) +
  facet_grid(cv_type ~ n_folds, scales = "free_x", labeller = labeller(cv_type = function(x) ifelse(x == "pixel_grouped", "", x))) +
  scale_color_manual(values = MODEL_COLOURS, name = "Model", drop = FALSE) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey40") +
  labs(
    title = if (show_titles) expression("Training-set size vs fold " * R^2) else NULL,
    # Each point is one model-fold-seed result
    x = "Training rows in fold",
    y = expression("Fold " * R^2)
  )

ggsave(
  filename = file.path(plots_dir, "sensitivity_train_size_vs_r2.png"),
  plot = p_train_size, width = 13, height = 8, dpi = dpi
)

# 4) Fold composition sensitivity: SS_total and y spread against fold R2
p_ss_fold <- ggplot(
  by_fold_df,
  aes(x = ss_total, y = r2, color = model)
) +
  geom_point(alpha = 0.55, size = 1.6) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.9) +
  facet_grid(cv_type ~ n_folds, scales = "free",  labeller = labeller(cv_type = function(x) ifelse(x == "pixel_grouped", "", x))) +
  scale_color_manual(values = MODEL_COLOURS, name = "Model", drop = FALSE) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey40") +
  labs(
    title = expression("A. Fold composition sensitivity: " * "SS"[total] * " vs fold " * R^2),
    x = expression("Fold " * "SS"[total]),
    y = expression("Fold " * R^2)
  )

ggsave(
  filename = file.path(plots_dir, "sensitivity_ss_total_vs_r2.png"),
  plot = p_ss_fold, width = 13, height = 8, dpi = dpi
)

p_yspread_fold <- ggplot(
  by_fold_df,
  aes(x = y_sd, y = r2, color = model)
) +
  geom_point(alpha = 0.55, size = 1.6) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.9) +
  facet_grid(cv_type ~ n_folds, scales = "free", labeller = labeller(cv_type = function(x) ifelse(x == "pixel_grouped", "", x))) +
  scale_color_manual(values = MODEL_COLOURS, name = "Model", drop = FALSE) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey40") +
  labs(
    title = expression("B. Fold composition sensitivity: SD of observed values vs fold " * R^2),
    x = "Fold SD of observed values",
    y = expression("Fold " * R^2)
  )

ggsave(
  filename = file.path(plots_dir, "sensitivity_y_sd_vs_r2.png"),
  plot = p_yspread_fold, width = 13, height = 8, dpi = dpi
)

combine_fold_comp_plots(
  ss_total_plot = p_ss_fold,
  y_sd_plot = p_yspread_fold,
  out_file = file.path(plots_dir, "sensitivity_fold_comp_r2_combined.png"),
  width = 12, height = 11, dpi = 300
)

# 5) Aggregated effect summaries (from pre-computed sensitivity tables)
p_train_effect <- ggplot(
  train_effect_df,
  aes(x = n_folds, y = cor_train_r2, color = model)
) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey40") +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.1) +
  facet_wrap(~ cv_type, scales = "free_y", labeller = labeller(cv_type = function(x) ifelse(x == "pixel_grouped", "", x))) +
  scale_color_manual(values = MODEL_COLOURS, name = "Model", drop = FALSE) +
  labs(
    title = if (show_titles) "Aggregate training-size effect" else NULL,
    x = "Number of folds",
    y = expression("Correlation between training size and fold " * R^2)
  )

ggsave(
  filename = file.path(plots_dir, "sensitivity_train_size_effect_summary.png"),
  plot = p_train_effect, width = 11, height = 6, dpi = dpi
)

p_ss_effect <- ggplot(
  ss_effect_df,
  aes(x = n_folds, y = cor_ss_total_r2, color = model)
) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey40") +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.1) +
  facet_wrap(~ cv_type, scales = "free_y", labeller = labeller(cv_type = function(x) ifelse(x == "pixel_grouped", "", x))) +
  scale_color_manual(values = MODEL_COLOURS, name = "Model", drop = FALSE) +
  labs(
    title = if (show_titles) "Aggregate fold-composition effect" else NULL,
    x = "Number of folds",
    y = expression("Correlation between fold " * "SS"[total] * " and fold " * R^2)
  )

ggsave(
  filename = file.path(plots_dir, "sensitivity_ss_total_effect_summary.png"),
  plot = p_ss_effect, width = 11, height = 6, dpi = dpi
)

# 6) Seed-count sensitivity: convergence and plateau
# Strip labels: metric is pooled_r2 vs mean_fold_r2 (see sensitivity_suite.R).
metric_strip_labels <- c(
  mean_fold_r2 = "Mean fold R\u00b2",
  pooled_r2 = "Pooled R\u00b2"
)
p_seed_conv <- ggplot(
  seed_conv_df,
  aes(x = n_seeds, y = cumulative_mean, color = model)
) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.9) +
  geom_errorbar(
    aes(ymin = cumulative_mean - cumulative_se, ymax = cumulative_mean + cumulative_se),
    width = 0.12, alpha = 0.5
  ) +
  facet_grid(
    metric ~ n_folds,
    scales = "free_y",
    labeller = labeller(
      metric = as_labeller(metric_strip_labels, default = label_value),
      n_folds = function(x) paste0("k = ", x)
    )
  ) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey40") +
  scale_color_manual(values = MODEL_COLOURS, name = "Model", drop = FALSE) +
  labs(
    title = if (show_titles) "Seed-count sensitivity: cumulative performance convergence" else NULL,
    # subtitle = "Error bars are cumulative standard errors across seeds",
    x = "Number of seeds included",
    y = "Cumulative mean (R\u00b2)"
  )

ggsave(
  filename = file.path(plots_dir, "sensitivity_seed_count_convergence.png"),
  plot = p_seed_conv, width = 13, height = 8, dpi = dpi
)

p_seed_plateau <- ggplot(
  seed_plateau_df,
  aes(x = n_folds, y = plateau_n_seeds, color = model, group = model)
) +
  geom_point(size = 2.2, alpha = 0.9) +
  geom_line(linewidth = 0.9, alpha = 0.9) +
  facet_wrap(
    ~ metric,
    scales = "free_y",
    labeller = labeller(metric = as_labeller(metric_strip_labels, default = label_value))
  ) +
  scale_color_manual(values = MODEL_COLOURS, name = "Model", drop = FALSE) +
  labs(
    title = if (show_titles) "Estimated seeds needed to reach performance plateau" else NULL,
    # subtitle = "Plateau defined by cumulative-mean changes below tolerance for consecutive steps",
    x = "Number of folds",
    y = "Estimated seeds-to-plateau"
  )

ggsave(
  filename = file.path(plots_dir, "sensitivity_seed_count_plateau.png"),
  plot = p_seed_plateau, width = 11, height = 6, dpi = dpi
)

# 7) Environmental summaries
p_env_pairdist_vs_r2 <- ggplot(
  by_fold_df,
  aes(x = env_mean_pairwise_dist, y = r2, color = model)
) +
  geom_point(alpha = 0.3, size = 1.6) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.9) +
  facet_grid(cv_type ~ n_folds, scales = "free", labeller = labeller(cv_type = function(x) ifelse(x == "pixel_grouped", "", x))) +
  scale_color_manual(values = MODEL_COLOURS, name = "Model", drop = FALSE) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey40") +
  labs(
    title = expression("A. Environmental spread vs fold " * R^2),
    # subtitle = "Environmental spread measured by mean pairwise distance in scaled predictor space",
    x = "Fold mean pairwise environmental distance",
    y = expression("Fold " * R^2)
  )

ggsave(
  filename = file.path(plots_dir, "sensitivity_env_pairdist_vs_r2.png"),
  plot = p_env_pairdist_vs_r2, width = 13, height = 8, dpi = dpi
)

p_env_predictor_sd_vs_r2 <- ggplot(
  by_fold_df,
  aes(x = env_mean_predictor_sd, y = r2, color = model)
) +
  geom_point(alpha = 0.3, size = 1.6) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.9) +
  facet_grid(cv_type ~ n_folds, scales = "free", labeller = labeller(cv_type = function(x) ifelse(x == "pixel_grouped", "", x))) +
  scale_color_manual(values = MODEL_COLOURS, name = "Model", drop = FALSE) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey40") +
  labs(
    title = expression("B. Within-fold predictor spread vs fold " * R^2),
    x = "Mean predictor SD within fold (scaled predictors)",
    y = expression("Fold " * R^2)
  )

ggsave(
  filename = file.path(plots_dir, "sensitivity_env_predictor_sd_vs_r2.png"),
  plot = p_env_predictor_sd_vs_r2, width = 13, height = 8, dpi = dpi
)

combine_env_r2_plots(
  pairdist_plot = p_env_pairdist_vs_r2,
  predictor_sd_plot = p_env_predictor_sd_vs_r2,
  out_file = file.path(plots_dir, "sensitivity_env_r2_combined.png"),
  width = 12, height = 11, dpi = 300
)

p_env_effect <- ggplot(
  env_effect_df,
  aes(x = n_folds, y = cor_env_pairdist_r2, color = model)
) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey40") +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.1) +
  facet_wrap(~ cv_type, scales = "free_y", labeller = labeller(cv_type = function(x) ifelse(x == "pixel_grouped", "", x))) +
  scale_color_manual(values = MODEL_COLOURS, name = "Model", drop = FALSE) +
  labs(
    title = if (show_titles) expression("Aggregate environmental-spread effect on " * R^2) else NULL,
    x = "Number of folds",
    y = expression("Correlation between fold environmental pairwise distance and fold " * R^2)
  )

ggsave(
  filename = file.path(plots_dir, "sensitivity_env_effect_summary.png"),
  plot = p_env_effect, width = 11, height = 6, dpi = dpi
)

env_dist_long <- env_dist_df %>%
  dplyr::select(
    cv_type, n_folds, model,
    p01_env_pairdist, p05_env_pairdist, p10_env_pairdist,
    median_env_pairdist, p90_env_pairdist, p95_env_pairdist, p99_env_pairdist
  ) %>%
  tidyr::pivot_longer(
    cols = c(
      p01_env_pairdist, p05_env_pairdist, p10_env_pairdist,
      median_env_pairdist, p90_env_pairdist, p95_env_pairdist, p99_env_pairdist
    ),
    names_to = "quantile",
    values_to = "env_pairdist"
  ) %>%
  dplyr::mutate(
    quantile = dplyr::recode(
      quantile,
      p01_env_pairdist = "1st",
      p05_env_pairdist = "5th",
      p10_env_pairdist = "10th",
      median_env_pairdist = "median",
      p90_env_pairdist = "90th",
      p95_env_pairdist = "95th",
      p99_env_pairdist = "99th"
    ),
    quantile = factor(quantile, levels = c("1st", "5th", "10th", "median", "90th", "95th", "99th"))
  )

p_env_dist <- ggplot(
  env_dist_long,
  aes(x = n_folds, y = env_pairdist,
      linetype = quantile,
      group = interaction(model, quantile))
) +
  geom_line(linewidth = 0.9, alpha = 0.9) +
  facet_wrap(~ cv_type, scales = "free_y", labeller = labeller(cv_type = function(x) ifelse(x == "pixel_grouped", "", x))) +
  scale_linetype_manual(
    values = c(
      "1st" = "dotted",
      "5th" = "dotdash",
      "10th" = "dashed",
      "median" = "solid",
      "90th" = "dashed",
      "95th" = "dotdash",
      "99th" = "dotted"
    ),
    name = "Quantile"
  ) +
  labs(
    title = if (show_titles) "Distribution of fold environmental spread across k" else NULL,
    # subtitle = "Percentile curves shown from the 1st to 99th percentile",
    x = "Number of folds",
    y = "Fold mean pairwise environmental distance",
    # color = "Model"
  )

ggsave(
  filename = file.path(plots_dir, "sensitivity_env_distribution_percentiles.png"),
  plot = p_env_dist, width = 11, height = 6, dpi = 300
)
}

# 8) Optional: tuning seed-count sweep performance plots
if (!is.null(tuning_seed_sweep_df) && nrow(tuning_seed_sweep_df) > 0L &&
    !is.na(tuning_seed_sweep_path) && nzchar(tuning_seed_sweep_path)) {
  source(file.path(project_root, "modelling/plots/plot_tuning_seed_sweep.R"))
  sweep_plot_dir <- dirname(tuning_seed_sweep_path)
  dir.create(sweep_plot_dir, recursive = TRUE, showWarnings = FALSE)
  plot_tuning_seed_sweep_summary(
    sweep_df = tuning_seed_sweep_df,
    out_dir = sweep_plot_dir,
    project_root = project_root,
    dpi = dpi,
    model_list = model_list
  )
} else {
  plot_log_sp("No tuning seed sweep summary found; skipping sweep R2/RMSE plots.")
}

plot_log_sp(paste0("Wrote sensitivity plots to: ", rel_path_sp(plots_dir)))