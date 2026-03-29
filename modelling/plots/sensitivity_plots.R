project_root <- here::here()
setwd(project_root)

library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
source(file.path(project_root, "modelling/config/pipeline_config.R"))
source(file.path(project_root, "modelling/R/plot_config.R"))

cfg <- get_pipeline_config()
apply_pipeline_defaults(
  cfg,
  c("cv_regime_name", "cv_type_label", "robust_fold_seed_list", "eval_fold_seed_list"),
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
  run_output_dir <- get0("run_output_dir", envir = .GlobalEnv, ifnotfound = NULL)
  if (!is.null(run_output_dir) && dir.exists(run_output_dir)) {
    return(run_output_dir)
  }

  env_seeds <- get0("robust_fold_seed_list", envir = .GlobalEnv, ifnotfound = NULL)
  eval_seeds <- get0("eval_fold_seed_list", envir = .GlobalEnv, ifnotfound = NULL)
  cv_type_label <- get0("cv_type_label", envir = .GlobalEnv, ifnotfound = "pixel_grouped")
  if (!is.null(env_seeds) && length(env_seeds) > 0L && !is.null(eval_seeds) && length(eval_seeds) > 0L) {
    target <- file.path(
      cv_pipeline_dir,
      build_seeded_run_folder_name(
        cv_type_label = cv_type_label,
        folder_type = "evaluation",
        repeat_seed_list = eval_seeds,
        robust_seed_list = env_seeds,
        include_seed_values = TRUE
      )
    )
    if (dir.exists(target)) return(target)
  }

  if (!is.null(env_seeds) && length(env_seeds) > 0L) {
    seeds_str <- paste(as.integer(env_seeds), collapse = "-")
    target <- file.path(
      cv_pipeline_dir,
      paste0("robust_pixel_grouped_evaluation_robustSeeds_", seeds_str)
    )
    if (dir.exists(target)) return(target)
  }
  candidates <- c(
    Sys.glob(file.path(cv_pipeline_dir, "pixel_grouped_evaluation_*x*_seeds*")),
    Sys.glob(file.path(cv_pipeline_dir, "robust_pixel_grouped_evaluation_robustSeeds_*"))
  )
  candidates <- unique(candidates)
  candidates <- candidates[dir.exists(candidates)]
  if (length(candidates) == 0L) {
    stop("No robust evaluation directory found in ", cv_pipeline_dir)
  }
  candidates[which.max(file.info(candidates)$mtime)]
}

robust_eval_dir <- resolve_robust_eval_dir()
sens_dir <- file.path(robust_eval_dir, "sensitivity_suite")
legacy_sens_dir <- file.path(cv_pipeline_dir, "sensitivity_suite")
dir.create(sens_dir, recursive = TRUE, showWarnings = FALSE)

resolve_input_file <- function(filename) {
  candidates <- c(file.path(sens_dir, filename), file.path(legacy_sens_dir, filename))
  for (p in candidates) {
    if (file.exists(p)) return(p)
  }
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
input_paths <- setNames(vapply(required_inputs, resolve_input_file, character(1)), required_inputs)
path_present <- !is.na(input_paths) & nzchar(input_paths) & file.exists(input_paths)
missing_inputs <- names(input_paths)[!path_present]
has_full_sensitivity_inputs <- all(path_present)

tuning_seed_sweep_path <- {
  cands <- c(file.path(sens_dir, "sensitivity_tuning_seed_sweep_summary.csv"),
             file.path(legacy_sens_dir, "sensitivity_tuning_seed_sweep_summary.csv"))
  hit <- cands[file.exists(cands)]
  if (length(hit) > 0L) hit[[1]] else NA_character_
}

if (!has_full_sensitivity_inputs && (is.na(tuning_seed_sweep_path) || !file.exists(tuning_seed_sweep_path))) {
  stop(
    "Missing full sensitivity inputs and no sweep summary found.\n",
    "Missing:\n  ", paste(missing_inputs, collapse = "\n  "),
    "\nSearched in:\n  ", sens_dir, "\n  ", legacy_sens_dir
  )
}
if (!has_full_sensitivity_inputs) {
  message(
    "Missing full sensitivity-suite CSVs; generating sweep plots only.\n",
    "Missing:\n  ", paste(missing_inputs, collapse = "\n  ")
  )
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
}
tuning_seed_sweep_df <- if (!is.na(tuning_seed_sweep_path) && file.exists(tuning_seed_sweep_path)) {
  readr::read_csv(tuning_seed_sweep_path, show_col_types = FALSE)
} else {
  NULL
}

theme_set(theme_minimal(base_size = 12))

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
  scale_color_discrete(name = "Model") +
  scale_linetype_manual(name = "R2 type", values = METRIC_LINESTYLES[c("Pooled R2", "Mean-fold R2")]) +
  scale_shape_manual(
    name = "R2 type",
    values = c("Pooled R2" = 16, "Mean-fold R2" = 17)
  ) +
  labs(
    title = "Sensitivity of R2 to fold count",
    x = "Number of folds",
    y = "R2"
  )

ggsave(
  filename = file.path(sens_dir, "sensitivity_r2_fold_count.png"),
  plot = p_fold_count, width = 11, height = 6, dpi = 220
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
  scale_color_discrete(name = "Model") +

  labs(
    title = ("Gap: pooled R2 - mean-fold R2"),
    x = "Number of folds",
    y = ("Change in R2 from pooled to mean-fold R2")
  )

ggsave(
  filename = file.path(sens_dir, "sensitivity_r2_gap_by_fold_count.png"),
  plot = p_gap, width = 11, height = 6, dpi = 220
)

# 3) Training size vs fold R2 (raw fold points + linear trend)
p_train_size <- ggplot(
  by_fold_df,
  aes(x = n_train_raw, y = r2, color = model)
) +
  geom_point(alpha = 0.55, size = 1.6) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.9) +
  facet_grid(cv_type ~ n_folds, scales = "free_x", labeller = labeller(cv_type = function(x) ifelse(x == "pixel_grouped", "", x))) +
  scale_color_discrete(name = "Model") +
  labs(
    title = "Training-set size vs fold R2",
    subtitle = "Each point is one model-fold-seed result",
    x = "Training rows in fold",
    y = "Fold R2"
  )

ggsave(
  filename = file.path(sens_dir, "sensitivity_train_size_vs_r2.png"),
  plot = p_train_size, width = 13, height = 8, dpi = 220
)

# 4) Fold composition sensitivity: SS_total and y spread against fold R2
p_ss_fold <- ggplot(
  by_fold_df,
  aes(x = ss_total, y = r2, color = model)
) +
  geom_point(alpha = 0.55, size = 1.6) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.9) +
  facet_grid(cv_type ~ n_folds, scales = "free",  labeller = labeller(cv_type = function(x) ifelse(x == "pixel_grouped", "", x))) +
  scale_color_discrete(name = "Model") +
  labs(
    title = "A. Fold composition sensitivity: SS_total vs fold R2",
    x = "Fold SS_total",
    y = ("Fold R2")
  )

ggsave(
  filename = file.path(sens_dir, "sensitivity_ss_total_vs_r2.png"),
  plot = p_ss_fold, width = 13, height = 8, dpi = 220
)

p_yspread_fold <- ggplot(
  by_fold_df,
  aes(x = y_sd, y = r2, color = model)
) +
  geom_point(alpha = 0.55, size = 1.6) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.9) +
  facet_grid(cv_type ~ n_folds, scales = "free", labeller = labeller(cv_type = function(x) ifelse(x == "pixel_grouped", "", x))) +
  scale_color_discrete(name = "Model") +
  labs(
    title = "B. Fold composition sensitivity: sd(y) vs fold R2",
    x = "Fold sd(observed)",
    y = ("Fold R2")
  )

ggsave(
  filename = file.path(sens_dir, "sensitivity_y_sd_vs_r2.png"),
  plot = p_yspread_fold, width = 13, height = 8, dpi = 220
)

combine_fold_comp_plots(
  ss_total_plot = p_ss_fold,
  y_sd_plot = p_yspread_fold,
  out_file = file.path(sens_dir, "sensitivity_fold_comp_r2_combined.png"),
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
  scale_color_discrete(name = "Model") +
  labs(
    title = "Aggregate training-size effect",
    subtitle = "Correlation between training size and fold R2",
    x = "Number of folds",
    y = "Correlation between training size and fold R2"
  )

ggsave(
  filename = file.path(sens_dir, "sensitivity_train_size_effect_summary.png"),
  plot = p_train_effect, width = 11, height = 6, dpi = 220
)

p_ss_effect <- ggplot(
  ss_effect_df,
  aes(x = n_folds, y = cor_ss_total_r2, color = model)
) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey40") +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.1) +
  facet_wrap(~ cv_type, scales = "free_y", labeller = labeller(cv_type = function(x) ifelse(x == "pixel_grouped", "", x))) +
  scale_color_discrete(name = "Model") +
  labs(
    title = "Aggregate fold-composition effect",
    subtitle = "Correlation between fold SS_total and fold R2 (averaged across seeds)",
    x = "Number of folds",
    y = "Correlation between fold SS_total and fold R2"
  )

ggsave(
  filename = file.path(sens_dir, "sensitivity_ss_total_effect_summary.png"),
  plot = p_ss_effect, width = 11, height = 6, dpi = 220
)

# 6) Seed-count sensitivity: convergence and plateau
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
  facet_grid(metric ~ n_folds, scales = "free_y", labeller = labeller(metric = function(x) ifelse(x == "cumulative_mean", "Cumulative mean", x))) +
  scale_color_discrete(name = "Model") +
  labs(
    title = "Seed-count sensitivity: cumulative performance convergence",
    subtitle = "Error bars are cumulative standard errors across seeds",
    x = "Number of seeds included",
    y = "Cumulative mean metric"
  )

ggsave(
  filename = file.path(sens_dir, "sensitivity_seed_count_convergence.png"),
  plot = p_seed_conv, width = 13, height = 8, dpi = 220
)

p_seed_plateau <- ggplot(
  seed_plateau_df,
  aes(x = n_folds, y = plateau_n_seeds, color = model, group = model)
) +
  geom_point(size = 2.2, alpha = 0.9) +
  geom_line(linewidth = 0.9, alpha = 0.9) +
  facet_wrap(~ metric, scales = "free_y", labeller = labeller(metric = function(x) ifelse(x == "cumulative_mean", "Cumulative mean", x))) +
  scale_color_discrete(name = "Model") +
  labs(
    title = "Estimated seeds needed to reach plateau",
    subtitle = "Plateau defined by cumulative-mean changes below tolerance for consecutive steps",
    x = "Number of folds",
    y = "Estimated seeds-to-plateau"
  )

ggsave(
  filename = file.path(sens_dir, "sensitivity_seed_count_plateau.png"),
  plot = p_seed_plateau, width = 11, height = 6, dpi = 220
)

# 7) Environmental summaries
p_env_pairdist_vs_r2 <- ggplot(
  by_fold_df,
  aes(x = env_mean_pairwise_dist, y = r2, color = model)
) +
  geom_point(alpha = 0.3, size = 1.6) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.9) +
  facet_grid(cv_type ~ n_folds, scales = "free", labeller = labeller(cv_type = function(x) ifelse(x == "pixel_grouped", "", x))) +
  scale_color_discrete(name = "Model") +
  labs(
    title = "A. Environmental spread vs fold R2",
    subtitle = "Environmental spread measured by mean pairwise distance in scaled predictor space",
    x = "Fold mean pairwise environmental distance",
    y = "Fold R2"
  )
p_env_pairdist_vs_r2
ggsave(
  filename = file.path(sens_dir, "sensitivity_env_pairdist_vs_r2.png"),
  plot = p_env_pairdist_vs_r2, width = 13, height = 8, dpi = 220
)

p_env_predictor_sd_vs_r2 <- ggplot(
  by_fold_df,
  aes(x = env_mean_predictor_sd, y = r2, color = model)
) +
  geom_point(alpha = 0.3, size = 1.6) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.9) +
  facet_grid(cv_type ~ n_folds, scales = "free", labeller = labeller(cv_type = function(x) ifelse(x == "pixel_grouped", "", x))) +
  scale_color_discrete(name = "Model") +
  labs(
    title = ("B. Within-fold predictor spread vs fold R2"),
    x = "Mean predictor SD within fold (scaled predictors)",
    y = ("Fold R2")
  )

ggsave(
  filename = file.path(sens_dir, "sensitivity_env_predictor_sd_vs_r2.png"),
  plot = p_env_predictor_sd_vs_r2, width = 13, height = 8, dpi = 220
)

combine_env_r2_plots(
  pairdist_plot = p_env_pairdist_vs_r2,
  predictor_sd_plot = p_env_predictor_sd_vs_r2,
  out_file = file.path(sens_dir, "sensitivity_env_r2_combined.png"),
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
  scale_color_discrete(name = "Model") +
  labs(
    title = ("Aggregate environmental-spread effect on R2"),
    subtitle = "Correlation between fold environmental pairwise distance and fold R2 (averaged across seeds)",
    x = "Number of folds",
    y = "Correlation between fold environmental pairwise distance and fold R2"
  )

ggsave(
  filename = file.path(sens_dir, "sensitivity_env_effect_summary.png"),
  plot = p_env_effect, width = 11, height = 6, dpi = 220
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
      group = interaction(model, quantile),
      color = model)
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
    title = "Distribution of fold environmental spread across k",
    subtitle = "Percentile curves shown from the 1st to 99th percentile",
    x = "Number of folds",
    y = "Fold mean pairwise environmental distance",
    color = "Model"
  )

ggsave(
  filename = file.path(sens_dir, "sensitivity_env_distribution_percentiles.png"),
  plot = p_env_dist, width = 11, height = 6, dpi = 300
)
}

# 8) Optional: tuning seed-count sweep performance plots
if (!is.null(tuning_seed_sweep_df) && nrow(tuning_seed_sweep_df) > 0L &&
    !is.na(tuning_seed_sweep_path) && nzchar(tuning_seed_sweep_path)) {
  source(file.path(project_root, "modelling/R/plot_tuning_seed_sweep.R"))
  sweep_plot_dir <- dirname(tuning_seed_sweep_path)
  dir.create(sweep_plot_dir, recursive = TRUE, showWarnings = FALSE)
  plot_tuning_seed_sweep_summary(tuning_seed_sweep_df, sweep_plot_dir)
} else {
  message("No tuning seed sweep summary found; skipping sweep R2/RMSE plots.")
}

cat("Wrote sensitivity plots to:\n", sens_dir, "\n", sep = "")