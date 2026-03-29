# Plots for tuning seed-count sweep summaries (shared by tuning_seed_sweep.R and sensitivity_plots.R).

#' @param sweep_df Rows from sensitivity_tuning_seed_sweep_summary.csv (or equivalent tibble).
#' @param out_dir Directory for PNG outputs (e.g. sensitivity_suite/).
#' @param project_root Repo root; defaults to here::here().
plot_tuning_seed_sweep_summary <- function(sweep_df, out_dir,
                                           project_root = here::here()) {
  if (is.null(sweep_df) || nrow(sweep_df) == 0L) {
    message("No tuning seed sweep rows; skipping R2/RMSE sweep plots.")
    return(invisible(NULL))
  }
  pc_env <- new.env(parent = baseenv())
  sys.source(
    file.path(project_root, "modelling/R/plot_config.R"),
    envir = pc_env
  )
  METRIC_LINESTYLES <- pc_env$METRIC_LINESTYLES

  sweep_long <- sweep_df |>
    dplyr::select(
      dplyr::all_of(c(
        "model", "n_tuning_seeds", "mean_mean_r2", "mean_pooled_r2",
        "mean_mean_rmse", "mean_pooled_rmse"
      ))
    ) |>
    tidyr::pivot_longer(
      cols = c(mean_mean_r2, mean_pooled_r2, mean_mean_rmse, mean_pooled_rmse),
      names_to = "metric_type",
      values_to = "value"
    ) |>
    dplyr::mutate(
      metric_type = dplyr::recode(
        metric_type,
        mean_mean_r2 = "Mean fold R2",
        mean_pooled_r2 = "Pooled R2",
        mean_mean_rmse = "Mean fold RMSE",
        mean_pooled_rmse = "Pooled RMSE"
      )
    )

  sweep_stats <- sweep_long |>
    dplyr::group_by(model, n_tuning_seeds, metric_type) |>
    dplyr::summarise(
      mean_value = mean(value, na.rm = TRUE),
      sd_value = sd(value, na.rm = TRUE),
      n_subsets = sum(is.finite(value)),
      se_value = ifelse(n_subsets > 1L, sd_value / sqrt(n_subsets), NA_real_),
      .groups = "drop"
    )

  r2_stats <- sweep_stats |> dplyr::filter(metric_type %in% c("Mean fold R2", "Pooled R2"))
  p_r2 <- ggplot2::ggplot(
    r2_stats,
    ggplot2::aes(x = n_tuning_seeds, y = mean_value, color = model, linetype = metric_type)
  ) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = mean_value - se_value, ymax = mean_value + se_value),
      width = 0.18,
      alpha = 0.7,
      na.rm = TRUE
    ) +
    ggplot2::scale_linetype_manual(
      values = c(
        "Pooled R2" = METRIC_LINESTYLES[["Pooled R2"]],
        "Mean fold R2" = METRIC_LINESTYLES[["Mean-fold R2"]]
      )
    ) +
    ggplot2::labs(
      title = expression("Performance vs number of tuning seeds (R"^2*" values)"),
      x = "Number of tuning seeds (robust_fold_seed_list size)",
      y = expression("R"^2*" value"),
      color = "Model",
      linetype = "Metric"
    )

  ggplot2::ggsave(
    filename = file.path(out_dir, "sensitivity_tuning_seed_sweep_r2.png"),
    plot = p_r2,
    width = 11,
    height = 6,
    dpi = 220
  )

  rmse_stats <- sweep_stats |> dplyr::filter(metric_type %in% c("Mean fold RMSE", "Pooled RMSE"))
  p_rmse <- ggplot2::ggplot(
    rmse_stats,
    ggplot2::aes(x = n_tuning_seeds, y = mean_value, color = model, linetype = metric_type)
  ) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = mean_value - se_value, ymax = mean_value + se_value),
      width = 0.18,
      alpha = 0.7,
      na.rm = TRUE
    ) +
    ggplot2::scale_linetype_manual(
      values = c(
        "Pooled RMSE" = METRIC_LINESTYLES[["Pooled RMSE"]],
        "Mean fold RMSE" = METRIC_LINESTYLES[["Mean-fold RMSE"]]
      )
    ) +
    ggplot2::labs(
      title = "Performance vs number of tuning seeds (RMSE)",
      x = "Number of tuning seeds (robust_fold_seed_list size)",
      y = "RMSE",
      color = "Model",
      linetype = "Metric"
    )

  ggplot2::ggsave(
    filename = file.path(out_dir, "sensitivity_tuning_seed_sweep_rmse.png"),
    plot = p_rmse,
    width = 11,
    height = 6,
    dpi = 220
  )

  cat(
    "Wrote tuning seed sweep plots to:\n  ",
    file.path(out_dir, "sensitivity_tuning_seed_sweep_r2.png"), "\n  ",
    file.path(out_dir, "sensitivity_tuning_seed_sweep_rmse.png"), "\n",
    sep = ""
  )
  invisible(NULL)
}

#' Pick the tuning_seed_list with best robust RMSE score at fixed n_tuning_seeds.
#' Robust score is mean(mean RMSE + lambda * sd RMSE) across models.
#'
#' @return A list with integer vector `robust_fold_seed_list`, `tuning_seed_list` string,
#'   and diagnostic columns suitable for one-row CSV.
pick_representative_tuning_seed_subset <- function(sweep_df,
                                                  n_tuning_seeds = 5L,
                                                  metric_col = "mean_mean_rmse",
                                                  metric_sd_col = "sd_mean_rmse",
                                                  rmse_lambda = 0.5) {
  n_seeds_target <- as.integer(n_tuning_seeds)[1L]
  if (!metric_col %in% names(sweep_df)) {
    stop("pick_representative_tuning_seed_subset: missing column ", metric_col)
  }
  if (!metric_sd_col %in% names(sweep_df)) {
    stop("pick_representative_tuning_seed_subset: missing column ", metric_sd_col)
  }
  rmse_lambda <- max(0, as.numeric(rmse_lambda))
  sub <- sweep_df |> dplyr::filter(.data$n_tuning_seeds == .env$n_seeds_target)
  if (nrow(sub) == 0L) {
    return(NULL)
  }
  by_list <- sub |>
    dplyr::group_by(.data$tuning_seed_list) |>
    dplyr::summarise(
      mean_metric_across_models = mean(.data[[metric_col]], na.rm = TRUE),
      mean_sd_metric_across_models = mean(.data[[metric_sd_col]], na.rm = TRUE),
      robust_rmse_score_across_models = mean(.data[[metric_col]] + .env$rmse_lambda * .data[[metric_sd_col]], na.rm = TRUE),
      .groups = "drop"
    )
  min_score <- min(by_list$robust_rmse_score_across_models, na.rm = TRUE)
  cand <- by_list |>
    dplyr::filter(.data$robust_rmse_score_across_models == min_score) |>
    dplyr::arrange(.data$mean_metric_across_models, .data$tuning_seed_list)
  chosen <- cand[1L, ]
  seed_str <- as.character(chosen$tuning_seed_list)
  seed_ints <- as.integer(strsplit(seed_str, "-", fixed = TRUE)[[1L]])
  list(
    robust_fold_seed_list = seed_ints,
    tuning_seed_list = seed_str,
    n_tuning_seeds = n_seeds_target,
    mean_metric_across_models = as.numeric(chosen$mean_metric_across_models),
    mean_sd_metric_across_models = as.numeric(chosen$mean_sd_metric_across_models),
    robust_rmse_score_across_models = as.numeric(chosen$robust_rmse_score_across_models),
    robust_rmse_lambda = rmse_lambda,
    selection_rule = "min mean(mean_mean_rmse + lambda * sd_mean_rmse) across models at fixed n_tuning_seeds"
  )
}