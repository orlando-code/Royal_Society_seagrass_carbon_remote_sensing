# Plots for tuning seed-count sweep summaries (shared by tuning_seed_sweep.R and sensitivity_plots.R).

#' @param sweep_df Rows from sensitivity_tuning_seed_sweep_summary.csv (or equivalent tibble).
#' @param out_dir Directory for PNG outputs (e.g. sensitivity_suite/).
#' @param project_root Repo root; defaults to here::here().
#' @param dpi Resolution for PNG outputs.
plot_tuning_seed_sweep_summary <- function(sweep_df, out_dir,
                                           project_root = here::here(),
                                           dpi = 220) {
  if (is.null(sweep_df) || nrow(sweep_df) == 0L) {
    message("No tuning seed sweep rows; skipping R2/RMSE sweep plots.")
    return(invisible(NULL))
  }
  pc_env <- new.env(parent = baseenv())
  sys.source(
    file.path(project_root, "modelling/plots/plot_config.R"),
    envir = pc_env
  )
  METRIC_LINESTYLES <- pc_env$METRIC_LINESTYLES
  MODEL_COLOURS <- pc_env$MODEL_COLOURS

  ggplot2::theme_set(pc_env$theme_paper())

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
    ) |>
    dplyr::mutate(
      pooling = dplyr::case_when(
        .data$metric_type %in% c("Pooled R2", "Pooled RMSE") ~ "Pooled",
        .data$metric_type %in% c("Mean fold R2", "Mean fold RMSE") ~ "Mean fold",
        TRUE ~ NA_character_
      ),
      pooling = factor(.data$pooling, levels = c("Pooled", "Mean fold"))
    )

  # Shared linetype levels (Pooled vs Mean fold) so patchwork can collect one guide on the combined plot.
  scale_pooling_linetype <- function() {
    ggplot2::scale_linetype_manual(
      name = "Type",
      values = c(
        "Pooled" = METRIC_LINESTYLES[["Pooled R2"]],
        "Mean fold" = METRIC_LINESTYLES[["Mean-fold R2"]]
      )
    )
  }

  r2_stats <- sweep_stats |> dplyr::filter(metric_type %in% c("Mean fold R2", "Pooled R2"))
  p_r2 <- ggplot2::ggplot(
    r2_stats,
    ggplot2::aes(x = n_tuning_seeds, y = mean_value, color = model, linetype = pooling)
  ) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = mean_value - se_value, ymax = mean_value + se_value),
      width = 0.18,
      alpha = 0.7,
      na.rm = TRUE,
      show.legend = FALSE
    ) +
    ggplot2::scale_color_manual(values = MODEL_COLOURS, name = "Model", drop = FALSE) +
    scale_pooling_linetype() +
    ggplot2::labs(
      # title = expression("Performance vs number of tuning seeds (R"^2*" values)"),
      x = "Number of model tuning seeds",
      y = expression("R"^2*" value"),
      color = "Model"
    )

  ggplot2::ggsave(
    filename = file.path(out_dir, "sensitivity_tuning_seed_sweep_r2.png"),
    plot = p_r2,
    width = 11,
    height = 6,
    dpi = dpi
  )

  rmse_stats <- sweep_stats |> dplyr::filter(metric_type %in% c("Mean fold RMSE", "Pooled RMSE"))
  p_rmse <- ggplot2::ggplot(
    rmse_stats,
    ggplot2::aes(x = n_tuning_seeds, y = mean_value, color = model, linetype = pooling)
  ) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = mean_value - se_value, ymax = mean_value + se_value),
      width = 0.18,
      alpha = 0.7,
      na.rm = TRUE,
      show.legend = FALSE
    ) +
    ggplot2::scale_color_manual(values = MODEL_COLOURS, name = "Model", drop = FALSE) +
    scale_pooling_linetype() +
    ggplot2::labs(
      # title = "Performance vs number of tuning seeds (RMSE)",
      x = "Number of model tuning seeds",
      y = "RMSE",
      color = "Model"
    )

  ggplot2::ggsave(
    filename = file.path(out_dir, "sensitivity_tuning_seed_sweep_rmse.png"),
    plot = p_rmse,
    width = 11,
    height = 6,
    dpi = dpi
  )

  combined_path <- NA_character_
  if (requireNamespace("patchwork", quietly = TRUE)) {
    # Remove x-axis label and ticks from the top plot (p_r2)
    p_r2_shared <- p_r2 +
      ggplot2::theme(
        axis.title.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank()
      )
    # p_rmse (bottom plot) keeps the x-axis
    p_combined <- p_r2_shared / p_rmse +
      patchwork::plot_layout(guides = "collect", heights = c(1, 1))
    combined_path <- file.path(out_dir, "sensitivity_tuning_seed_sweep_rmse_r2_combined.png")
    ggplot2::ggsave(
      filename = combined_path,
      plot = p_combined,
      width = 11,
      height = 7,
      dpi = dpi
    )
  } else {
    message("Install package 'patchwork' to write the combined RMSE/R2 tuning sweep figure.")
  }

  cat(
    "Wrote tuning seed sweep plots to:\n  ",
    file.path(out_dir, "sensitivity_tuning_seed_sweep_r2.png"), "\n  ",
    file.path(out_dir, "sensitivity_tuning_seed_sweep_rmse.png"), "\n",
    sep = ""
  )
  if (!is.na(combined_path)) {
    cat("  ", combined_path, "\n", sep = "")
  }
  invisible(NULL)
}

#' Pick a representative tuning_seed_list at fixed n_tuning_seeds.
#' Representative is defined as the subset whose mean RMSE across models is
#' closest to the median mean RMSE across candidate subsets at that n.
#' Robust score is still reported as mean(mean RMSE + lambda * sd RMSE).
#'
#' @return A list with integer vector `robust_fold_seed_list`, `tuning_seed_list` string,
#'   and diagnostic columns suitable for one-row CSV.
pick_representative_tuning_seed_subset <- function(sweep_df,
                                                  n_tuning_seeds = 5L,
                                                  metric_col = "mean_mean_rmse",
                                                  metric_sd_col = "sd_mean_rmse") {
  n_seeds_target <- as.integer(n_tuning_seeds)[1L]
  if (!metric_col %in% names(sweep_df)) {
    stop("pick_representative_tuning_seed_subset: missing column ", metric_col)
  }
  if (!metric_sd_col %in% names(sweep_df)) {
    stop("pick_representative_tuning_seed_subset: missing column ", metric_sd_col)
  }
  sub <- sweep_df |> dplyr::filter(.data$n_tuning_seeds == .env$n_seeds_target)
  if (nrow(sub) == 0L) {
    return(NULL)
  }
  by_list <- sub |>
    dplyr::group_by(.data$tuning_seed_list) |>
    dplyr::summarise(
      mean_metric_across_models = mean(.data[[metric_col]], na.rm = TRUE),
      mean_sd_metric_across_models = mean(.data[[metric_sd_col]], na.rm = TRUE),
      .groups = "drop"
    )
  median_metric <- stats::median(by_list$mean_metric_across_models, na.rm = TRUE)
  by_list$abs_deviation_from_median <- abs(by_list$mean_metric_across_models - median_metric)
  cand <- by_list |>
    dplyr::filter(.data$abs_deviation_from_median == min(.data$abs_deviation_from_median, na.rm = TRUE)) |>
    dplyr::arrange(.data$mean_sd_metric_across_models, .data$tuning_seed_list)
  chosen <- cand[1L, ]
  seed_str <- as.character(chosen$tuning_seed_list)
  seed_ints <- as.integer(strsplit(seed_str, "-", fixed = TRUE)[[1L]])
  list(
    robust_fold_seed_list = seed_ints,
    tuning_seed_list = seed_str,
    n_tuning_seeds = n_seeds_target,
    mean_metric_across_models = as.numeric(chosen$mean_metric_across_models),
    mean_sd_metric_across_models = as.numeric(chosen$mean_sd_metric_across_models),
    selection_rule = "min abs(mean_mean_rmse_across_models - median(mean_mean_rmse_across_models)) across subsets at fixed n_tuning_seeds (tie-break: lower mean_sd_rmse, then seed string)"
  )
}
