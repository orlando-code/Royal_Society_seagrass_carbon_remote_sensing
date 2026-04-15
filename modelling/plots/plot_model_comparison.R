## Plot helpers for model_comparison.R outputs.
##
## Sources at runtime:
## - modelling/R/init_repo.R (bootstrap; via seagrass_init_repo)
## - modelling/pipeline_config.R (current model_list/cv_output_dir)
## - modelling/plots/plot_config.R (model colours/theme constants)
##
## Required input CSVs:
## - <analysis_out_dir>/model_comparison_summary.csv
## - <analysis_out_dir>/model_comparison_pooled_summary.csv
## - <analysis_out_dir>/model_comparison_by_fold.csv
##
## Written PNGs:
## - <analysis_out_dir>/model_comparison_mean_rmse_barplot.png
## - <analysis_out_dir>/model_comparison_mean_r2_barplot.png
## - <analysis_out_dir>/model_comparison_combined_plots.png

MODEL_COMPARISON_REQUIRED_CSVS <- c(
  "model_comparison_summary.csv",
  "model_comparison_pooled_summary.csv",
  "model_comparison_by_fold.csv"
)

MODEL_COMPARISON_OUTPUT_PNGS <- c(
  "model_comparison_mean_rmse_barplot.png",
  "model_comparison_mean_r2_barplot.png",
  "model_comparison_combined_plots.png"
)

plot_log_mc <- function(message_text) {
  message(sprintf("[plot:plot_model_comparison] %s", message_text))
}

rel_path_mc <- function(path) {
  p <- normalizePath(path, winslash = "/", mustWork = FALSE)
  root <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  sub(paste0("^", root, "/?"), "", p)
}

bootstrap_plot_model_comparison <- function() {
  if (!exists("seagrass_init_repo", mode = "function", inherits = TRUE)) {
    init_path <- file.path("modelling", "R", "init_repo.R")
    if (!file.exists(init_path)) {
      ff <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
      if (!length(ff)) stop("Run from repo root or with: Rscript /path/to/plot_model_comparison.R", call. = FALSE)
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
  list(
    project_root = project_root,
    cfg = get_pipeline_config()
  )
}

plot_model_comparison_outputs <- function(summary_tbl,
                                          pooled_summary_tbl,
                                          by_fold,
                                          model_list,
                                          out_dir,
                                          dpi = 300,
                                          model_colours = NULL) {
  collapse_approach_plot <- function(df, mean_col, sd_col) {
    df %>%
      dplyr::mutate(
        approach_plot = dplyr::case_when(
          approach == "species_mean_baseline" ~ "Species\nMean Baseline",
          TRUE ~ model
        )
      ) %>%
      dplyr::group_by(approach_plot) %>%
      dplyr::summarise(
        mean_val = mean(.data[[mean_col]], na.rm = TRUE),
        sd_val = mean(.data[[sd_col]], na.rm = TRUE),
        .groups = "drop"
      )
  }

  rmse_summary_fold <- summary_tbl %>%
    dplyr::select(model, approach, mean_rmse, sd_rmse) %>%
    collapse_approach_plot("mean_rmse", "sd_rmse") %>%
    dplyr::rename(mean_rmse = mean_val, sd_rmse = sd_val) %>%
    dplyr::mutate(aggregation = "Mean fold")

  rmse_summary_pooled <- pooled_summary_tbl %>%
    dplyr::select(model, approach, mean_rmse, sd_rmse) %>%
    collapse_approach_plot("mean_rmse", "sd_rmse") %>%
    dplyr::rename(mean_rmse = mean_val, sd_rmse = sd_val) %>%
    dplyr::mutate(aggregation = "Pooled")

  rmse_order <- rmse_summary_fold %>%
    dplyr::arrange(mean_rmse) %>%
    dplyr::pull(approach_plot) %>%
    unique()

  rmse_both <- dplyr::bind_rows(rmse_summary_fold, rmse_summary_pooled) %>%
    dplyr::mutate(
      approach_plot = factor(.data$approach_plot, levels = rmse_order),
      aggregation = factor(.data$aggregation, levels = c("Mean fold", "Pooled"))
    )

  r2_summary_fold <- summary_tbl %>%
    dplyr::select(model, approach, mean_r2, sd_r2) %>%
    collapse_approach_plot("mean_r2", "sd_r2") %>%
    dplyr::rename(mean_r2 = mean_val, sd_r2 = sd_val) %>%
    dplyr::mutate(aggregation = "Mean fold")

  r2_summary_pooled <- pooled_summary_tbl %>%
    dplyr::select(model, approach, mean_r2, sd_r2) %>%
    collapse_approach_plot("mean_r2", "sd_r2") %>%
    dplyr::rename(mean_r2 = mean_val, sd_r2 = sd_val) %>%
    dplyr::mutate(aggregation = "Pooled")

  r2_both <- dplyr::bind_rows(r2_summary_fold, r2_summary_pooled) %>%
    dplyr::mutate(
      approach_plot = factor(.data$approach_plot, levels = rmse_order),
      aggregation = factor(.data$aggregation, levels = c("Mean fold", "Pooled"))
    )

  if (is.null(model_colours)) {
    model_colours <- if (exists("MODEL_COLOURS", inherits = TRUE)) get("MODEL_COLOURS", inherits = TRUE) else list()
  }
  plot_colours <- c("Species\nMean Baseline" = "#8c8c8c")
  for (m in model_list) {
    col_m <- model_colours[[m]]
    if (is.null(col_m) || !nzchar(as.character(col_m))) col_m <- "#333333"
    plot_colours[m] <- as.character(col_m)
  }

  dodge_pos <- ggplot2::position_dodge(width = 0.88)

  p_rmse_bar <- ggplot2::ggplot(
    rmse_both,
    ggplot2::aes(x = mean_rmse, y = approach_plot, fill = approach_plot, group = aggregation)
  ) +
    ggplot2::geom_col(
      ggplot2::aes(alpha = aggregation, colour = aggregation),
      position = dodge_pos,
      width = 0.38,
      linewidth = 0.55
    ) +
    ggplot2::geom_errorbarh(
      ggplot2::aes(
        xmin = mean_rmse - sd_rmse,
        xmax = mean_rmse + sd_rmse,
        group = aggregation
      ),
      position = dodge_pos,
      width = 0.18,
      linewidth = 0.35,
      colour = "#444444"
    ) +
    ggplot2::scale_fill_manual(values = plot_colours, drop = FALSE) +
    ggplot2::scale_alpha_manual(
      values = c("Mean fold" = 1, "Pooled" = 0.55),
      name = "CV summary",
      labels = c(
        "Mean fold" = "Mean of fold-wise metrics",
        "Pooled" = "Pooled (all held-out points across all seeds)"
      )
    ) +
    ggplot2::scale_colour_manual(
      values = c("Mean fold" = "#1a1a1a", "Pooled" = "#6a6a6a"),
      guide = "none"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    ggplot2::labs(
      x = "RMSE (mean ± 1 SD)",
      y = "Model"
    ) +
    ggplot2::guides(fill = "none") +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold"),
      legend.position = "bottom"
    )

  ggplot2::ggsave(
    file.path(out_dir, "model_comparison_mean_rmse_barplot.png"),
    p_rmse_bar, width = 9, height = 4.8, dpi = dpi
  )
  plot_log_mc(paste0("Wrote: ", rel_path_mc(file.path(out_dir, "model_comparison_mean_rmse_barplot.png"))))

  p_r2_bar <- ggplot2::ggplot(
    r2_both,
    ggplot2::aes(x = mean_r2, y = approach_plot, fill = approach_plot, group = aggregation)
  ) +
    ggplot2::geom_col(
      ggplot2::aes(alpha = aggregation, colour = aggregation),
      position = dodge_pos,
      width = 0.38,
      linewidth = 0.55
    ) +
    ggplot2::geom_errorbarh(
      ggplot2::aes(
        xmin = mean_r2 - sd_r2,
        xmax = mean_r2 + sd_r2,
        group = aggregation
      ),
      position = dodge_pos,
      width = 0.18,
      linewidth = 0.35,
      colour = "#444444"
    ) +
    ggplot2::scale_fill_manual(values = plot_colours, drop = FALSE) +
    ggplot2::scale_alpha_manual(
      values = c("Mean fold" = 1, "Pooled" = 0.55),
      name = "CV summary",
      labels = c(
        "Mean fold" = "Mean of fold-wise metrics",
        "Pooled" = "Pooled (all held-out points across all seeds)"
      )
    ) +
    ggplot2::scale_colour_manual(
      values = c("Mean fold" = "#1a1a1a", "Pooled" = "#6a6a6a"),
      guide = "none"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = expression(paste(R^2, " (mean ± 1 SD)")),
      y = "Model"
    ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    ggplot2::guides(fill = "none") +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold"),
      legend.position = "bottom"
    )

  ggplot2::ggsave(
    file.path(out_dir, "model_comparison_mean_r2_barplot.png"),
    p_r2_bar, width = 9, height = 4.8, dpi = dpi
  )
  plot_log_mc(paste0("Wrote: ", rel_path_mc(file.path(out_dir, "model_comparison_mean_r2_barplot.png"))))

  p_combined <- p_rmse_bar + p_r2_bar + patchwork::plot_layout(nrow = 2, guides = "collect") &
    ggplot2::theme(legend.position = "bottom")
  ggplot2::ggsave(
    file.path(out_dir, "model_comparison_combined_plots.png"),
    p_combined, width = 9, height = 7.2, dpi = dpi
  )
  plot_log_mc(paste0("Wrote: ", rel_path_mc(file.path(out_dir, "model_comparison_combined_plots.png"))))

  invisible(list(
    p_rmse_bar = p_rmse_bar,
    p_r2_bar = p_r2_bar,
    p_combined = p_combined
  ))
}

plot_model_comparison_from_csv <- function(out_dir,
                                           model_list = get0("model_list", envir = .GlobalEnv, ifnotfound = character()),
                                           dpi = 300L) {
  out_dir <- as.character(out_dir)
  dpi <- as.integer(dpi)
  if (!is.finite(dpi) || dpi <= 0L) dpi <- 300L
  if (!dir.exists(out_dir)) {
    stop("analysis_out_dir does not exist: ", out_dir, call. = FALSE)
  }
  summary_csv <- file.path(out_dir, "model_comparison_summary.csv")
  pooled_summary_csv <- file.path(out_dir, "model_comparison_pooled_summary.csv")
  by_fold_csv <- file.path(out_dir, "model_comparison_by_fold.csv")
  if (!file.exists(summary_csv) || !file.exists(pooled_summary_csv) || !file.exists(by_fold_csv)) {
    stop(
      "Missing required CSVs in analysis_out_dir. Expected:\n  ",
      paste(file.path("<analysis_out_dir>", MODEL_COMPARISON_REQUIRED_CSVS), collapse = "\n  "),
      call. = FALSE
    )
  }

  summary_tbl <- readr::read_csv(summary_csv, show_col_types = FALSE)
  pooled_summary_tbl <- readr::read_csv(pooled_summary_csv, show_col_types = FALSE)
  by_fold <- readr::read_csv(by_fold_csv, show_col_types = FALSE)

  requested_models <- unique(as.character(model_list))
  requested_models <- requested_models[nzchar(requested_models)]
  if (length(requested_models) == 0L) {
    stop("model_list from current config is empty; nothing to plot.", call. = FALSE)
  }
  available_models <- summary_tbl %>%
    dplyr::filter(.data$approach != "species_mean_baseline") %>%
    dplyr::pull(.data$model) %>%
    as.character() %>%
    unique()
  missing_requested <- setdiff(requested_models, available_models)
  if (length(missing_requested) > 0L) {
    stop(
      "Requested model_list contains models not present in model comparison CSVs.\n",
      "Requested: ", paste(requested_models, collapse = ", "), "\n",
      "Available: ", paste(available_models, collapse = ", "), "\n",
      "Missing: ", paste(missing_requested, collapse = ", "),
      call. = FALSE
    )
  }
  extra_available <- setdiff(available_models, requested_models)
  if (length(extra_available) > 0L) {
    plot_log_mc(paste0("Input has extra models that will be excluded: ", paste(extra_available, collapse = ", ")))
  }
  summary_tbl <- summary_tbl[summary_tbl$model %in% requested_models, , drop = FALSE]
  pooled_summary_tbl <- pooled_summary_tbl[pooled_summary_tbl$model %in% requested_models, , drop = FALSE]
  by_fold <- by_fold[by_fold$model %in% requested_models, , drop = FALSE]

  plot_model_comparison_outputs(
    summary_tbl = summary_tbl,
    pooled_summary_tbl = pooled_summary_tbl,
    by_fold = by_fold,
    model_list = requested_models,
    out_dir = out_dir,
    dpi = dpi
  )
  invisible(list(
    summary_tbl = summary_tbl,
    pooled_summary_tbl = pooled_summary_tbl,
    by_fold = by_fold
  ))
}

resolve_model_comparison_out_dir <- function(project_root, cfg, run_tag = NULL) {
  robust_seed_str <- paste(as.integer(cfg$robust_fold_seed_list), collapse = "-")
  preferred_seed_root <- file.path(project_root, "output", paste0(cfg$cv_regime_name, "_", robust_seed_str))

  run_output_dir <- get0("run_output_dir", envir = .GlobalEnv, ifnotfound = NA_character_)
  run_output_dir <- as.character(run_output_dir)
  if (length(run_output_dir) != 1L || !nzchar(run_output_dir) || is.na(run_output_dir)) {
    run_output_dir <- NA_character_
  }
  if (!is.na(run_output_dir) && !dir.exists(run_output_dir)) {
    run_output_dir <- NA_character_
  }

  # Prefer the seed-derived directory from current config when it exists.
  analysis_root <- if (dir.exists(preferred_seed_root)) preferred_seed_root else run_output_dir
  if (is.na(analysis_root)) {
    cv_output_dir <- as.character(cfg$cv_output_dir)
    if (!nzchar(cv_output_dir)) {
      stop("cv_output_dir missing in current config.", call. = FALSE)
    }
    candidate <- file.path(project_root, cv_output_dir)
    if (dir.exists(candidate)) analysis_root <- candidate
  }
  if (is.na(analysis_root)) {
    stop("Could not resolve analysis root from run_output_dir/cv_output_dir/current robust seeds.", call. = FALSE)
  }

  model_comparison_eval_run_dir <- get0("model_comparison_eval_run_dir", envir = .GlobalEnv, ifnotfound = NA_character_)
  model_comparison_eval_run_dir <- as.character(model_comparison_eval_run_dir)
  if (length(model_comparison_eval_run_dir) != 1L || !nzchar(model_comparison_eval_run_dir) || is.na(model_comparison_eval_run_dir)) {
    model_comparison_eval_run_dir <- NA_character_
  }
  if (is.null(run_tag) || !nzchar(as.character(run_tag))) {
    run_tag <- if (!is.na(model_comparison_eval_run_dir)) basename(normalizePath(model_comparison_eval_run_dir, winslash = "/", mustWork = FALSE)) else "current_config"
  }
  file.path(analysis_root, "analysis", as.character(run_tag))
}

plot_model_comparison_from_current_config <- function(dpi = NULL, run_tag = NULL) {
  boot <- bootstrap_plot_model_comparison()
  project_root <- boot$project_root
  cfg <- boot$cfg
  dpi_use <- if (is.null(dpi)) as.integer(cfg$dpi) else as.integer(dpi)
  if (!is.finite(dpi_use) || dpi_use <= 0L) dpi_use <- 300L
  model_list <- as.character(cfg$model_list)
  model_list <- model_list[nzchar(model_list)]
  if (length(model_list) == 0L) {
    stop("model_list from current config is empty; nothing to plot.", call. = FALSE)
  }
  out_dir <- resolve_model_comparison_out_dir(project_root, cfg, run_tag = run_tag)
  plot_log_mc("Sources: modelling/R/init_repo.R, modelling/pipeline_config.R, modelling/plots/plot_config.R")
  plot_log_mc(paste0("Input directory: ", rel_path_mc(out_dir)))
  plot_log_mc(paste0("Requested model_list: ", paste(model_list, collapse = ", ")))
  plot_log_mc(paste0("Output directory: ", rel_path_mc(out_dir)))
  plot_model_comparison_from_csv(out_dir = out_dir, model_list = model_list, dpi = dpi_use)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 1L && nzchar(args[[1L]])) {
    out_dir_arg <- args[[1L]]
    dpi_arg <- if (length(args) >= 2L && nzchar(args[[2L]])) suppressWarnings(as.integer(args[[2L]])) else NA_integer_
    boot <- bootstrap_plot_model_comparison()
    cfg <- boot$cfg
    model_list <- as.character(cfg$model_list)
    model_list <- model_list[nzchar(model_list)]
    if (!grepl("^/", out_dir_arg)) {
      out_dir_arg <- file.path(boot$project_root, out_dir_arg)
    }
    plot_log_mc("Sources: modelling/R/init_repo.R, modelling/pipeline_config.R, modelling/plots/plot_config.R")
    plot_log_mc(paste0("Input directory: ", rel_path_mc(out_dir_arg)))
    plot_log_mc(paste0("Requested model_list: ", paste(model_list, collapse = ", ")))
    plot_log_mc(paste0("Output directory: ", rel_path_mc(out_dir_arg)))
    plot_model_comparison_from_csv(
      out_dir = out_dir_arg,
      model_list = model_list,
      dpi = if (is.finite(dpi_arg) && dpi_arg > 0L) dpi_arg else as.integer(cfg$dpi)
    )
  } else {
    plot_model_comparison_from_current_config()
  }
} else {
  disable_autorun <- isTRUE(get0("plot_model_comparison_disable_autorun_on_source", envir = .GlobalEnv, ifnotfound = FALSE))
  if (interactive() && !disable_autorun) {
    plot_log_mc("Sourced in interactive session; auto-running plot generation from current config.")
    plot_model_comparison_from_current_config()
  } else {
    resolved_out_dir <- tryCatch({
      boot <- bootstrap_plot_model_comparison()
      resolve_model_comparison_out_dir(boot$project_root, boot$cfg)
    }, error = function(e) NA_character_)
    resolved_out_dir_display <- if (!is.na(resolved_out_dir) && nzchar(resolved_out_dir)) {
      rel_path_mc(resolved_out_dir)
    } else {
      "<could not resolve from current config>"
    }
    required_display <- if (!is.na(resolved_out_dir) && nzchar(resolved_out_dir)) {
      paste(file.path(resolved_out_dir_display, MODEL_COMPARISON_REQUIRED_CSVS), collapse = "\n  ")
    } else {
      paste(file.path("<analysis_out_dir>", MODEL_COMPARISON_REQUIRED_CSVS), collapse = "\n  ")
    }
    written_display <- if (!is.na(resolved_out_dir) && nzchar(resolved_out_dir)) {
      paste(file.path(resolved_out_dir_display, MODEL_COMPARISON_OUTPUT_PNGS), collapse = "\n  ")
    } else {
      paste(file.path("<analysis_out_dir>", MODEL_COMPARISON_OUTPUT_PNGS), collapse = "\n  ")
    }
    message(
      "[plot:plot_model_comparison] Loaded.\n",
      "Sources: modelling/R/init_repo.R, modelling/pipeline_config.R, modelling/plots/plot_config.R\n",
      "Resolved output directory from current config:\n  ",
      resolved_out_dir_display, "\n",
      "Required CSVs in analysis output dir:\n  ",
      required_display,
      "\nWritten PNGs:\n  ",
      written_display,
      "\n",
      "Run interactively with:\n",
      "  plot_model_comparison_from_current_config()\n",
      "  # or: plot_model_comparison_from_csv('output/.../analysis/current_config')\n",
      "Or from shell with:\n",
      "  Rscript modelling/plots/plot_model_comparison.R\n",
      "  # optional override: Rscript modelling/plots/plot_model_comparison.R <analysis_out_dir> [dpi]"
    )
  }
}

