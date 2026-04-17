## Plot tuning-seed sweep outputs from a saved summary CSV.
##
## Sources at runtime:
## - modelling/R/init_repo.R (bootstrap; via seagrass_init_repo)
## - modelling/pipeline_config.R (current model_list/tuning sweep run id)
## - modelling/plots/plot_config.R (theme constants loaded by plot helper)
## - modelling/plots/plot_tuning_seed_sweep.R (plotting function definitions)
##
## Required input file:
## - <sweep_out_dir>/sensitivity_tuning_seed_sweep_summary.csv
##
## Usage:
##   source("modelling/plots/plot_tuning_seed_sweep.R")
##   plot_tuning_seed_sweep()  # auto-resolves from current config
##   # optional:
##   plot_tuning_seed_sweep("output/tuning_seed_sweep_runs/sweep_000")

TUNING_SWEEP_REQUIRED_CSV <- "sensitivity_tuning_seed_sweep_summary.csv"
TUNING_SWEEP_OUTPUT_PNGS <- c(
  "sensitivity_tuning_seed_sweep_r2.png",
  "sensitivity_tuning_seed_sweep_rmse.png",
  "sensitivity_tuning_seed_sweep_rmse_r2_combined.png"
)

plot_log_tss <- function(message_text) {
  message(sprintf("[plot:plot_tuning_seed_sweep] %s", message_text))
}

rel_path_tss <- function(path) {
  p <- normalizePath(path, winslash = "/", mustWork = FALSE)
  root <- normalizePath(project_root, winslash = "/", mustWork = FALSE)
  sub(paste0("^", root, "/?"), "", p)
}

if (!exists("seagrass_init_repo", mode = "function", inherits = TRUE)) {
  init_path <- file.path("modelling", "R", "init_repo.R")
  if (!file.exists(init_path)) {
    ff <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
    if (!length(ff)) stop("Run from repo root or with: Rscript /path/to/plot_tuning_seed_sweep.R", call. = FALSE)
    script_path <- normalizePath(sub("^--file=", "", ff[[1]]), winslash = "/", mustWork = FALSE)
    init_path <- normalizePath(file.path(dirname(script_path), "..", "R", "init_repo.R"), winslash = "/", mustWork = FALSE)
  }
  if (!file.exists(init_path)) stop("Missing bootstrap helper: modelling/R/init_repo.R", call. = FALSE)
  sys.source(init_path, envir = .GlobalEnv)
}

project_root <- seagrass_init_repo(
  packages = c("here", "readr", "dplyr", "tidyr", "ggplot2", "patchwork"),
  source_files = c("modelling/pipeline_config.R", "modelling/plots/plot_config.R")
)
cfg <- get_pipeline_config()

resolve_sweep_out_dir <- function(project_root, cfg, sweep_out_dir = NULL) {
  if (!is.null(sweep_out_dir) && nzchar(as.character(sweep_out_dir))) {
    p <- as.character(sweep_out_dir)
    if (!grepl("^/", p)) p <- file.path(project_root, p)
    return(p)
  }
  run_id <- as.character(cfg$tuning_seed_sweep_run_id)
  if (!is.na(run_id) && nzchar(run_id)) {
    candidate <- file.path(project_root, "output", "tuning_seed_sweep_runs", paste0("sweep_", run_id))
    if (dir.exists(candidate)) return(candidate)
  }
  sweep_dirs <- Sys.glob(file.path(project_root, "output", "tuning_seed_sweep_runs", "sweep_*"))
  sweep_dirs <- sweep_dirs[dir.exists(sweep_dirs)]
  if (length(sweep_dirs) == 0L) {
    stop("No sweep_* directory found under output/tuning_seed_sweep_runs.", call. = FALSE)
  }
  sweep_dirs[which.max(file.info(sweep_dirs)$mtime)]
}

validate_tss_model_list <- function(sweep_df, model_list) {
  available <- sort(unique(as.character(sweep_df$model)))
  requested <- as.character(model_list)
  missing <- setdiff(requested, available)
  if (length(missing) > 0L) {
    stop(
      "Requested model_list not found in sweep summary CSV: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  requested
}

plot_tuning_seed_sweep_summary <- function(
  sweep_df,
  out_dir,
  project_root,
  dpi = cfg$dpi,
  model_list = cfg$model_list
) {
  required_cols <- c(
    "model", "n_tuning_seeds",
    "mean_pooled_r2", "mean_mean_r2",
    "mean_pooled_rmse", "mean_mean_rmse"
  )
  missing_cols <- setdiff(required_cols, names(sweep_df))
  if (length(missing_cols) > 0L) {
    stop(
      "Sweep summary missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  if (!("sd_pooled_r2" %in% names(sweep_df))) sweep_df$sd_pooled_r2 <- NA_real_
  if (!("sd_mean_r2" %in% names(sweep_df))) sweep_df$sd_mean_r2 <- NA_real_
  if (!("sd_pooled_rmse" %in% names(sweep_df))) sweep_df$sd_pooled_rmse <- NA_real_
  if (!("sd_mean_rmse" %in% names(sweep_df))) sweep_df$sd_mean_rmse <- NA_real_

  model_list <- validate_tss_model_list(sweep_df, model_list)

  df <- sweep_df %>%
    dplyr::mutate(
      model = as.character(model),
      n_tuning_seeds = as.integer(n_tuning_seeds)
    ) %>%
    dplyr::filter(model %in% model_list)

  # Aggregate across sweep repeats so each model has one pooled and one mean-fold
  # point per n_tuning_seeds; error bars are standard errors across repeats.
  summarise_by_repeat <- function(long_df) {
    long_df %>%
      dplyr::group_by(model, n_tuning_seeds, Type) %>%
      dplyr::summarise(
        n_rep = sum(is.finite(Value)),
        Value = mean(Value, na.rm = TRUE),
        SE = dplyr::if_else(
          n_rep > 1L,
          stats::sd(Value_raw, na.rm = TRUE) / sqrt(n_rep),
          NA_real_
        ),
        .groups = "drop"
      )
  }

  r2_long <- dplyr::bind_rows(
    df %>%
      dplyr::transmute(
        model,
        n_tuning_seeds,
        Type = "Pooled",
        Value_raw = as.numeric(mean_pooled_r2)
      ),
    df %>%
      dplyr::transmute(
        model,
        n_tuning_seeds,
        Type = "Mean fold",
        Value_raw = as.numeric(mean_mean_r2)
      )
  ) %>%
    dplyr::mutate(Value = Value_raw) %>%
    summarise_by_repeat()

  rmse_long <- dplyr::bind_rows(
    df %>%
      dplyr::transmute(
        model,
        n_tuning_seeds,
        Type = "Pooled",
        Value_raw = as.numeric(mean_pooled_rmse)
      ),
    df %>%
      dplyr::transmute(
        model,
        n_tuning_seeds,
        Type = "Mean fold",
        Value_raw = as.numeric(mean_mean_rmse)
      )
  ) %>%
    dplyr::mutate(Value = Value_raw) %>%
    summarise_by_repeat()

  r2_long <- r2_long %>%
    dplyr::mutate(
      model = factor(model, levels = model_list),
      Type = factor(Type, levels = c("Pooled", "Mean fold"))
    )
  rmse_long <- rmse_long %>%
    dplyr::mutate(
      model = factor(model, levels = model_list),
      Type = factor(Type, levels = c("Pooled", "Mean fold"))
    )

  linestyle_values <- c(
    "Pooled" = unname(METRIC_LINESTYLES["Pooled R2"]),
    "Mean fold" = unname(METRIC_LINESTYLES["Mean-fold R2"])
  )

  p_r2 <- ggplot2::ggplot(
    r2_long,
    ggplot2::aes(
      x = n_tuning_seeds,
      y = Value,
      color = model,
      linetype = Type,
      group = interaction(model, Type)
    )
  ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = Value - SE, ymax = Value + SE),
      width = 0.22,
      alpha = 0.75,
      na.rm = TRUE
    ) +
    ggplot2::scale_color_manual(values = MODEL_COLOURS, breaks = model_list, drop = FALSE) +
    ggplot2::scale_linetype_manual(values = linestyle_values, drop = FALSE) +
    ggplot2::labs(
      x = "Number of model tuning seeds",
      y = expression(R^2*" value"),
      color = "Model",
      linetype = "Type"
    ) +
    theme_paper()

  p_rmse <- ggplot2::ggplot(
    rmse_long,
    ggplot2::aes(
      x = n_tuning_seeds,
      y = Value,
      color = model,
      linetype = Type,
      group = interaction(model, Type)
    )
  ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = Value - SE, ymax = Value + SE),
      width = 0.22,
      alpha = 0.75,
      na.rm = TRUE
    ) +
    ggplot2::scale_color_manual(values = MODEL_COLOURS, breaks = model_list, drop = FALSE) +
    ggplot2::scale_linetype_manual(values = linestyle_values, drop = FALSE) +
    ggplot2::labs(
      x = "Number of model tuning seeds",
      y = "RMSE",
      color = "Model",
      linetype = "Type"
    ) +
    theme_paper()

  combined <- p_r2 / p_rmse + patchwork::plot_layout(guides = "collect")

  out_r2 <- file.path(out_dir, TUNING_SWEEP_OUTPUT_PNGS[[1]])
  out_rmse <- file.path(out_dir, TUNING_SWEEP_OUTPUT_PNGS[[2]])
  out_combined <- file.path(out_dir, TUNING_SWEEP_OUTPUT_PNGS[[3]])

  ggplot2::ggsave(out_r2, p_r2, width = 9.5, height = 4.6, dpi = dpi)
  ggplot2::ggsave(out_rmse, p_rmse, width = 9.5, height = 4.6, dpi = dpi)
  ggplot2::ggsave(out_combined, combined, width = 10.5, height = 8.4, dpi = dpi)

  plot_log_tss(paste0("Wrote: ", rel_path_tss(out_r2)))
  plot_log_tss(paste0("Wrote: ", rel_path_tss(out_rmse)))
  plot_log_tss(paste0("Wrote: ", rel_path_tss(out_combined)))
  invisible(combined)
}

plot_tuning_seed_sweep <- function(sweep_out_dir = NULL,
                                            model_list = cfg$model_list,
                                            dpi = cfg$dpi) {
  sweep_out_dir <- resolve_sweep_out_dir(project_root, cfg, sweep_out_dir = sweep_out_dir)
  if (!dir.exists(sweep_out_dir)) {
    stop("sweep_out_dir does not exist: ", sweep_out_dir)
  }
  sweep_csv <- file.path(sweep_out_dir, TUNING_SWEEP_REQUIRED_CSV)
  if (!file.exists(sweep_csv)) {
    stop("Missing sweep summary CSV: ", sweep_csv, call. = FALSE)
  }

  plot_log_tss("Sources: modelling/R/init_repo.R, modelling/pipeline_config.R, modelling/plots/plot_config.R, modelling/plots/plot_tuning_seed_sweep.R")
  plot_log_tss(paste0("Input CSV: ", rel_path_tss(sweep_csv)))
  plot_log_tss(paste0("Requested model_list: ", paste(as.character(model_list), collapse = ", ")))
  plot_log_tss(paste0("Output directory: ", rel_path_tss(sweep_out_dir)))
  sweep_df <- readr::read_csv(sweep_csv, show_col_types = FALSE)
  if (nrow(sweep_df) == 0L) {
    warning("Sweep summary CSV is empty: ", sweep_csv, call. = FALSE)
    return(invisible(NULL))
  }
  plot_tuning_seed_sweep_summary(
    sweep_df = sweep_df,
    out_dir = sweep_out_dir,
    project_root = project_root,
    dpi = dpi,
    model_list = model_list
  )
  invisible(sweep_df)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  sweep_out_dir <- if (length(args) >= 1L && nzchar(args[[1L]])) args[[1L]] else NULL
  dpi <- if (length(args) >= 2L && nzchar(args[[2L]])) suppressWarnings(as.integer(args[[2L]])) else cfg$dpi
  plot_tuning_seed_sweep(sweep_out_dir = sweep_out_dir, model_list = cfg$model_list, dpi = dpi)
} else {
  disable_autorun <- isTRUE(get0("plot_tuning_seed_sweep_disable_autorun_on_source", envir = .GlobalEnv, ifnotfound = FALSE))
  if (interactive() && !disable_autorun) {
    plot_log_tss("Sourced in interactive session; auto-running plot generation from current config.")
    plot_tuning_seed_sweep(sweep_out_dir = NULL, model_list = cfg$model_list, dpi = cfg$dpi)
  } else {
    message(
      "[plot:plot_tuning_seed_sweep] Loaded.\n",
      "Sources: modelling/R/init_repo.R, modelling/pipeline_config.R, modelling/plots/plot_config.R, modelling/plots/plot_tuning_seed_sweep.R\n",
      "Required CSV:\n  ", file.path("<sweep_out_dir>", TUNING_SWEEP_REQUIRED_CSV), "\n",
      "Written PNGs:\n  ", paste(file.path("<sweep_out_dir>", TUNING_SWEEP_OUTPUT_PNGS), collapse = "\n  "), "\n",
      "Run interactively with:\n",
      "  plot_tuning_seed_sweep()\n",
      "  # optional override: plot_tuning_seed_sweep('output/tuning_seed_sweep_runs/sweep_001')\n",
      "Or from shell with:\n",
      "  Rscript modelling/plots/plot_tuning_seed_sweep.R\n",
      "  # optional override: Rscript modelling/plots/plot_tuning_seed_sweep.R <sweep_out_dir> [dpi]"
    )
  }
}
