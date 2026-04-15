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
  packages = c("here", "readr", "dplyr", "tidyr", "ggplot2"),
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

  source(file.path(project_root, "modelling/plots/plot_tuning_seed_sweep.R"))
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
