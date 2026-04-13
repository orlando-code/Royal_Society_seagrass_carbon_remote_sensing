## Plot robust multi-seed SHAP importance with error bars.
##
## Inputs:
##   output/<cv_regime_name>/covariate_selection/robust_pixel_grouped/
##     shap_importance_summary_robust_pixel_grouped_seeds_<...>.csv
##
## Outputs:
##   output/<cv_regime_name>/covariate_selection/robust_pixel_grouped/
##     shap_importance_robust_errorbars_by_model.png
##     shap_importance_robust_errorbars_combined.png

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
source(file.path(project_root, "modelling/plots/plot_config.R"))
source(file.path(project_root, "modelling/pipeline_config.R"))

cfg <- get_pipeline_config()
apply_pipeline_defaults(
  cfg,
  c(
    "dpi", "show_titles", "robust_fold_seed_list", "cv_regime_name", "robust_shap_plot_top_n"
  ),
  envir = .GlobalEnv
)

load_packages(c("dplyr", "ggplot2", "patchwork", "readr", "here", "ggtext"))
# show_titles <- get0("show_titles", envir = .GlobalEnv, ifnotfound = FALSE)
# cv_regime_name <- get0("cv_regime_name", envir = .GlobalEnv, ifnotfound = "pixel_grouped")
# robust_fold_seed_list <- as.integer(get0("robust_fold_seed_list", envir = .GlobalEnv, ifnotfound = integer()))
if (length(robust_fold_seed_list) == 0L) {
  stop("robust_fold_seed_list is missing in .GlobalEnv; run from robust pipeline driver.")
}

seeds_str <- paste(robust_fold_seed_list, collapse = "-")
run_output_dir <- get0("run_output_dir", envir = .GlobalEnv, ifnotfound = NA_character_)
if (is.na(run_output_dir) || !nzchar(as.character(run_output_dir))) {
  stop("run_output_dir must be set in .GlobalEnv for robust SHAP importance plots.")
}
run_output_dir <- as.character(run_output_dir)
subset_work_root <- if (basename(run_output_dir) == "evaluation") dirname(run_output_dir) else run_output_dir
robust_cov_dir <- if (basename(run_output_dir) == "evaluation") {
  file.path(subset_work_root, "covariates")
} else {
  file.path(run_output_dir, "covariate_selection", "robust_pixel_grouped")
}
in_summary <- file.path(
  robust_cov_dir,
  paste0("shap_importance_summary_robust_pixel_grouped_seeds_", seeds_str, ".csv")
)
if (!file.exists(in_summary)) {
  stop("Missing SHAP summary CSV: ", in_summary)
}

imp <- readr::read_csv(in_summary, show_col_types = FALSE)
req <- c("model", "variable", "shap_importance_mean", "shap_importance_sd")
if (!all(req %in% names(imp))) {
  stop("SHAP summary CSV missing required columns: ", paste(setdiff(req, names(imp)), collapse = ", "))
}

top_n <- as.integer(get0("robust_shap_plot_top_n", envir = .GlobalEnv, ifnotfound = 15L))
imp_top <- imp %>%
  dplyr::group_by(model) %>%
  dplyr::slice_max(order_by = shap_importance_mean, n = top_n, with_ties = FALSE) %>%
  dplyr::arrange(model, dplyr::desc(shap_importance_mean)) %>%
  dplyr::ungroup()

shared_counts <- imp_top %>%
  dplyr::distinct(model, variable) %>%
  dplyr::count(variable, name = "n_models")

imp_top <- imp_top %>%
  dplyr::left_join(shared_counts, by = "variable") %>%
  dplyr::mutate(
    variable_label = label_vars(variable),
    variable_styled = dplyr::case_when(
      n_models >= 3L ~ paste0("**", variable_label, "**"),
      n_models == 1L ~ paste0("*", variable_label, "*"),
      TRUE ~ variable_label
    ),
    ymin = pmax(0, shap_importance_mean - shap_importance_sd),
    ymax = shap_importance_mean + shap_importance_sd
  )

make_model_plot <- function(dfm, m, legend_models) {
  dfm <- dplyr::arrange(dfm, shap_importance_mean)
  dfm$model <- factor(dfm$model, levels = legend_models)
  ggplot2::ggplot(
    dfm,
    ggplot2::aes(
      x = shap_importance_mean,
      y = reorder(variable_styled, shap_importance_mean),
      fill = model
    )
  ) +
    ggplot2::geom_col(alpha = 0.85) +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = ymin, xmax = ymax), width = 0.2, colour = "black", linewidth = 0.35) +
    ggplot2::scale_fill_manual(
      values = MODEL_COLOURS[legend_models],
      breaks = legend_models,
      limits = legend_models,
      name = "Model",
      drop = FALSE
    ) +
    ggplot2::labs(
      title = if (show_titles) m else NULL,
      x = "Mean SHAP",
      y = NULL
    ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      axis.text.y = ggtext::element_markdown(),
      axis.text.y.left = ggtext::element_markdown()
    )
}

models <- unique(imp_top$model)

plots <- lapply(models, function(m) {
  make_model_plot(
    imp_top %>% dplyr::filter(model == m),
    m,
    models
  )
})

names(plots) <- models

combined <- patchwork::wrap_plots(
  plots,
  ncol = length(plots)
) +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(
    title = if (show_titles)
      "Robust multi-seed SHAP importance (mean ± SD)" else NULL,
    subtitle = if (show_titles)
      paste0("Top ", top_n, " variables per model; seeds=", seeds_str) else NULL
  ) &
  ggplot2::theme(
    legend.position = "bottom",
    legend.direction = "horizontal"
  )

out_combined <- file.path(robust_cov_dir, "shap_importance_robust_errorbars_combined.png")
ggsave(out_combined, combined, width = max(10, 4 * length(plots)), height = 7, dpi = dpi)

cat("Saved robust SHAP error-bar plots:\n")
cat("  ", out_combined, "\n", sep = "")
