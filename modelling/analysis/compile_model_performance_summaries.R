## Build communication-ready summaries of model performance.
##
## Writes:
##   output/<cv_regime>/cv_pipeline/model_performance_across_methods_summary.csv
##   output/<cv_regime>/cv_pipeline/deployment_robustness_pixel_grouped.csv

project_root <- here::here()
setwd(project_root)

source(file.path(project_root, "modelling/R/helpers.R"))
load_packages(c("here", "dplyr", "readr"))

cv_regime_name <- get0("cv_regime_name", envir = .GlobalEnv, ifnotfound = "pixel_grouped")

base_dir <- file.path("output", cv_regime_name, "cv_pipeline")
if (!dir.exists(base_dir)) stop("Missing directory: ", base_dir)

sum_path <- file.path(base_dir, "post_tuning_cv_results_summary.csv")
det_path <- file.path(base_dir, "post_tuning_cv_results_detailed.csv")
if (!file.exists(sum_path)) stop("Missing: ", sum_path)
if (!file.exists(det_path)) stop("Missing: ", det_path)

post_summary <- readr::read_csv(sum_path, show_col_types = FALSE)
post_detailed <- readr::read_csv(det_path, show_col_types = FALSE)

neg_counts <- post_detailed %>%
  group_by(method, model) %>%
  summarise(
    neg_r2_fraction = mean(r2 < 0, na.rm = TRUE),
    min_r2 = min(r2, na.rm = TRUE),
    .groups = "drop"
  )

methods_tbl <- post_summary %>%
  left_join(neg_counts, by = c("method", "model")) %>%
  arrange(method, mean_rmse)

write.csv(
  methods_tbl,
  file.path(base_dir, "model_performance_across_methods_summary.csv"),
  row.names = FALSE
)

# ---------------------------------------------------------------------------
# Deployment robustness check (multi-seed robust evaluation outputs)
# ---------------------------------------------------------------------------
robust_eval_candidates <- c(
  Sys.glob(file.path(base_dir, "pixel_grouped_evaluation_*x_*_seeds*")),
  Sys.glob(file.path(base_dir, "robust_pixel_grouped_evaluation_robustSeeds_*"))
)
robust_eval_candidates <- unique(robust_eval_candidates)
robust_eval_candidates <- robust_eval_candidates[dir.exists(robust_eval_candidates)]
legacy_repeated_dir <- file.path(base_dir, "repeated_pixel_grouped_seed_check")

if (length(robust_eval_candidates) > 0L) {
  # Prefer the newest robust run directory if multiple exist.
  dir_info <- file.info(robust_eval_candidates)
  robust_dir <- robust_eval_candidates[which.max(dir_info$mtime)]
} else {
  robust_dir <- legacy_repeated_dir
}

rob_sum_path <- file.path(robust_dir, "across_seeds_summary.csv")
rob_detail_path <- file.path(robust_dir, "by_seed_detailed.csv")

if (file.exists(rob_sum_path) && file.exists(rob_detail_path)) {
  rob_sum <- readr::read_csv(rob_sum_path, show_col_types = FALSE)
  rob_detail <- readr::read_csv(rob_detail_path, show_col_types = FALSE)

  rob_neg <- rob_detail %>%
    group_by(model) %>%
    summarise(neg_r2_fraction = mean(r2 < 0, na.rm = TRUE), min_r2 = min(r2, na.rm = TRUE), .groups = "drop")

  rob_tbl <- rob_sum %>%
    left_join(rob_neg, by = "model") %>%
    arrange(mean_mean_rmse)

  write.csv(
    rob_tbl,
    file.path(base_dir, "deployment_robustness_pixel_grouped.csv"),
    row.names = FALSE
  )
}

cat("Wrote summaries to:\n",
    file.path(base_dir, "model_performance_across_methods_summary.csv"), "\n",
    file.path(base_dir, "deployment_robustness_pixel_grouped.csv"), "\n",
    sep = "")

