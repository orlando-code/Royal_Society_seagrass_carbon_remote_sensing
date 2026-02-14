# GPR-only predictor selection: correlation filter + backward elimination
#
# 1. Apply correlation filter (|r| > 0.8 between predictors; keep one per pair).
# 2. Baseline = mean 1 km spatial CV R² for GPR with that non-correlated set.
# 3. Backward elimination: remove one variable at a time only if mean R² stays >= baseline.
#
# Output: figures/covariate_selection/pruned_variables_to_include_gpr.csv (variable list)
#         figures/cv_pipeline_output/gpr_r2_pruning_summary.csv (baseline R², final R², n_removed)
#
# Requires: spatial_folds_cache.rds from cv_pipeline (run paper figures / cv_pipeline first).

# rm(list = ls())
setwd(here::here())
set.seed(42)

source("modelling/R/helpers.R")
source("modelling/R/extract_covariates_from_rasters.R")
load_packages(c("here", "tidyverse", "blockCV", "sf"))

target_var <- "median_carbon_density_100cm"
out_dir <- "figures/cv_pipeline_output"
covariate_dir <- "figures/covariate_selection"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(covariate_dir, recursive = TRUE, showWarnings = FALSE)

cor_threshold <- 0.8
block_size_m <- 1000
block_size_name <- paste0("spatial_block_", block_size_m, "m")

# ---------------------------------------------------------------------------
# 1. Load core-level data (same as model_permutation_pruning.R)
# ---------------------------------------------------------------------------
dat <- read_rds("data/all_extracted_new.rds")
predictor_vars_full <- raster_covariates[raster_covariates %in% colnames(dat)]
if (length(predictor_vars_full) == 0) {
  stop("No predictor_vars found in dat. Check all_extracted_new.rds and raster config.")
}

if ("random_core_variable" %in% colnames(dat)) {
  core_data <- dat %>%
    dplyr::group_by(random_core_variable) %>%
    dplyr::summarise(
      median_carbon_density = median(.data[[target_var]], na.rm = TRUE),
      dplyr::across(dplyr::all_of(c(predictor_vars_full, "longitude", "latitude")), first),
      .groups = "drop"
    ) %>%
    dplyr::filter(!is.na(median_carbon_density))
} else {
  core_data <- dat %>%
    dplyr::mutate(median_carbon_density = .data[[target_var]]) %>%
    dplyr::select(longitude, latitude, median_carbon_density, dplyr::all_of(predictor_vars_full))
}

core_data_complete <- core_data %>%
  dplyr::select(longitude, latitude, median_carbon_density, dplyr::all_of(predictor_vars_full)) %>%
  dplyr::filter(complete.cases(.))
predictor_vars_full <- predictor_vars_full[predictor_vars_full %in% colnames(core_data_complete)]
n_cores <- nrow(core_data_complete)

cat("Cores (complete cases):", n_cores, "\n")
cat("Predictors (all):", length(predictor_vars_full), "\n\n")

# ---------------------------------------------------------------------------
# 2. Correlation filter (keep one per |r| > cor_threshold pair)
# ---------------------------------------------------------------------------
predictor_vars <- predictor_vars_full
if (length(predictor_vars) >= 2) {
  cat("Applying correlation filter (|r| >", cor_threshold, ")...\n")
  dat_complete <- core_data_complete[, c("median_carbon_density", predictor_vars), drop = FALSE]
  dat_complete <- dat_complete[complete.cases(dat_complete), ]
  if (nrow(dat_complete) >= 3) {
    cor_with_target <- sapply(predictor_vars, function(v) {
      if (v %in% names(dat_complete)) {
        stats::cor(dat_complete[[v]], dat_complete[["median_carbon_density"]], use = "complete.obs")
      } else {
        NA_real_
      }
    })
    cor_df <- data.frame(
      variable = predictor_vars,
      correlation = cor_with_target,
      abs_correlation = abs(cor_with_target),
      stringsAsFactors = FALSE
    )
    cor_df <- cor_df[!is.na(cor_df$correlation), ]
    vars_for_cor <- cor_df$variable

    if (length(vars_for_cor) >= 2) {
      covar_data <- dat_complete[, vars_for_cor, drop = FALSE]
      cor_matrix <- stats::cor(covar_data, use = "complete.obs")
      high_cor_pairs <- which(abs(cor_matrix) > cor_threshold & abs(cor_matrix) < 1, arr.ind = TRUE)
      vars_removed <- character()

      if (nrow(high_cor_pairs) > 0) {
        pair_df <- data.frame(
          var1 = rownames(cor_matrix)[high_cor_pairs[, 1]],
          var2 = colnames(cor_matrix)[high_cor_pairs[, 2]],
          correlation = cor_matrix[high_cor_pairs],
          stringsAsFactors = FALSE
        )
        pair_df <- pair_df[pair_df$var1 < pair_df$var2, ]
        pair_df <- pair_df[order(abs(pair_df$correlation), decreasing = TRUE), ]

        for (i in seq_len(nrow(pair_df))) {
          v1 <- pair_df$var1[i]
          v2 <- pair_df$var2[i]
          if (v1 %in% vars_removed || v2 %in% vars_removed) next
          c1 <- cor_df$abs_correlation[cor_df$variable == v1]
          c2 <- cor_df$abs_correlation[cor_df$variable == v2]
          if (length(c1) == 0 || length(c2) == 0) next
          if (round(c1, 3) == round(c2, 3)) next
          if (c1 < c2) {
            vars_removed <- c(vars_removed, v1)
          } else {
            vars_removed <- c(vars_removed, v2)
          }
        }
      }

      vars_removed <- unique(vars_removed)
      if (length(vars_removed) > 0) {
        cat("  Removed", length(vars_removed), "predictors:", paste(vars_removed, collapse = ", "), "\n")
      }
      predictor_vars <- setdiff(predictor_vars, vars_removed)
    }
  }
  cat("Predictors after correlation filter:", length(predictor_vars), "\n\n")
}

# # ---------------------------------------------------------------------------
# # 3. Load 1 km spatial folds
# # ---------------------------------------------------------------------------
# spatial_folds_cache_path <- file.path(out_dir, "spatial_folds_cache.rds")
# if (!file.exists(spatial_folds_cache_path)) {
#   stop("spatial_folds_cache.rds not found. Run cv_pipeline / run_paper_figures first.")
# }

# cached <- tryCatch(readRDS(spatial_folds_cache_path), error = function(e) NULL)
# if (is.null(cached) || !"spatial_strategies" %in% names(cached)) {
#   stop("spatial_folds_cache.rds does not contain spatial_strategies.")
# }

# strategy_block_size <- NULL
# for (s in cached$spatial_strategies) {
#   if (s$method == block_size_name) {
#     strategy_block_size <- s
#     break
#   }
# }
# if (is.null(strategy_block_size)) {
#   stop(paste0("Strategy ", block_size_name, " not found in spatial_folds_cache.rds."))
# }

# folds_ids <- strategy_block_size$folds
# if (length(folds_ids) != nrow(core_data_complete)) {
#   stop("Folds length != nrow(core_data_complete).")
# }
# cat("Using", length(unique(folds_ids)), "folds for", block_size_name, "\n\n")

# core_sf <- sf::st_as_sf(
#   core_data_complete,
#   coords = c("longitude", "latitude"),
#   crs = 4326,
#   remove = FALSE
# )

# # ---------------------------------------------------------------------------
# # 4. Helper: run 1 km spatial CV (GPR only), return mean R²
# # ---------------------------------------------------------------------------
# run_gpr_cv_r2 <- function(core_data_use, predictor_vars_use) {
#   if (length(predictor_vars_use) < 1) return(NA_real_)
#   res <- run_cv(
#     cv_method_name = block_size_name,
#     fold_indices = folds_ids,
#     core_data = core_data_use,
#     predictor_vars = predictor_vars_use,
#     tune_hyperparams = FALSE,
#     nested_tuning = TRUE,
#     verbose = FALSE,
#     buffer_m = NA_real_,
#     core_sf = core_sf,
#     return_predictions = FALSE,
#     skip_inla = TRUE,
#     models = "GPR"
#   )
#   gpr <- res[res$model == "GPR", ]
#   if (nrow(gpr) == 0) return(NA_real_)
#   mean(gpr$r2, na.rm = TRUE)
# }

# # ---------------------------------------------------------------------------
# # 5. Baseline R² (non-correlated set only)
# # ---------------------------------------------------------------------------
# cat("Computing baseline GPR R² (correlation-filtered set)...\n")
# baseline_r2 <- run_gpr_cv_r2(core_data_complete, predictor_vars)
# cat("Baseline mean R² (GPR):", round(baseline_r2, 4), "\n\n")

# # ---------------------------------------------------------------------------
# # 6. Backward elimination: remove one at a time only if R² >= baseline
# # ---------------------------------------------------------------------------
# current_vars <- predictor_vars
# n_removed <- 0
# elimination_log <- list()

# while (length(current_vars) > 1) {
#   best_removal <- NULL
#   best_r2_after <- -Inf

#   for (v in current_vars) {
#     try_vars <- setdiff(current_vars, v)
#     r2_after <- run_gpr_cv_r2(core_data_complete, try_vars)
#     if (is.finite(r2_after) && r2_after >= baseline_r2 && r2_after > best_r2_after) {
#       best_r2_after <- r2_after
#       best_removal <- v
#     }
#   }

#   if (is.null(best_removal)) break

#   current_vars <- setdiff(current_vars, best_removal)
#   n_removed <- n_removed + 1
#   elimination_log[[length(elimination_log) + 1]] <- data.frame(
#     removed = best_removal,
#     n_remaining = length(current_vars),
#     mean_r2 = best_r2_after,
#     stringsAsFactors = FALSE
#   )
#   cat("  Removed:", best_removal, "| remaining:", length(current_vars), "| mean R²:", round(best_r2_after, 4), "\n")
# }

# final_r2 <- if (length(current_vars) > 0) {
#   run_gpr_cv_r2(core_data_complete, current_vars)
# } else {
#   NA_real_
# }

# cat("\nBackward elimination done. Removed", n_removed, "variables.\n")
# cat("Final predictor count:", length(current_vars), "\n")
# cat("Final mean R² (GPR):", round(final_r2, 4), "\n")

# ---------------------------------------------------------------------------
# 7. Write outputs
# ---------------------------------------------------------------------------
current_vars <- predictor_vars

pruned_csv <- file.path(covariate_dir, "pruned_variables_to_include_gpr.csv")
write.csv(data.frame(variable = current_vars), pruned_csv, row.names = FALSE)
cat("Wrote:", pruned_csv, "\n")
if (exists("combine_pruned_model_variables", mode = "function")) {
  out_combined <- combine_pruned_model_variables(covariate_dir)
  if (!is.null(out_combined)) cat("Updated combined file:", out_combined, "\n")
}

summary_df <- data.frame(
  baseline_r2 = baseline_r2,
  final_r2 = final_r2,
  n_after_cor_filter = length(predictor_vars),
  n_removed_backward = n_removed,
  n_final = length(current_vars),
  stringsAsFactors = FALSE
)
summary_path <- file.path(out_dir, "gpr_r2_pruning_summary.csv")
write.csv(summary_df, summary_path, row.names = FALSE)
cat("Wrote:", summary_path, "\n")

if (length(elimination_log) > 0) {
  elim_path <- file.path(out_dir, "gpr_r2_pruning_elimination_log.csv")
  write.csv(dplyr::bind_rows(elimination_log), elim_path, row.names = FALSE)
  cat("Wrote:", elim_path, "\n")
}
