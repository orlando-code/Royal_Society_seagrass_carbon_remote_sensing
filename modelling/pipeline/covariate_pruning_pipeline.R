# Covariate Pruning Pipeline
#
# Global variables (set by run_paper_figures.R):
#   use_correlation_filter, correlation_filter_threshold
#   use_fair_comparison        – TRUE = one covariate set for all models
#   fair_comparison_max_vars   – max vars retained after permutation importance
#   fair_comparison_reference_model – "RF" | "GAM" | "XGB" | "GPR"
#   model_list                 – character vector of model names
#
# Outputs -> figures/covariate_selection/
#   pruned_variables_to_include.csv
#   pruned_variables_to_include_<model>.csv  (one per model in model_list)

source("modelling/R/helpers.R")
source("modelling/R/gpr_funs.R")
source("modelling/R/assign_region_from_latlon.R")
load_packages(c("dplyr", "readr", "blockCV", "sf"))

# ---------------------------------------------------------------------------
# Data and region exclusion
# ---------------------------------------------------------------------------
dat <- readr::read_rds("data/all_extracted_new.rds")
exclude_regions <- get0("exclude_regions", envir = .GlobalEnv, ifnotfound = character(0))
if (length(exclude_regions) > 0L) {
  if (!"region" %in% names(dat)) dat <- assign_region_from_latlon(dat)
  dat <- dat[is.na(dat$region) | !dat$region %in% exclude_regions, ]
  cat("Excluded regions:", paste(exclude_regions, collapse = ", "), "\n")
}

if (!exists("model_list")) model_list <- c("GPR", "GAM", "XGB")
target_var <- "median_carbon_density_100cm"
predictor_vars <- setdiff(colnames(dat),
  c("latitude", "longitude", "number_id_final_version",
    "seagrass_species", "region", target_var))

# ---------------------------------------------------------------------------
# Step i: Correlation filter
# ---------------------------------------------------------------------------
if (isTRUE(get0("use_correlation_filter"))) {
  cat("Pruning by correlation (threshold:", correlation_filter_threshold, ")...\n")
  predictor_vars <- prune_by_correlation(dat, predictor_vars, target_var,
    cor_threshold = correlation_filter_threshold)
  cat("Variables after correlation filter:", length(predictor_vars), "\n")
}

# ---------------------------------------------------------------------------
# Step ii: Fair comparison – shared covariate set via reference-model permutation
# ---------------------------------------------------------------------------
use_fair_comparison <- get0("use_fair_comparison", ifnotfound = TRUE)
fair_comparison_max_vars <- get0("fair_comparison_max_vars", ifnotfound = 15L)
fair_comparison_reference_model <- get0("fair_comparison_reference_model", ifnotfound = "RF")

if (use_fair_comparison && length(predictor_vars) > fair_comparison_max_vars) {
  cat("Fair comparison: selecting top ", fair_comparison_max_vars, " vars via ",
      fair_comparison_reference_model, " permutation importance (1 km spatial CV).\n", sep = "")

  fair_folds_cache <- "figures/cv_pipeline_output/fair_comparison_1km_folds.rds"
  fold_cache_key <- list(n_cores = nrow(dat), exclude_regions = exclude_regions, block_size = 1000L, k = 5L)
  fold_indices <- NULL
  if (file.exists(fair_folds_cache)) {
    cached <- tryCatch(readRDS(fair_folds_cache), error = function(e) NULL)
    if (!is.null(cached) && identical(cached$cache_key, fold_cache_key)) {
      fold_indices <- cached$fold_indices
      cat("  Using cached 1 km spatial folds.\n")
    }
  }
  if (is.null(fold_indices)) {
    cat("  Creating 1 km spatial folds...\n")
    core_sf <- sf::st_as_sf(dat, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)
    fold_indices <- blockCV::cv_spatial(
      x = core_sf, k = 5L, size = 1000L, selection = "random",
      iteration = 50, progress = FALSE, biomod2 = FALSE, hexagon = FALSE, plot = FALSE
    )$folds_ids
    dir.create(dirname(fair_folds_cache), recursive = TRUE, showWarnings = FALSE)
    saveRDS(list(cache_key = fold_cache_key, fold_indices = fold_indices), fair_folds_cache)
    cat("  Cached 1 km folds to:", fair_folds_cache, "\n")
  }

  cat("  Permutation importance (", length(predictor_vars), " variables):\n", sep = "")
  imp_df <- permutation_importance_spatial_cv(
    dat, predictor_vars, fair_comparison_reference_model,
    fold_indices, n_permutations = 1L, verbose = TRUE)
  cum_frac <- cumsum(pmax(imp_df$rmse_increase, 0))
  cum_frac <- cum_frac / max(cum_frac, na.rm = TRUE)
  # top N by importance, capped at 99% cumulative coverage (min 2)
  n_keep <- max(2L, min(fair_comparison_max_vars, sum(cum_frac <= 0.99, na.rm = TRUE) + 1L, nrow(imp_df)))
  predictor_vars <- imp_df$variable[seq_len(n_keep)]
  cat("  Selected", length(predictor_vars), "variables for fair comparison.\n")
}

# ---------------------------------------------------------------------------
# Write pruned variable sets
# ---------------------------------------------------------------------------
covariate_dir <- "figures/covariate_selection"
dir.create(covariate_dir, recursive = TRUE, showWarnings = FALSE)
pv_df <- data.frame(variable = predictor_vars)
write.csv(pv_df, file.path(covariate_dir, "pruned_variables_to_include.csv"), row.names = FALSE)
cat("  Saved: pruned_variables_to_include.csv (", length(predictor_vars), " vars)\n", sep = "")
for (model in model_list) {
  out_path <- file.path(covariate_dir, paste0("pruned_variables_to_include_", tolower(model), ".csv"))
  write.csv(pv_df, out_path, row.names = FALSE)
  cat("  Saved:", out_path, "\n")
}

out_combined <- combine_pruned_model_variables(covariate_dir)
if (!is.null(out_combined)) {
  n_files <- length(list.files(covariate_dir, pattern = "^pruned_variables_to_include_.+\\.csv$"))
  cat("Combined", n_files, "pruned variable files ->", out_combined, "\n")
} else {
  cat("No pruned_variables_to_include_*.csv found; skipping pruned_model_variables.csv\n")
}
cat("Covariate pruning complete.\n")