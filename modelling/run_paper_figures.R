# =============================================================================
# Paper pipeline driver: exploration, variable choice, prediction, supplement
#
# Runs the full pipeline in order to regenerate all paper and supplement figures.
# Uses cached data where possible (spatial folds, tuning results) to avoid
# re-running heavy steps. All outputs go to figures/cv_pipeline_output,
# figures/interpolation, and figures/supplement.
#
# Usage: setwd(project_root); source("modelling/run_paper_figures.R")
#        Or: Rscript modelling/run_paper_figures.R
# =============================================================================

if (!requireNamespace("here", quietly = TRUE)) install.packages("here", quiet = TRUE)
setwd(here::here())
set.seed(42)

cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("PAPER PIPELINE: EXPLORATION -> CHOICE -> PREDICTION -> SUPPLEMENT\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# -----------------------------------------------------------------------------
# -1. 
# -----------------------------------------------------------------------------
cat("Step -1: Confirming paper pipeline configuration ...\n")
# pruning steps
correlation_filter_on <- TRUE
correlation_filter_threshold <- 0.8
model_permutation_pruning_on <- TRUE
model_permutation_pruning_keep_frac <- 0.99
cat("  correlation_filter_on:", correlation_filter_on, "\n")
cat("  correlation_filter_threshold:", correlation_filter_threshold, "\n")
cat("  model_permutation_pruning_on:", model_permutation_pruning_on, "\n")
cat("  model_permutation_pruning_keep_frac:", model_permutation_pruning_keep_frac, "\n")
cat("\n")


# -----------------------------------------------------------------------------
# 0. Data (build if missing)
# -----------------------------------------------------------------------------
if (!file.exists("data/all_extracted_new.rds")) {
  cat("Step 0: Building all_extracted_new.rds ...\n")
  source("modelling/pipeline/build_all_extracted_new.R")
  cat("\n")
} else {
  cat("Step 0: Using existing data/all_extracted_new.rds\n\n")
}

# -----------------------------------------------------------------------------
# 1. CV pipeline (folds, CV results, inline plots) -> figures/cv_pipeline_output
# -----------------------------------------------------------------------------
cat("Step 1: CV pipeline (spatial folds, comparison, results)\n")
source("modelling/pipeline/cv_pipeline.R")
cat("\n")

# -----------------------------------------------------------------------------
# 2. Model permutation pruning (variable selection) -> figures/cv_pipeline_output
#    Uses correlation_filter_on, correlation_filter_threshold, model_permutation_pruning_on, model_permutation_pruning_keep_frac
# -----------------------------------------------------------------------------
config_path <- "figures/cv_pipeline_output/paper_pipeline_config.rds"
dir.create("figures/cv_pipeline_output", recursive = TRUE, showWarnings = FALSE)
saveRDS(list(
  correlation_filter_on = correlation_filter_on,
  correlation_filter_threshold = correlation_filter_threshold,
  model_permutation_pruning_on = model_permutation_pruning_on,
  model_permutation_pruning_keep_frac = model_permutation_pruning_keep_frac
), config_path)

# Run Step 2 when we need correlation filtering and/or permutation pruning
if (correlation_filter_on || model_permutation_pruning_on) {
  cat("Step 2: Predictor selection (correlation filter:", correlation_filter_on, "| permutation pruning:", model_permutation_pruning_on, ")\n")
  source("modelling/pipeline/model_permutation_pruning.R")
  pruning_csv <- list.files("figures/cv_pipeline_output", pattern = "^cv_model_permutation_pruning_.*\\.csv$", full.names = TRUE)
  if (length(pruning_csv) > 0L) {
    pr <- read.csv(pruning_csv[1L], stringsAsFactors = FALSE)
    keep_ok <- (pr$keep == TRUE) | (tolower(as.character(pr$keep)) %in% c("true", "t", "1"))
    if (all(c("model", "variable", "keep") %in% names(pr))) {
      dir.create("figures/covariate_selection", recursive = TRUE, showWarnings = FALSE)
      gpr_vars <- unique(pr$variable[pr$model == "GPR" & keep_ok])
      if (length(gpr_vars) > 0L) {
        write.csv(data.frame(variable = gpr_vars), "figures/covariate_selection/pruned_variables_to_include_gpr.csv", row.names = FALSE)
      }
      gam_vars <- unique(pr$variable[pr$model == "GAM" & keep_ok])
      if (length(gam_vars) > 0L) {
        write.csv(data.frame(variable = gam_vars), "figures/covariate_selection/pruned_variables_to_include.csv", row.names = FALSE)
      }
    }
  }
} else {
  cat("Step 2: Skipped (correlation_filter_on and model_permutation_pruning_on both FALSE)\n")
}
cat("\n")

# -----------------------------------------------------------------------------
# 3. GPR hyperparameter tuning -> figures/cv_pipeline_output
# -----------------------------------------------------------------------------
cat("Step 3: GPR hyperparameter tuning\n")
source("modelling/pipeline/gpr_hyperparameter_tuning.R")
cat("\n")

# -----------------------------------------------------------------------------
# 4. GPR categorical and spatial investigations (feed into parameter summary plot)
# -----------------------------------------------------------------------------
cat("Step 4a: GPR categorical investigation (species/region)\n")
source("modelling/gpr/gpr_categorical_investigation.R")
cat("\n")
cat("Step 4b: GPR spatial terms investigation (lat/lon, region)\n")
source("modelling/gpr/gpr_spatial_investigation.R")
cat("\n")
cat("Step 4c: GPR parameter search summary (tuning plot)\n")
source("modelling/gpr/gpr_parameter_search_summary.R")
cat("\n")

# -----------------------------------------------------------------------------
# 5. Plot-only: CV figures and permutation importance figures
# -----------------------------------------------------------------------------
cat("Step 5: CV and importance plots\n")
source("modelling/plots/cv_pipeline_plots.R")
cat("\n")
source("modelling/plots/model_permutation_pruning_plots.R")
cat("\n")

# -----------------------------------------------------------------------------
# 6. Best GPR feature importance (standalone plot)
# -----------------------------------------------------------------------------
cat("Step 6: Best GPR feature importance plot\n")
source("modelling/plots/gpr_feature_importance_plot.R")
cat("\n")

# -----------------------------------------------------------------------------
# 7. GPR predictions (maps, PDPs, importance) -> figures/interpolation
# -----------------------------------------------------------------------------
cat("Step 7: GPR predictions and maps\n")
source("modelling/gpr/gpr_predictions.R")
cat("\n")

# -----------------------------------------------------------------------------
# 8. Region collage (spatial predictions by region) -> figures/interpolation
# -----------------------------------------------------------------------------
cat("Step 8: GPR region collage\n")
source("modelling/plots/gpr_region_collages.R")  # defines create_gpr_region_collage (may rm ls)
if (file.exists("figures/interpolation/gpr_prediction_grid.rds")) {
  pred_grid <- readRDS("figures/interpolation/gpr_prediction_grid.rds")
  collage <- create_gpr_region_collage(pred_grid)
  ggplot2::ggsave("figures/interpolation/gpr_region_collage.png", collage, width = 14, height = 10, dpi = 150)
  cat("Saved figures/interpolation/gpr_region_collage.png\n")
} else {
  cat("  (gpr_prediction_grid.rds not found; run Step 7 first or create collage manually from gpr_result$prediction_grid)\n")
}
cat("\n")

# -----------------------------------------------------------------------------
# 9. Supplement: region outlines map -> figures/supplement
# -----------------------------------------------------------------------------
cat("Step 9: Supplement – region outlines\n")
source("modelling/plots/plot_region_outlines.R")
cat("\n")

# -----------------------------------------------------------------------------
# Done
# -----------------------------------------------------------------------------
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("PAPER PIPELINE COMPLETE\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")
cat("Outputs:\n")
cat("  figures/cv_pipeline_output/  – CV results, pruning, tuning, importance plots\n")
cat("  figures/interpolation/       – GPR maps, PDPs, region collage\n")
cat("  figures/supplement/          – region_shapes.png\n")
cat("\n")
