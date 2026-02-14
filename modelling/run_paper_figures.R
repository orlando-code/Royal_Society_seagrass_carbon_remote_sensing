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
cat("\t\tStep -1: Configuring paper pipeline ...\n")

# Formatting
dpi <- 150
show_titles <- TRUE
assign("dpi", dpi, envir = .GlobalEnv)
assign("show_titles", show_titles, envir = .GlobalEnv)

# Pruning configuration (run before CV pipeline)
use_correlation_filter <- TRUE
correlation_filter_threshold <- 0.8
use_nested_selection <- FALSE
nested_selection_max_vars <- 15
model_list <- c("GPR", "GAM", "XGB") # Quick test: only GPR, GAM, XGBoost
n_folds <- 5

assign("model_list", model_list, envir = .GlobalEnv)
assign("use_correlation_filter", use_correlation_filter, envir = .GlobalEnv)
assign("correlation_filter_threshold", correlation_filter_threshold, envir = .GlobalEnv)
assign("use_nested_selection", use_nested_selection, envir = .GlobalEnv)
assign("nested_selection_max_vars", nested_selection_max_vars, envir = .GlobalEnv)
assign("n_folds", n_folds, envir = .GlobalEnv)
if (exists("raster_covariates")) assign("raster_covariates", raster_covariates, envir = .GlobalEnv)

# generate the raster dir information and ensure globals
source("modelling/R/extract_covariates_from_rasters.R")
cat("\n")
cat("  use_correlation_filter:", use_correlation_filter, "\n")
cat("  use_nested_selection:", use_nested_selection, "\n")
cat("  correlation_filter_threshold:", correlation_filter_threshold, "\n")
cat("  nested_selection_max_vars:", nested_selection_max_vars, "\n")
cat("  model_list:", paste(model_list, collapse = ", "), "\n")
cat("  n_folds:", n_folds, "\n")
cat("  dpi:", dpi, "\n")
cat("  show_titles:", show_titles, "\n")
cat("\n")


# -----------------------------------------------------------------------------
# 0. Data (build if missing)
# -----------------------------------------------------------------------------
if (!file.exists("data/all_extracted_new.rds")) {
  cat("\t\tStep 0: Building all_extracted_new.rds ...\n")
  source("modelling/pipeline/build_all_extracted_new.R")
  cat("\n")
} else {
  cat("\t\tStep 0: Using existing data/all_extracted_new.rds\n\n")
}


# -----------------------------------------------------------------------------
# 1. Covariate pruning (correlation and/or nested selection)
# -----------------------------------------------------------------------------
cat("\t\tStep 1: Covariate pruning (correlation and/or nested selection)\n")
source("modelling/pipeline/covariate_pruning_pipeline.R")
cat("\n")
# TODO: look further into model-specific and/or generic pruning methods to allow direct comparison between model performance

# -----------------------------------------------------------------------------
# 2. CV pipeline (folds, CV results, inline plots) -> figures/cv_pipeline_output
# -----------------------------------------------------------------------------
cat("\t\tStep 2: CV pipeline (spatial folds, comparison, results)\n")
source("modelling/pipeline/cv_pipeline.R")
cat("\n")

# -----------------------------------------------------------------------------
# 3. GPR hyperparameter tuning -> figures/cv_pipeline_output
# -----------------------------------------------------------------------------
cat("\t\tStep 3: GPR hyperparameter tuning\n")
source("modelling/pipeline/gpr_hyperparameter_tuning.R")
cat("\n")

# -----------------------------------------------------------------------------
# 4. GPR categorical and spatial investigations (feed into parameter summary plot)
# -----------------------------------------------------------------------------
cat("\t\tStep 4a: GPR categorical investigation (species/region)\n")
source("modelling/gpr/gpr_categorical_investigation.R")
cat("\n")
cat("\t\tStep 4b: GPR spatial terms investigation (lat/lon, region)\n")
source("modelling/gpr/gpr_spatial_investigation.R")
cat("\n")
cat("\t\tStep 4c: GPR parameter search summary (tuning plot)\n")
source("modelling/gpr/gpr_parameter_search_summary.R")
cat("\n")

# -----------------------------------------------------------------------------
# 5. Plot-only: CV figures and permutation importance figures
# -----------------------------------------------------------------------------
cat("\t\tStep 5: CV and importance plots\n")
source("modelling/plots/cv_pipeline_plots.R")
cat("\n")
# source("modelling/plots/model_permutation_pruning_plots.R")
  # cat("\n")# TODO: this will depend on feature pruning regime.
# -----------------------------------------------------------------------------
# 6. Best GPR feature importance (standalone plot)
# -----------------------------------------------------------------------------
cat("\t\tStep 6: Best GPR feature importance plot\n")
source("modelling/plots/gpr_feature_importance_plot.R")
cat("\n")

# -----------------------------------------------------------------------------
# 7. GPR predictions (maps, PDPs, importance) -> figures/interpolation
# -----------------------------------------------------------------------------
cat("\t\tStep 7: GPR predictions and maps\n")
source("modelling/gpr/gpr_predictions.R")
cat("\n")

# -----------------------------------------------------------------------------
# 8. Region collage (spatial predictions by region) -> figures/interpolation
# -----------------------------------------------------------------------------
cat("\t\tStep 8: GPR region collage\n")
source("modelling/plots/gpr_region_collages.R") # defines create_gpr_region_collage (may rm ls)
if (file.exists("figures/interpolation/gpr_prediction_grid.rds")) {
  pred_grid <- readRDS("figures/interpolation/gpr_prediction_grid.rds")
  collage <- create_gpr_region_collage(pred_grid)
  ggplot2::ggsave("figures/interpolation/gpr_region_collage.png", collage, width = 14, height = 10, dpi = dpi)
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
cat("  figures/covariate_selection/ – Covariate pruning results\n")
cat("  figures/cv_pipeline_output/  – CV results, pruning, tuning, importance plots\n")
cat("  figures/predictions/       – GPR maps, PDPs, region collage\n")
cat("  figures/supplement/          – region_shapes.png\n")
cat("\n")
