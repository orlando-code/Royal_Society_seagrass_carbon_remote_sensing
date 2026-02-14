# =============================================================================
# GPR feature importance: best hyperparameters + GPR-pruned covariate set
# Fits the best-performing GPR (from gpr_best_config.rds) on the GPR-pruned
# predictor set (from model permutation pruning or pruned_variables_to_include_gpr.csv),
# computes permutation importance, and saves one clear plot + CSV.
# =============================================================================

setwd(here::here())
set.seed(42)

source("modelling/R/helpers.R")
source("modelling/R/assign_region_from_latlon.R")
source("modelling/R/plot_config.R")
load_packages(c("here", "tidyverse", "ggplot2", "GauPro"))

# Config
target_var <- "median_carbon_density_100cm"
include_species <- FALSE
include_region <- FALSE
include_spatial_smoothing <- TRUE
n_permutations <- 5
n_val <- 500

# Output
out_dir <- "figures/cv_pipeline_output"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (exists("dpi", envir = .GlobalEnv)) dpi <- get("dpi", envir = .GlobalEnv) else dpi <- 150

# Best config path: check figures then modelling
best_config_path <- "figures/cv_pipeline_output/gpr_best_config.rds"
if (!file.exists(best_config_path)) best_config_path <- "modelling/cv_pipeline_output/gpr_best_config.rds"
if (!file.exists(best_config_path)) stop("Run gpr_hyperparameter_tuning.R first to create gpr_best_config.rds")

cat("Loading best GPR config from:", best_config_path, "\n")
best_config_list <- readRDS(best_config_path)
gpr_kernel <- best_config_list$kernel
cat("Kernel:", gpr_kernel, "\n\n")

# Data
cat("Loading data...\n")
raster_covariates_dat <- readRDS("data/all_extracted_new.rds")
raster_covariates_dat <- process_rs_covariates(raster_covariates_dat)
if (!"region" %in% names(raster_covariates_dat) && all(c("longitude", "latitude") %in% names(raster_covariates_dat))) {
  raster_covariates_dat <- assign_region_from_latlon(raster_covariates_dat)
}
data_cols <- names(raster_covariates_dat)

# GPR-pruned covariate set (same as gpr_predictions.R)
gam_pruned_vars <- get_pruned_predictors_for_model("GPR", data_cols = data_cols)
if (is.null(gam_pruned_vars) || length(gam_pruned_vars) < 2L) {
  pruned_file_gpr <- "figures/covariate_selection/pruned_variables_to_include_gpr.csv"
  pruned_file <- "figures/covariate_selection/pruned_variables_to_include.csv"
  if (file.exists(pruned_file_gpr)) {
    gam_pruned_vars <- read_csv(pruned_file_gpr, show_col_types = FALSE)$variable
  } else if (file.exists(pruned_file)) {
    gam_pruned_vars <- read_csv(pruned_file, show_col_types = FALSE)$variable
  } else {
    stop("Pruned variables not found. Run model_permutation_pruning.R or covariate_pruning scripts.")
  }
}
gam_pruned_vars <- gam_pruned_vars[gam_pruned_vars %in% data_cols]
cat("GPR-pruned predictor set (", length(gam_pruned_vars), " vars): ", paste(gam_pruned_vars, collapse = ", "), "\n\n", sep = "")

predictor_vars <- gam_pruned_vars
if (include_species && "seagrass_species" %in% data_cols) predictor_vars <- c(predictor_vars, "seagrass_species")
if (include_region && "Region" %in% data_cols) predictor_vars <- c(predictor_vars, "Region")
predictor_vars <- setdiff(predictor_vars, setdiff(predictor_vars, data_cols))
if (!target_var %in% data_cols) stop("Target variable not in data.")
if (!all(c("longitude", "latitude") %in% data_cols)) stop("longitude/latitude not in data.")

# Minimal prediction grid (one complete row) so fit_gaussian_process_regression can run
need_cols <- unique(c("longitude", "latitude", predictor_vars))
ok <- complete.cases(raster_covariates_dat[, c(target_var, need_cols), drop = FALSE])
if (!any(ok)) stop("No complete cases for target and predictors.")
minimal_grid <- raster_covariates_dat[which(ok)[1], need_cols, drop = FALSE]

source("modelling/R/gpr_funs.R")
pruned_formula_reg <- as.formula(paste(target_var, "~", paste(predictor_vars, collapse = " + ")))

cat("Fitting GPR (best kernel + GPR-pruned set)...\n")
gpr_result <- fit_gaussian_process_regression(
  dat = raster_covariates_dat,
  value_var = target_var,
  coords = c("longitude", "latitude"),
  predictor_vars = predictor_vars,
  formula = pruned_formula_reg,
  prediction_grid = minimal_grid,
  include_spatial = include_spatial_smoothing,
  kernel = gpr_kernel
)
cat("Fitted. n_train =", gpr_result$n_train, "\n\n")

# Permutation importance
cat("Computing permutation importance (n_permutations =", n_permutations, ", n_val =", n_val, ")...\n")
train_data_imp <- raster_covariates_dat[, c(target_var, predictor_vars), drop = FALSE]
train_data_imp <- train_data_imp[complete.cases(train_data_imp), ]
importance_results <- gpr_permutation_importance(
  gpr_result$model, train_data_imp, target_var, predictor_vars,
  n_permutations = n_permutations, n_val = n_val,
  scale_params = gpr_result$scale_params
)
cat("Done.\n\n")

# Save CSV
write.csv(importance_results, file.path(out_dir, "gpr_feature_importance_best_model.csv"), row.names = FALSE)
cat("Saved:", file.path(out_dir, "gpr_feature_importance_best_model.csv"), "\n")

# Plot: best fitted model importance; sqrt scale so all 14 bars are visible (not dominated by top 2-3)
importance_results$variable_label <- label_vars(importance_results$variable)
p <- ggplot(importance_results, aes(x = reorder(variable_label, rmse_increase), y = rmse_increase)) +
  geom_col(fill = "steelblue", alpha = 0.85) +
  coord_flip() +
  scale_y_sqrt() +
  labs(
    x = NULL,
    y = "RMSE increase (permutation importance; sqrt scale)"
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.y = element_text(size = rel(0.9)))
print(p)
ggsave(file.path(out_dir, "gpr_feature_importance_best_model.png"), p,
       width = 9, height = max(5, nrow(importance_results) * 0.35), dpi = dpi)
cat("Saved:", file.path(out_dir, "gpr_feature_importance_best_model.png"), " (best fitted GPR; sqrt scale)\n")
