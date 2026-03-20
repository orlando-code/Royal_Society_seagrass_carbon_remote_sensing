# SHAP (Shapley) feature importance for each model after hyperparameter tuning
#
# Uses the best hyperparameters and best covariate set per model (from
# hyperparameter_tuning_pipeline and pruned_model_variables_shap/perm).
# Fits each model with best config on full data and computes mean |phi| per
# variable via iml::Shapley, then saves CSV and bar plot per model.
#
# Requires: hyperparameter_tuning_pipeline.R has been run (best_config_*.rds),
#           pruned_model_variables_shap.csv or pruned_model_variables_perm.csv
# Outputs:  output/<cv_regime>/cv_pipeline/importance_shap_<model>.csv and .png
#
# Usage: source after step 3 in run_paper.R, or run standalone.

setwd(here::here())
set.seed(42)

source("modelling/R/helpers.R")
source("modelling/R/extract_covariates_from_rasters.R")
source("modelling/R/assign_region_from_latlon.R")
load_packages(c("here", "mgcv", "tidyverse", "ggplot2", "GauPro", "xgboost", "sf", "iml"))

cv_out  <- get0("cv_output_dir", envir = .GlobalEnv, ifnotfound = "output")
out_dir <- file.path(cv_out, "cv_pipeline")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cov_dir    <- file.path(cv_out, "covariate_selection")
config_dir <- out_dir

target_var <- "median_carbon_density_100cm"
log_response <- isTRUE(get0("log_transform_target", envir = .GlobalEnv, ifnotfound = TRUE))
exclude_regions <- get0("exclude_regions", envir = .GlobalEnv, ifnotfound = character(0))
model_list <- get0("model_list", envir = .GlobalEnv, ifnotfound = c("GPR", "GAM", "XGB"))
n_shap_points <- as.integer(get0("n_shap_points", envir = .GlobalEnv, ifnotfound = 30L))
dpi <- get0("dpi", envir = .GlobalEnv, ifnotfound = 150)

# Data (same as permutation_importance_final / fit_final_models)
dat <- readr::read_rds("data/all_extracted_new.rds")
if (length(exclude_regions) > 0L) {
  if (!"region" %in% names(dat)) dat <- assign_region_from_latlon(dat)
  dat <- dat[is.na(dat$region) | !dat$region %in% exclude_regions, ]
}
predictor_vars_all <- raster_covariates[raster_covariates %in% colnames(dat)]
core_data <- dat %>%
  dplyr::mutate(median_carbon_density = .data[[target_var]]) %>%
  dplyr::select(
    longitude, latitude, median_carbon_density,
    dplyr::all_of(predictor_vars_all),
    dplyr::all_of(intersect(c("seagrass_species", "region"), names(dat)))
  ) %>%
  dplyr::filter(complete.cases(.))
core_data <- as.data.frame(core_data)

# Per-model covariate set
use_shap_per_model <- isTRUE(get0("use_shap_per_model", envir = .GlobalEnv, ifnotfound = FALSE))
per_model_vars <- get_per_model_vars(cov_dir, colnames(core_data), use_shap_first = use_shap_per_model)

# Helper: load best config (same as permutation_importance_final)
load_best_config <- function(model_name) {
  base <- file.path(config_dir, "best_config")
  paths <- switch(model_name,
    XGB = c(paste0(base, "_xgb.rds"), file.path(config_dir, "xgb_best_config.rds")),
    GAM = c(paste0(base, "_gam.rds"), file.path(config_dir, "gam_best_config.rds")),
    GPR = c(paste0(base, "_gpr.rds"), file.path(config_dir, "gpr_best_config.rds")),
    NULL
  )
  if (is.null(paths)) return(NULL)
  for (p in paths) if (file.exists(p)) return(readRDS(p))
  NULL
}

cat("\n========================================\n")
cat("SHAP IMPORTANCE (best config + best vars per model)\n")
cat("========================================\n\n")

for (model_name in model_list) {
  pvars <- tryCatch(load_model_vars(model_name, per_model_vars, use_shap_first = use_shap_per_model), error = function(e) NULL)
  if (is.null(pvars) || length(pvars) < 2L) {
    cat("Skip", model_name, ": no predictor set.\n")
    next
  }
  hp <- load_best_config(model_name)
  if (model_name == "GPR" && is.list(hp)) {
    hp <- list(kernel = hp$kernel, nug.min = hp$nug.min, nug.max = hp$nug.max, nug.est = TRUE)
  }
  if (model_name == "XGB" && is.list(hp)) {
    hp <- list(nrounds = hp$nrounds, max_depth = hp$max_depth, learning_rate = hp$learning_rate %||% 0.1,
               subsample = hp$subsample %||% 0.8, colsample_bytree = hp$colsample_bytree %||% 0.8,
               min_child_weight = hp$min_child_weight %||% 1L,
               min_split_loss = hp$min_split_loss %||% 0,
               reg_reg_lambda = hp$reg_reg_lambda %||% 1)
  }
  if (model_name == "GAM" && is.list(hp)) {
    hp <- list(k_covariate = hp$k_covariate %||% hp$k_spatial %||% 6L)
  }
  cat("Model:", model_name, "| Predictors:", length(pvars), "| Hyperparams:", if (length(hp)) "from tuning" else "defaults", "| n_points:", n_shap_points, "\n")
  imp <- tryCatch(
    compute_shap_importance(
      core_data, "median_carbon_density", pvars, model_name,
      n_points = n_shap_points, log_response = log_response, hyperparams = hp
    ),
    error = function(e) { cat("  Error:", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(imp) || nrow(imp) == 0L) next
  imp$model <- model_name
  write.csv(imp, file.path(out_dir, paste0("importance_shap_", tolower(model_name), ".csv")), row.names = FALSE)
  cat("  Saved", paste0("importance_shap_", tolower(model_name), ".csv"), "\n")
  p <- ggplot(imp, aes(x = reorder(variable, shap_importance), y = shap_importance)) +
    geom_col(fill = "darkgreen", alpha = 0.85) +
    coord_flip() +
    labs(x = NULL, y = "Mean |SHAP| (Shapley importance)", title = paste(model_name, "SHAP importance")) +
    theme_minimal(base_size = 12)
  ggsave(file.path(out_dir, paste0("importance_shap_", tolower(model_name), ".png")), p,
         width = 9, height = max(5, nrow(imp) * 0.35), dpi = dpi)
  cat("  Saved", paste0("importance_shap_", tolower(model_name), ".png"), "\n\n")
}

cat("SHAP importance complete. Outputs in", out_dir, "\n")
