# Spatial / lat-lon / region effect for all models (GPR, GAM, XGB)
#
# For each model in model_list, runs CV with configurations:
#   Env only, + latitude, + longitude, + lat+lon, + region (if present), + lat+lon+region
# Uses the same CV method (cv_type) and folds as the rest of the pipeline.
# Produces a bar chart of mean R² by configuration and model (like gpr_parameter_search_summary).
#
# Expects from .GlobalEnv: model_list, cv_type, n_folds, cv_blocksize, exclude_regions,
#   target_var, log_transform_target, use_shap_per_model (for load_model_vars logic in data prep)

setwd(here::here())
set.seed(42)

source("modelling/R/helpers.R")
source("modelling/R/extract_covariates_from_rasters.R")
source("modelling/R/assign_region_from_latlon.R")
load_packages(c("here", "dplyr", "ggplot2", "tidyr", "sf"))

target_var   <- get0("target_var", envir = .GlobalEnv, ifnotfound = "median_carbon_density_100cm")
model_list   <- get0("model_list", envir = .GlobalEnv, ifnotfound = c("GPR", "GAM", "XGB"))
n_folds      <- get0("n_folds", envir = .GlobalEnv, ifnotfound = 5L)
cv_type      <- get0("cv_type", envir = .GlobalEnv, ifnotfound = "spatial")
cv_blocksize <- get0("cv_blocksize", envir = .GlobalEnv, ifnotfound = 5000L)
exclude_regions <- get0("exclude_regions", envir = .GlobalEnv, ifnotfound = character(0))
log_response  <- isTRUE(get0("log_transform_target", envir = .GlobalEnv, ifnotfound = TRUE))
out_dir      <- "output/cv_pipeline"
fig_dir      <- out_dir
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dpi          <- get0("dpi", envir = .GlobalEnv, ifnotfound = 150)
show_titles  <- isTRUE(get0("show_titles", envir = .GlobalEnv, ifnotfound = TRUE))

cat("\n========================================\n")
cat("SPATIAL/CATEGORICAL EFFECT (ALL MODELS)\n")
cat("========================================\n\n")

# -----------------------------------------------------------------------------
# Data (core-level, same as cv_pipeline)
# -----------------------------------------------------------------------------
dat <- readr::read_rds("data/all_extracted_new.rds")
if (length(exclude_regions) > 0L) {
  if (!"region" %in% names(dat)) dat <- assign_region_from_latlon(dat)
  dat <- dat[is.na(dat$region) | !dat$region %in% exclude_regions, ]
}
predictor_vars <- raster_covariates[raster_covariates %in% colnames(dat)]
if (length(predictor_vars) == 0) stop("No predictor_vars in dat.")
need_cols <- c("longitude", "latitude", target_var, predictor_vars)
if ("region" %in% names(dat)) need_cols <- c(need_cols, "region")
complete_dat <- dat %>%
  dplyr::select(dplyr::all_of(intersect(need_cols, names(dat)))) %>%
  dplyr::filter(complete.cases(.))
complete_dat$median_carbon_density <- complete_dat[[target_var]]
predictor_vars <- predictor_vars[predictor_vars %in% colnames(complete_dat)]
has_region <- "region" %in% names(complete_dat)
n_cores <- nrow(complete_dat)
cat("Cores:", n_cores, " | region in data:", has_region, "\n")

# Per-model base (env-only) predictor sets: prefer SHAP when use_shap_per_model
cov_dir            <- "output/covariate_selection"
use_shap_per_model <- isTRUE(get0("use_shap_per_model", envir = .GlobalEnv, ifnotfound = FALSE))
shap_combined      <- file.path(cov_dir, "pruned_model_variables_shap.csv")
perm_combined      <- file.path(cov_dir, "pruned_model_variables_perm.csv")
src_file           <- if (use_shap_per_model && file.exists(shap_combined)) shap_combined else if (file.exists(perm_combined)) perm_combined else NULL
per_model_vars     <- list()
if (!is.null(src_file)) {
  df <- read.csv(src_file, stringsAsFactors = FALSE)
  if (all(c("model", "variable") %in% names(df))) {
    for (m in unique(df$model)) {
      v <- intersect(df$variable[df$model == m], colnames(complete_dat))
      if (length(v) >= 2) per_model_vars[[m]] <- v
    }
  }
}
# Fallback: pruned_variables_to_include_<model>.csv
for (m in model_list) {
  if (m %in% names(per_model_vars)) next
  path <- file.path(cov_dir, paste0("pruned_variables_to_include_", tolower(m), ".csv"))
  if (file.exists(path)) {
    v <- intersect(read.csv(path, stringsAsFactors = FALSE)$variable, colnames(complete_dat))
    if (length(v) >= 2) per_model_vars[[m]] <- v
  }
}
if (length(per_model_vars) == 0) {
  shared <- file.path(cov_dir, "pruned_variables_to_include.csv")
  if (file.exists(shared))
    per_model_vars <- setNames(rep(list(intersect(read.csv(shared, stringsAsFactors = FALSE)$variable, colnames(complete_dat))), length(model_list)), model_list)
  else
    per_model_vars <- setNames(rep(list(predictor_vars), length(model_list)), model_list)
}
# Env-only: drop lat, lon, region for "Env only" base
env_only_by_model <- lapply(per_model_vars, function(v) setdiff(v, c("latitude", "longitude", "region")))

# -----------------------------------------------------------------------------
# Folds (same as pipeline: cv_type)
# -----------------------------------------------------------------------------
if (identical(cv_type, "spatial") && length(cv_blocksize) > 0) {
  fold_info <- get_cached_spatial_folds(
    dat = complete_dat, block_size = cv_blocksize, n_folds = n_folds,
    cache_tag = "spatial_categorical_effect", exclude_regions = exclude_regions, progress = TRUE
  )
  fold_indices <- fold_info$fold_indices
  cv_method_name <- paste0("spatial_block_", cv_blocksize, "m")
  cat("Using spatial folds (", cv_blocksize, " m).\n", sep = "")
} else {
  fold_indices <- sample(rep(seq_len(n_folds), length.out = n_cores))
  cv_method_name <- "random_split"
  cat("Using random folds.\n")
}
core_sf <- sf::st_as_sf(complete_dat, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)

# -----------------------------------------------------------------------------
# Configurations: Env only, + lat, + lon, + lat+lon, + region, + lat+lon+region
# -----------------------------------------------------------------------------
configs <- list("Env only" = NULL)  # NULL = use env_only_by_model per model
configs[["+ latitude"]]     <- "latitude"
configs[["+ longitude"]]    <- "longitude"
configs[["+ lat+lon"]]      <- c("latitude", "longitude")
if (has_region) {
  configs[["+ region"]]       <- "region"
  configs[["+ lat+lon+region"]] <- c("latitude", "longitude", "region")
}

results_list <- list()
for (model_name in intersect(model_list, c("GPR", "GAM", "XGB"))) {
  base_vars <- env_only_by_model[[model_name]]
  if (is.null(base_vars) || length(base_vars) < 2) base_vars <- predictor_vars
  for (config_name in names(configs)) {
    extra <- configs[[config_name]]
    pvars <- if (is.null(extra)) base_vars else c(extra, base_vars)
    pvars <- unique(pvars[pvars %in% colnames(complete_dat)])
    if (length(pvars) < 2) next
    cat("  ", model_name, " — ", config_name, " (", length(pvars), " vars) ... ", sep = "")
    res <- tryCatch(
      run_cv(
        cv_method_name = cv_method_name,
        fold_indices   = fold_indices,
        core_data      = complete_dat,
        predictor_vars = pvars,
        models         = model_name,
        verbose        = FALSE,
        core_sf        = core_sf,
        log_response   = log_response
      ),
      error = function(e) { cat("error:", e$message, "\n"); NULL }
    )
    if (is.null(res) || nrow(res) == 0) { cat("no results\n"); next }
    agg <- res %>% dplyr::group_by(model) %>% dplyr::summarise(
      mean_r2 = mean(r2, na.rm = TRUE), sd_r2 = sd(r2, na.rm = TRUE),
      mean_rmse = mean(rmse, na.rm = TRUE), sd_rmse = sd(rmse, na.rm = TRUE),
      n_folds_valid = sum(!is.na(r2)), .groups = "drop"
    )
    if (nrow(agg) > 0) {
      results_list[[length(results_list) + 1]] <- data.frame(
        model = model_name,
        configuration = config_name,
        mean_r2 = agg$mean_r2[1],
        sd_r2 = agg$sd_r2[1],
        mean_rmse = agg$mean_rmse[1],
        sd_rmse = agg$sd_rmse[1],
        n_folds_valid = agg$n_folds_valid[1],
        stringsAsFactors = FALSE
      )
      cat("R² =", round(agg$mean_r2[1], 4), "\n")
    } else cat("no results\n")
  }
}

if (length(results_list) == 0) {
  cat("No results. Check data and model list.\n")
} else {
  summary_tab <- dplyr::bind_rows(results_list)
  write.csv(summary_tab, file.path(out_dir, "spatial_categorical_effect_all_models.csv"), row.names = FALSE)
  cat("Saved:", file.path(out_dir, "spatial_categorical_effect_all_models.csv"), "\n")

  # Bar chart: configuration vs mean_r2, faceted by model (or grouped)
  summary_tab <- summary_tab %>%
    dplyr::mutate(
      configuration = factor(configuration, levels = names(configs)),
      model = factor(model, levels = intersect(c("GPR", "GAM", "XGB"), model_list))
    )
  p <- ggplot(summary_tab, aes(x = configuration, y = mean_r2, fill = model)) +
    geom_col(position = position_dodge(0.9), alpha = 0.9) +
    geom_errorbar(aes(ymin = mean_r2 - sd_r2, ymax = mean_r2 + sd_r2),
                  position = position_dodge(0.9), width = 0.2, linewidth = 0.4, color = "gray20") +
    scale_fill_manual(values = c(GPR = "#e41a1c", GAM = "#4daf4a", XGB = "#377eb8"), name = "Model") +
    labs(
      x = "",
      y = paste0("Mean R² (", cv_method_name, ")"),
      title = if (show_titles) "Effect of latitude, longitude, and region on CV R² (all models)" else NULL,
      subtitle = if (show_titles) "Same folds and env predictors per model; adding spatial terms can lower R² under spatial CV." else NULL
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 9, color = "gray40"),
      legend.position = "right",
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  ggsave(file.path(out_dir, "spatial_categorical_effect_all_models.png"), p, width = 9, height = 5, dpi = dpi)
  if (fig_dir != out_dir) ggsave(file.path(fig_dir, "spatial_categorical_effect_all_models.png"), p, width = 9, height = 5, dpi = dpi)
  cat("Figure saved: spatial_categorical_effect_all_models.png\n")
}

cat("\nDone.\n")
