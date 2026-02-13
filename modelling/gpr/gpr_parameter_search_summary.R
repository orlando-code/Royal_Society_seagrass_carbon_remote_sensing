# GPR parameter search: sequential investigation, single bar chart
#
# Order: 1. Nugget → 2. Kernel → 3. Spatial → 4. Categorical (regularization first, then model expansion).
# Each step uses the best from the previous; baseline is carried forward so the plot never shows a drop.
#
# Why can R² *decrease* when adding spatial or species?
# - Spatial block CV is strict: test blocks are spatially separate, so we're predicting in
#   new locations. Adding lat/lon or region can overfit to training geography and generalize
#   worse to held-out blocks. More predictors also increase variance of the CV estimate.
# - Species/region can be redundant with env (e.g. env already captures regional differences),
#   or add noise in folds where a level is rare in the test set. So "more information"
#   does not guarantee higher R² in spatial CV.

rm(list = ls())
setwd(here::here())
set.seed(42)

source("modelling/R/helpers.R")
source("modelling/R/extract_covariates_from_rasters.R")
source("modelling/R/gpr_funs.R")
load_packages(c("here", "tidyverse", "ggplot2", "dplyr"))

target_var <- "median_carbon_density_100cm"
out_dir <- "figures/cv_pipeline_output"
fig_dir <- out_dir
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

n_folds <- 5
kernels_to_test <- c("matern52", "matern32", "gaussian")

cat("\n========================================\n")
cat("GPR PARAMETER SEARCH (SEQUENTIAL)\n")
cat("========================================\n\n")

# -----------------------------------------------------------------------------
# Load data and env vars
# -----------------------------------------------------------------------------
dat <- readr::read_rds("data/all_extracted_new.rds")
dat <- process_rs_covariates(dat)

# Prefer generic model permutation pruning for GPR, then GPR-specific, then GAM list
env_vars <- get_pruned_predictors_for_model("GPR", data_cols = names(dat))
if (is.null(env_vars) || length(env_vars) < 2L) {
  pruned_file_gpr <- "figures/covariate_selection/pruned_variables_to_include_gpr.csv"
  pruned_file <- "figures/covariate_selection/pruned_variables_to_include.csv"
  if (file.exists(pruned_file_gpr)) {
    env_vars <- readr::read_csv(pruned_file_gpr, show_col_types = FALSE)$variable
  } else if (file.exists(pruned_file)) {
    env_vars <- readr::read_csv(pruned_file, show_col_types = FALSE)$variable
  } else {
    stop("Pruned variables not found. Run model_permutation_pruning.R or covariate pruning first.")
  }
}
env_vars <- env_vars[env_vars %in% names(dat)]

has_region <- "region" %in% names(dat)
has_species <- "seagrass_species" %in% names(dat)
need_vars <- c("longitude", "latitude", target_var, env_vars)
if (has_species) need_vars <- c(need_vars, "seagrass_species")
if (has_region) need_vars <- c(need_vars, "region")

dat_complete <- as.data.frame(dat[, intersect(need_vars, names(dat)), drop = FALSE])
dat_complete <- dat_complete[complete.cases(dat_complete), ]

cat("Complete cases:", nrow(dat_complete), "\n")
cat("Env predictors:", length(env_vars), " | region:", has_region, " | Species:", has_species, "\n\n")

# -----------------------------------------------------------------------------
# Load 1km spatial folds (both cache locations)
# -----------------------------------------------------------------------------
folds_ids <- NULL
for (base in c("figures/cv_pipeline_output", "modelling/cv_pipeline_output")) {
  candidate <- file.path(base, "spatial_folds_cache.rds")
  if (file.exists(candidate)) {
    cached <- tryCatch(readRDS(candidate), error = function(e) NULL)
    if (!is.null(cached) && "spatial_strategies" %in% names(cached)) {
      for (strategy in cached$spatial_strategies) {
        if (strategy$method == "spatial_block_1000m") {
          if (length(strategy$folds) == nrow(dat_complete)) {
            folds_ids <- strategy$folds
          } else if ("random_core_variable" %in% names(dat)) {
            dat_full <- readr::read_rds("data/all_extracted_new.rds")
            dat_full <- process_rs_covariates(dat_full)
            core_data_cv <- dat_full %>%
              dplyr::group_by(random_core_variable) %>%
              dplyr::summarise(longitude = dplyr::first(longitude), latitude = dplyr::first(latitude), .groups = "drop") %>%
              dplyr::filter(!is.na(longitude) & !is.na(latitude)) %>%
              dplyr::arrange(random_core_variable)
            if (nrow(core_data_cv) == length(strategy$folds)) {
              core_locations <- paste(round(core_data_cv$longitude, 8), round(core_data_cv$latitude, 8), sep = "_")
              fold_lookup <- data.frame(location = core_locations, fold = strategy$folds,
                                       longitude = core_data_cv$longitude, latitude = core_data_cv$latitude, stringsAsFactors = FALSE)
              point_locations <- paste(round(dat_complete$longitude, 8), round(dat_complete$latitude, 8), sep = "_")
              folds_ids <- fold_lookup$fold[match(point_locations, fold_lookup$location)]
              unmatched <- is.na(folds_ids)
              if (sum(unmatched) > 0) {
                for (i in which(unmatched)) {
                  lon <- dat_complete$longitude[i]; lat <- dat_complete$latitude[i]
                  dists <- sqrt((fold_lookup$longitude - lon)^2 + (fold_lookup$latitude - lat)^2)
                  folds_ids[i] <- fold_lookup$fold[which.min(dists)]
                }
              }
            }
          }
          if (!is.null(folds_ids)) break
        }
      }
    }
    if (!is.null(folds_ids)) break
  }
}
if (is.null(folds_ids)) {
  for (base in c("modelling/cv_pipeline_output", "figures/cv_pipeline_output")) {
    path_1km <- file.path(base, "spatial_folds_1km_cache.rds")
    if (file.exists(path_1km)) {
      spatial_cv <- tryCatch(readRDS(path_1km), error = function(e) NULL)
      if (!is.null(spatial_cv) && "folds_ids" %in% names(spatial_cv) && length(spatial_cv$folds_ids) == nrow(dat_complete)) {
        folds_ids <- spatial_cv$folds_ids
        break
      }
    }
  }
}
if (is.null(folds_ids)) {
  stop("1km spatial folds not found. Run gpr_hyperparameter_tuning.R or cv_pipeline.R first.")
}
cat("Using", length(unique(folds_ids)), "spatial folds.\n\n")

# Helper: run 5-fold CV for one predictor set and kernel (GPR parameter search only). Not the same as helpers::run_cv.
run_gpr_cv_parameter_search <- function(predictor_vars, kernel, nug.min = 1e-8, nug.est = TRUE) {
  pv <- intersect(predictor_vars, names(dat_complete))
  if (length(pv) == 0) return(list(mean_r2 = NA, sd_r2 = NA, mean_rmse = NA, sd_rmse = NA))
  r2_vec <- numeric(n_folds)
  rmse_vec <- numeric(n_folds)
  for (fold in 1:n_folds) {
    train_data <- dat_complete[folds_ids != fold, ]
    test_data <- dat_complete[folds_ids == fold, ]
    res <- fit_gpr_cv(train_data, test_data, target_var, pv, kernel = kernel, nug.min = nug.min, nug.est = nug.est)
    r2_vec[fold] <- res$r2
    rmse_vec[fold] <- res$rmse
  }
  valid <- !is.na(r2_vec)
  list(
    mean_r2 = mean(r2_vec[valid]),
    sd_r2 = if (sum(valid) > 1) sd(r2_vec[valid]) else NA_real_,
    mean_rmse = mean(rmse_vec[valid]),
    sd_rmse = if (sum(valid) > 1) sd(rmse_vec[valid]) else NA_real_,
    n_folds_valid = sum(valid)
  )
}

results_list <- list()
default_kernel <- "matern52"   # used for step 1 (nugget) only
nugget_mins <- c(1e-8, 1e-6, 1e-4)

# -----------------------------------------------------------------------------
# Step 1: Nugget (env only, default kernel; set regularization before expanding model)
# -----------------------------------------------------------------------------
cat("Step 1: Nugget (env only, kernel = ", default_kernel, ")...\n", sep = "")
step1_results <- list()
for (nug_min in nugget_mins) {
  out <- run_gpr_cv_parameter_search(env_vars, default_kernel, nug.min = nug_min, nug.est = TRUE)
  step1_results[[as.character(nug_min)]] <- out
  results_list[[length(results_list) + 1]] <- data.frame(
    step = "1. Nugget",
    configuration = paste0("nug.min = ", format(nug_min, scientific = TRUE)),
    mean_r2 = out$mean_r2,
    sd_r2 = out$sd_r2,
    mean_rmse = out$mean_rmse,
    sd_rmse = out$sd_rmse,
    n_folds_valid = out$n_folds_valid,
    stringsAsFactors = FALSE
  )
  cat("  nug.min ", format(nug_min, scientific = TRUE), ": R² =", round(out$mean_r2, 4), "\n")
}
best_nug_min <- nugget_mins[which.max(sapply(step1_results, function(x) x$mean_r2))]
cat("  Best nug.min:", format(best_nug_min, scientific = TRUE), "\n\n")

# -----------------------------------------------------------------------------
# Step 2: Kernel (carry forward step 1 best, then try other kernels; all use best_nug_min)
# -----------------------------------------------------------------------------
cat("Step 2: Kernel (nug.min = ", format(best_nug_min, scientific = TRUE), ")...\n", sep = "")
best_step1 <- step1_results[[as.character(best_nug_min)]]
results_list[[length(results_list) + 1]] <- data.frame(
  step = "2. Kernel",
  configuration = "Env only",
  mean_r2 = best_step1$mean_r2,
  sd_r2 = best_step1$sd_r2,
  mean_rmse = best_step1$mean_rmse,
  sd_rmse = best_step1$sd_rmse,
  n_folds_valid = best_step1$n_folds_valid,
  stringsAsFactors = FALSE
)
cat("  Env only (from step 1): R² =", round(best_step1$mean_r2, 4), "\n")

step2_results <- list("Env only" = best_step1)
for (k in kernels_to_test) {
  if (k == default_kernel) next   # already have Env only with default kernel
  out <- run_gpr_cv_parameter_search(env_vars, k, nug.min = best_nug_min, nug.est = TRUE)
  step2_results[[k]] <- out
  results_list[[length(results_list) + 1]] <- data.frame(
    step = "2. Kernel",
    configuration = k,
    mean_r2 = out$mean_r2,
    sd_r2 = out$sd_r2,
    mean_rmse = out$mean_rmse,
    sd_rmse = out$sd_rmse,
    n_folds_valid = out$n_folds_valid,
    stringsAsFactors = FALSE
  )
  cat("  ", k, ": R² =", round(out$mean_r2, 4), "\n")
}
best_step2_name <- names(step2_results)[which.max(sapply(step2_results, function(x) x$mean_r2))]
best_kernel <- if (best_step2_name == "Env only") default_kernel else best_step2_name
best_step2 <- step2_results[[best_step2_name]]
cat("  Best kernel:", best_kernel, "\n\n")

# -----------------------------------------------------------------------------
# Step 3: Spatial (best kernel + best nugget; carry forward, then + lat/lon, + region, etc.)
# -----------------------------------------------------------------------------
cat("Step 3: Spatial (kernel = ", best_kernel, ")...\n", sep = "")
results_list[[length(results_list) + 1]] <- data.frame(
  step = "3. Spatial",
  configuration = "Env only",
  mean_r2 = best_step2$mean_r2,
  sd_r2 = best_step2$sd_r2,
  mean_rmse = best_step2$mean_rmse,
  sd_rmse = best_step2$sd_rmse,
  n_folds_valid = best_step2$n_folds_valid,
  stringsAsFactors = FALSE
)
cat("  Env only (from step 2): R² =", round(best_step2$mean_r2, 4), "\n")

# Spatial: lat/lon and/or region (region included when present in data)
spatial_configs <- list(
  "+ lat/lon" = c("longitude", "latitude", env_vars)
)
if (has_region) {
  spatial_configs[["+ region"]] <- c(env_vars, "region")
  spatial_configs[["+ lat/lon + region"]] <- c("longitude", "latitude", env_vars, "region")
}

step3_results <- list("Env only" = best_step2)
for (config_name in names(spatial_configs)) {
  preds <- spatial_configs[[config_name]]
  out <- run_gpr_cv_parameter_search(preds, best_kernel, nug.min = best_nug_min, nug.est = TRUE)
  step3_results[[config_name]] <- out
  results_list[[length(results_list) + 1]] <- data.frame(
    step = "3. Spatial",
    configuration = config_name,
    mean_r2 = out$mean_r2,
    sd_r2 = out$sd_r2,
    mean_rmse = out$mean_rmse,
    sd_rmse = out$sd_rmse,
    n_folds_valid = out$n_folds_valid,
    stringsAsFactors = FALSE
  )
  cat("  ", config_name, ": R² =", round(out$mean_r2, 4), "\n")
}
best_spatial_name <- names(step3_results)[which.max(sapply(step3_results, function(x) x$mean_r2))]
best_spatial_predictors <- if (best_spatial_name == "Env only") env_vars else spatial_configs[[best_spatial_name]]
cat("  Best spatial:", best_spatial_name, "\n\n")

# -----------------------------------------------------------------------------
# Step 4: Categorical (carry forward step 3 best, then + species)
# -----------------------------------------------------------------------------
cat("Step 4: Categorical (kernel = ", best_kernel, ", spatial = ", best_spatial_name, ")...\n", sep = "")
best_step3 <- step3_results[[best_spatial_name]]
results_list[[length(results_list) + 1]] <- data.frame(
  step = "4. Categorical",
  configuration = "Baseline (spatial)",
  mean_r2 = best_step3$mean_r2,
  sd_r2 = best_step3$sd_r2,
  mean_rmse = best_step3$mean_rmse,
  sd_rmse = best_step3$sd_rmse,
  n_folds_valid = best_step3$n_folds_valid,
  stringsAsFactors = FALSE
)
cat("  Baseline (spatial) (from step 3): R² =", round(best_step3$mean_r2, 4), "\n")

step4_results <- list("Baseline (spatial)" = best_step3)
if (has_species) {
  out <- run_gpr_cv_parameter_search(c(best_spatial_predictors, "seagrass_species"), best_kernel, nug.min = best_nug_min, nug.est = TRUE)
  step4_results[["+ species"]] <- out
  results_list[[length(results_list) + 1]] <- data.frame(
    step = "4. Categorical",
    configuration = "+ species",
    mean_r2 = out$mean_r2,
    sd_r2 = out$sd_r2,
    mean_rmse = out$mean_rmse,
    sd_rmse = out$sd_rmse,
    n_folds_valid = out$n_folds_valid,
    stringsAsFactors = FALSE
  )
  cat("  + species: R² =", round(out$mean_r2, 4), "\n")
}
best_cat_name <- names(step4_results)[which.max(sapply(step4_results, function(x) x$mean_r2))]
best_final_predictors <- if (best_cat_name == "Baseline (spatial)") best_spatial_predictors else c(best_spatial_predictors, "seagrass_species")
cat("  Best categorical:", best_cat_name, "\n\n")

# -----------------------------------------------------------------------------
# Build table and save
# -----------------------------------------------------------------------------
summary_tab <- do.call(rbind, results_list)
summary_tab <- summary_tab %>%
  mutate(step = factor(step, levels = c("1. Nugget", "2. Kernel", "3. Spatial", "4. Categorical")))

write.csv(summary_tab, file.path(out_dir, "gpr_parameter_search_summary.csv"), row.names = FALSE)
cat("Saved:", file.path(out_dir, "gpr_parameter_search_summary.csv"), "\n")

# -----------------------------------------------------------------------------
# Single bar chart: step 1 at top, step 4 at bottom; within each step, lowest R² to highest
# Use unique plot_label per row (step + configuration) so repeated names don't overlap.
# -----------------------------------------------------------------------------
plot_tab <- summary_tab %>%
  mutate(plot_label = paste(step, configuration, sep = " — "))
# Order: step 4 at bottom → step 1 at top; within each step descending mean_r2 so that
# top-most bar = lowest R², bottom = highest R² (R² increases reading downwards)
plot_label_order <- plot_tab %>%
  arrange(desc(step), desc(mean_r2)) %>%
  pull(plot_label)
plot_tab <- plot_tab %>%
  mutate(plot_label = factor(plot_label, levels = plot_label_order))

step_colours <- c(
  "1. Nugget"      = "#e41a1c",
  "2. Kernel"      = "#4daf4a",
  "3. Spatial"     = "#377eb8",
  "4. Categorical" = "#984ea3"
)

# Display only configuration on axis (step is in the legend)
strip_step_prefix <- function(x) {
  vapply(strsplit(as.character(x), " — ", fixed = TRUE), function(z) z[2L], character(1L))
}

p <- ggplot(plot_tab, aes(x = plot_label, y = mean_r2, fill = step)) +
  geom_col(alpha = 0.9) +
  geom_errorbar(aes(ymin = mean_r2 - sd_r2, ymax = mean_r2 + sd_r2),
                width = 0.2, linewidth = 0.4, color = "gray20") +
  coord_flip() +
  scale_x_discrete(labels = strip_step_prefix) +
  scale_fill_manual(values = step_colours, name = "Step") +
  labs(
    x = "",
    y = "Mean R² (1 km spatial CV)",
    title = "GPR parameter search (sequential: 1. Nugget → 2. Kernel → 3. Spatial → 4. Categorical)",
    subtitle = "Regularization first, then model expansion. Adding spatial/species can lower R² in strict spatial CV (overfitting, redundancy)."
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "gray40"),
    legend.position = "right",
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = rel(0.9))
  )

ggsave(file.path(out_dir, "gpr_parameter_search_summary.png"), p, width = 9, height = 7, dpi = 150)
ggsave(file.path(out_dir, "gpr_parameter_search_summary.pdf"), p, width = 9, height = 7)
if (fig_dir != out_dir) {
  ggsave(file.path(fig_dir, "gpr_parameter_search_summary.png"), p, width = 9, height = 7, dpi = 150)
  ggsave(file.path(fig_dir, "gpr_parameter_search_summary.pdf"), p, width = 9, height = 7)
}
cat("Figure saved: gpr_parameter_search_summary.png (and .pdf)\n")

# Print summary
cat("\n========================================\n")
cat("SEQUENTIAL PARAMETER SEARCH TABLE\n")
cat("========================================\n\n")
print(summary_tab)
cat("\nBest: nug.min = ", format(best_nug_min, scientific = TRUE), ", kernel = ", best_kernel, ", spatial = ", best_spatial_name, ", categorical = ", best_cat_name, "\n", sep = "")
cat("Done.\n")
