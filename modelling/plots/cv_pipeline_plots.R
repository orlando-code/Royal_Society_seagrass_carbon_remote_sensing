# CV Pipeline: Plotting Only
# Run this script after cv_pipeline.R to regenerate all plots without re-running CV.
# Usage: source("modelling/cv_pipeline_plots.R") or Rscript modelling/cv_pipeline_plots.R

if (requireNamespace("here", quietly = TRUE)) setwd(here::here())

if (!requireNamespace("here", quietly = TRUE)) install.packages("here", quiet = TRUE)
library(here)
library(ggplot2)
library(dplyr)
library(patchwork)

# Load plot configuration (including model_colours)
# This file should define `model_colours`
source("modelling/R/plot_config.R")

# Prefer figures/cv_pipeline_output; fallback to modelling/cv_pipeline_output
out_dir <- "figures/cv_pipeline_output"
if (!dir.exists(out_dir)) out_dir <- "modelling/cv_pipeline_output"
if (exists("dpi", envir = .GlobalEnv)) dpi <- get("dpi", envir = .GlobalEnv) else dpi <- 300
if (exists("show_titles", envir = .GlobalEnv)) show_titles <- get("show_titles", envir = .GlobalEnv) else show_titles <- TRUE

required_files <- c(
  "cv_results_detailed.csv",
  "cv_results_summary_by_predictor_set.csv"
)
for (f in required_files) {
  if (!file.exists(file.path(out_dir, f))) {
    stop("Required file not found: ", file.path(out_dir, f), "\nRun cv_pipeline.R first to generate CV results.")
  }
}

cat("\n========================================\n")
cat("LOADING CV RESULTS FOR PLOTTING\n")
cat("========================================\n\n")

# Load results
results_df <- read.csv(file.path(out_dir, "cv_results_detailed.csv"), stringsAsFactors = FALSE)
summary_by_method_model <- read.csv(file.path(out_dir, "cv_results_summary_by_predictor_set.csv"), stringsAsFactors = FALSE)

# Load predictions if available (optional)
predictions_df <- NULL
if (file.exists(file.path(out_dir, "cv_predictions.csv"))) {
  predictions_df <- read.csv(file.path(out_dir, "cv_predictions.csv"), stringsAsFactors = FALSE)
  cat("Loaded predictions from cv_predictions.csv\n")
} else {
  cat("No cv_predictions.csv found; skipping plots that require predictions.\n")
}

# Convert internal method name to display label (for plots; handles e.g. spatial_block_1e+05m)
method_to_display <- function(m) {
  if (grepl("^spatial_block_.+m$", m)) {
    size_str <- sub("^spatial_block_(.+)m$", "\\1", m)
    size_m <- as.numeric(size_str)
    size_km <- if (!is.na(size_m) && size_m >= 1000) round(size_m / 1000) else size_m
    if (!is.na(size_km)) return(paste("Spatial block", size_km, "km"))
  }
  if (grepl("^distance_buffer_", m)) {
    km <- sub("distance_buffer_|km", "", m)
    return(paste("Distance buffer", km, "km"))
  }
  switch(m,
    random_split = "Random split",
    env_cluster = "Env. cluster",
    paste0(toupper(substring(m, 1, 1)), substring(m, 2))
  )
}

# Extract block size in km from method name (for ordering: 0 = random_split, NA = other, numeric = block size)
method_block_size_km <- function(m) {
  if (m == "random_split") return(0)
  if (grepl("^spatial_block_.+m$", m)) {
    size_str <- sub("^spatial_block_(.+)m$", "\\1", m)
    size_m <- as.numeric(size_str)
    if (!is.na(size_m)) return(round(size_m / 1000))
  }
  NA_real_
}

# Order methods: random_split first, then spatial blocks by size (1km, 5km, ..., 100km), then others
order_methods_by_size <- function(methods) {
  method_sizes <- vapply(methods, method_block_size_km, numeric(1))
  # Order: random_split (0), then spatial blocks (sorted by size), then others (NA)
  order_idx <- order(
    ifelse(is.na(method_sizes), 999, method_sizes),  # Put NA (non-spatial) last
    na.last = TRUE
  )
  methods[order_idx]
}

# Use only model colours for models that appear in the data (keep order for legend)
models_present <- unique(results_df$model)
model_colours <- model_colours[names(model_colours) %in% models_present]

cat("Generating plots...\n\n")

# Summary plot: R² and RMSE by method and model (only these two metrics)
# Crop RMSE y-axis to exclude extreme outliers that obscure the main data
metric_labeller <- c(r2 = "R²", rmse = "RMSE")
rmse_vals <- results_df$rmse[!is.na(results_df$rmse)]
rmse_upper <- if (length(rmse_vals) > 0) {
  q95 <- quantile(rmse_vals, 0.95, na.rm = TRUE)
  max(q95, 0.1, na.rm = TRUE)  # At least show up to 0.1
} else {
  0.1
}
p_metrics_data <- results_df %>%
  dplyr::select(method, model, r2, rmse) %>%
  mutate(
    method_label = vapply(method, method_to_display, character(1)),
    rmse_capped = pmin(rmse, rmse_upper)  # Cap extreme RMSE values for plotting
  )
# Order methods: random_split first, then spatial blocks by ascending size (1km to 100km)
method_order <- order_methods_by_size(unique(p_metrics_data$method))
method_label_order <- vapply(method_order, method_to_display, character(1))
p_metrics_data$method_label <- factor(p_metrics_data$method_label, levels = method_label_order)
# Order models by performance (highest mean R² first = leftmost in each dodged group)
model_order_r2 <- p_metrics_data %>%
  group_by(model) %>%
  summarise(m = mean(r2, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(m)) %>%
  pull(model)
p_metrics_data$model <- factor(p_metrics_data$model, levels = model_order_r2)
p_metrics <- p_metrics_data %>%
  tidyr::pivot_longer(cols = c(r2, rmse_capped), names_to = "metric", values_to = "value") %>%
  mutate(
    metric = ifelse(metric == "rmse_capped", "rmse", metric),
    # Order metrics: R² first (top), then RMSE (bottom)
    metric = factor(metric, levels = c("r2", "rmse"))
  ) %>%
  ggplot(aes(x = method_label, y = value, fill = model)) +
  geom_boxplot(outlier.size = 0.8, position = position_dodge(0.8), width = 0.7) +
  facet_wrap(~metric, scales = "free_y", ncol = 1, labeller = as_labeller(metric_labeller)) +
  scale_fill_manual(values = model_colours, na.value = "gray75") +
  labs(x = "CV method", y = NULL, title = if (show_titles) "Cross-validation performance by method and model" else NULL) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5)
  )
ggsave(file.path(out_dir, "cv_metrics_by_method_model.png"), p_metrics, width = 12, height = 8, dpi = dpi)
cat("Saved", file.path(out_dir, "cv_metrics_by_method_model.png"), "\n")
if (any(results_df$rmse > rmse_upper, na.rm = TRUE)) {
  n_capped <- sum(results_df$rmse > rmse_upper, na.rm = TRUE)
  cat("(RMSE values >", round(rmse_upper, 4), "were capped for display;", n_capped, "values affected)\n")
}

# Box plot using each model's best predictor set (best R² at 1 km spatial CV) across all methods
best_1km_for_plot <- NULL
if ("predictor_set" %in% names(results_df) && "block_size_m" %in% names(results_df)) {
  spatial_1km <- results_df %>%
    filter(grepl("^spatial_block_", method), block_size_m == 1000)
  if (nrow(spatial_1km) > 0) {
    best_1km_for_plot <- spatial_1km %>%
      group_by(model, predictor_set) %>%
      summarise(mean_r2_1km = mean(r2, na.rm = TRUE), .groups = "drop") %>%
      group_by(model) %>%
      slice_max(mean_r2_1km, n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      select(model, predictor_set)
  }
}
if (!is.null(best_1km_for_plot) && nrow(best_1km_for_plot) > 0) {
  results_best <- results_df %>%
    inner_join(best_1km_for_plot, by = c("model", "predictor_set"))
  method_order_best <- order_methods_by_size(unique(results_best$method))
  method_label_order_best <- vapply(method_order_best, method_to_display, character(1))
  model_order_best <- results_best %>%
    group_by(model) %>%
    summarise(m = mean(r2, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(m)) %>%
    pull(model)
  p_metrics_best <- results_best %>%
    mutate(
      method_label = vapply(method, method_to_display, character(1)),
      method_label = factor(method_label, levels = method_label_order_best),
      model = factor(model, levels = model_order_best),
      rmse_capped = pmin(rmse, rmse_upper)
    ) %>%
    tidyr::pivot_longer(cols = c(r2, rmse_capped), names_to = "metric", values_to = "value") %>%
    mutate(metric = ifelse(metric == "rmse_capped", "rmse", metric), metric = factor(metric, levels = c("r2", "rmse"))) %>%
    ggplot(aes(x = method_label, y = value, fill = model)) +
    geom_boxplot(outlier.size = 0.8, position = position_dodge(0.8), width = 0.7) +
    facet_wrap(~ metric, scales = "free_y", ncol = 1, labeller = as_labeller(metric_labeller)) +
    scale_fill_manual(values = model_colours, na.value = "gray75") +
    labs(x = "CV method", y = NULL, title = if (show_titles) "Cross-validation performance by method and model (best predictor set per model)" else NULL) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), legend.position = "bottom", plot.title = element_text(hjust = 0.5))
  ggsave(file.path(out_dir, "cv_metrics_by_method_model_best_predictors.png"), p_metrics_best, width = 12, height = 8, dpi = dpi)
  cat("Saved", file.path(out_dir, "cv_metrics_by_method_model_best_predictors.png"), "\n")
}

# Scatter: observed vs predicted for all methods (one panel per model)
if (!is.null(predictions_df) && nrow(predictions_df) > 0) {
  pred_all <- predictions_df %>% filter(if ("predictor_set" %in% names(predictions_df)) predictor_set == "all" else TRUE)
  if (nrow(pred_all) > 0) {
    pred_all <- pred_all %>% mutate(method_label = vapply(method, method_to_display, character(1)))
    p_scatter <- ggplot(pred_all, aes(x = predicted, y = observed)) +
      geom_point(alpha = 0.4, size = 1) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "red", linewidth = 0.5) +
      facet_wrap(~model, scales = "free", ncol = 3) +
      labs(x = "Predicted", y = "Observed", title = if (show_titles) "Observed vs predicted by model (all CV folds)" else NULL) +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5))
    ggsave(file.path(out_dir, "cv_observed_vs_predicted.png"), p_scatter, width = 12, height = 10, dpi = dpi)
    cat("Saved", file.path(out_dir, "cv_observed_vs_predicted.png"), "\n")
  }
}


# Direct comparison: GAM vs ML methods under the same CV (mean R² by model and CV method)
# Filter to full predictor set for clean comparison
results_all <- results_df %>%
  filter(if ("predictor_set" %in% names(results_df)) predictor_set == "all" else TRUE)
summary_plot <- results_all %>%
  group_by(method, model) %>%
  summarise(mean_r2 = mean(r2, na.rm = TRUE), sd_r2 = sd(r2, na.rm = TRUE), .groups = "drop") %>%
  mutate(method_label = vapply(method, method_to_display, character(1)))
# Order methods: random_split first, then spatial blocks by size
method_order_plot <- order_methods_by_size(unique(summary_plot$method))
method_label_order_plot <- vapply(method_order_plot, method_to_display, character(1))
summary_plot$method_label <- factor(summary_plot$method_label, levels = method_label_order_plot)
# Order models by performance (highest R² at left; with coord_flip, first level is at bottom = left)
p_gam_vs_ml <- ggplot(summary_plot, aes(x = reorder(model, -mean_r2), y = mean_r2, fill = model)) +
  geom_col() +
  geom_errorbar(aes(ymin = pmax(0, mean_r2 - sd_r2), ymax = pmin(1, mean_r2 + sd_r2)), width = 0.2, linewidth = 0.4) +
  facet_wrap(~method_label, scales = "free_x", ncol = 2) +
  scale_fill_manual(values = model_colours, na.value = "gray75") +
  coord_flip() +
  labs(x = NULL, y = "Mean R² (± SD across folds)", title = if (show_titles) "GAM vs ML methods: same CV test" else NULL) +
  theme_minimal() +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
ggsave(file.path(out_dir, "cv_gam_vs_ml_by_method.png"), p_gam_vs_ml, width = 10, height = 4 + 1.2 * length(unique(summary_plot$method)), dpi = dpi, limitsize = FALSE)
cat("Saved", file.path(out_dir, "cv_gam_vs_ml_by_method.png"), "\n")

# Two-panel: Random split vs Spatial block (first available) – direct side-by-side
spatial_methods <- unique(summary_plot$method[grepl("^spatial_block_", summary_plot$method)])
# Order spatial methods by size and take the first (smallest)
if (length(spatial_methods) > 0) {
  spatial_methods_ordered <- order_methods_by_size(spatial_methods)
  spatial_methods_ordered <- spatial_methods_ordered[spatial_methods_ordered != "random_split"]
  first_spatial <- if (length(spatial_methods_ordered) > 0) spatial_methods_ordered[1] else NULL
} else {
  first_spatial <- NULL
}
compare_methods <- c("random_split", first_spatial)
if (length(compare_methods) >= 2 && !is.null(first_spatial)) {
  summary_compare <- summary_plot %>% filter(method %in% compare_methods)
  # Order: random_split first, then spatial block
  compare_methods_ordered <- order_methods_by_size(compare_methods)
  compare_labels_ordered <- vapply(compare_methods_ordered, method_to_display, character(1))
  summary_compare$method_label <- factor(
    vapply(summary_compare$method, method_to_display, character(1)),
    levels = compare_labels_ordered
  )
  # Order models by mean R² (highest at left; with coord_flip, first level at bottom = left)
  model_order_compare <- summary_compare %>%
    group_by(model) %>%
    summarise(m = mean(mean_r2, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(m)) %>%
    pull(model)
  summary_compare$model <- factor(summary_compare$model, levels = model_order_compare)
  cv_method_colours <- setNames(c("#377eb8", "#e41a1c"), levels(summary_compare$method_label))
  p_random_vs_spatial <- ggplot(summary_compare, aes(x = model, y = mean_r2, fill = method_label)) +
    geom_col(position = position_dodge(0.9), width = 0.8) +
    geom_errorbar(aes(ymin = pmax(0, mean_r2 - sd_r2), ymax = pmin(1, mean_r2 + sd_r2)),
                  position = position_dodge(0.9), width = 0.2, linewidth = 0.4) +
    scale_fill_manual(values = cv_method_colours, name = NULL) +
    coord_flip() +
    labs(x = NULL, y = "Mean R² (± SD)", title = if (show_titles) "Random split vs spatial block CV: GAM and ML methods" else NULL) +
    theme_minimal() +
    theme(legend.title = element_blank(), legend.position = "bottom", plot.title = element_text(hjust = 0.5))
  ggsave(file.path(out_dir, "cv_random_vs_spatial_gam_ml.png"), p_random_vs_spatial, width = 9, height = 5, dpi = dpi)
  cat("Saved", file.path(out_dir, "cv_random_vs_spatial_gam_ml.png"), "\n")
}

# Bar summary: mean R² and mean RMSE by method (averaged across models)
grp_fold <- c("method", "fold")
if ("predictor_set" %in% names(results_df)) grp_fold <- c(grp_fold, "predictor_set")
grp_method <- c("method")
if ("predictor_set" %in% names(results_df)) grp_method <- c(grp_method, "predictor_set")

summary_method <- results_df %>%
  group_by(across(all_of(grp_fold))) %>%
  summarise(
    mean_r2 = mean(r2, na.rm = TRUE),
    mean_rmse = mean(rmse, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(across(all_of(grp_method))) %>%
  summarise(
    r2 = mean(mean_r2, na.rm = TRUE),
    sd_r2 = sd(mean_r2, na.rm = TRUE),
    rmse = mean(mean_rmse, na.rm = TRUE),
    sd_rmse = sd(mean_rmse, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(method_label = vapply(method, method_to_display, character(1)))
# Order methods: random_split first, then spatial blocks by size
method_order_method <- order_methods_by_size(unique(summary_method$method))
method_label_order_method <- vapply(method_order_method, method_to_display, character(1))
summary_method$method_label <- factor(summary_method$method_label, levels = method_label_order_method)
p_method <- summary_method %>%
  tidyr::pivot_longer(cols = c(r2, rmse), names_to = "metric", values_to = "value") %>%
  ggplot(aes(x = method_label, y = value, fill = metric)) +
  geom_col(position = "dodge") +
  labs(x = "CV method", y = "Mean (across folds and models)", title = if (show_titles) "Ensemble performance by CV method" else NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
if ("predictor_set" %in% names(summary_method)) {
  p_method <- p_method + facet_wrap(~predictor_set, ncol = 1)
}
ggsave(file.path(out_dir, "cv_ensemble_by_method.png"), p_method, width = 10, height = 5, dpi = dpi)
cat("Saved", file.path(out_dir, "cv_ensemble_by_method.png"), "\n")

# Performance over block sizes (1, 5, 25, 50, 100 km)
# Use, for each model, the predictor set that performs best at 1 km spatial block CV,
# then track that same predictor set across all block sizes.
spatial_results <- results_df %>%
  filter(grepl("^spatial_block_", method), is.na(block_size_m) == FALSE)
if (nrow(spatial_results) > 0 && "block_size_m" %in% names(spatial_results)) {
  best_1km <- NULL
  if ("predictor_set" %in% names(spatial_results)) {
    best_1km <- spatial_results %>%
      filter(block_size_m == 1000) %>%
      group_by(model, predictor_set) %>%
      summarise(mean_r2_1km = mean(r2, na.rm = TRUE),
                mean_rmse_1km = mean(rmse, na.rm = TRUE),
                .groups = "drop") %>%
      group_by(model) %>%
      slice_max(mean_r2_1km, n = 1, with_ties = FALSE) %>%
      ungroup()

    if (nrow(best_1km) > 0) {
      write.csv(
        best_1km,
        file.path(out_dir, "cv_best_predictor_set_by_model_1km.csv"),
        row.names = FALSE
      )
      spatial_results <- spatial_results %>%
        dplyr::inner_join(best_1km %>% dplyr::select(model, predictor_set),
                          by = c("model", "predictor_set"))
    }
  }

  spatial_best <- spatial_results %>%
    group_by(block_size_m, model) %>%
    summarise(mean_r2 = mean(r2, na.rm = TRUE),
              mean_rmse = mean(rmse, na.rm = TRUE),
              .groups = "drop")

  sum_block <- spatial_best %>%
    mutate(block_size_km = round(block_size_m / 1000)) %>%
    group_by(block_size_km, model) %>%
    summarise(mean_r2 = mean(mean_r2, na.rm = TRUE),
              mean_rmse = mean(mean_rmse, na.rm = TRUE),
              .groups = "drop")
  km_levels <- sort(unique(sum_block$block_size_km))
  sum_block <- sum_block %>%
    mutate(block_km = factor(block_size_km, levels = km_levels))
  # Order models by mean R² (highest first) for consistent legend and line order
  model_order_block <- sum_block %>%
    group_by(model) %>%
    summarise(m = mean(mean_r2, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(m)) %>%
    pull(model)
  sum_block$model <- factor(sum_block$model, levels = model_order_block)
  p_block_r2 <- ggplot(sum_block, aes(x = block_km, y = mean_r2, colour = model, group = model)) +
    geom_line(linewidth = 1) + geom_point(size = 3) +
    scale_colour_manual(values = model_colours, na.value = "gray75") +
    labs(x = "Block size (km)", y = "Mean R²", title = if (show_titles) "Performance over spatial block sizes" else NULL) +
    scale_x_discrete(labels = paste0(km_levels, " km")) +
    theme_minimal() + theme(legend.position = "bottom", plot.title = element_text(size = 10, hjust = 0.5))
  # Cap RMSE y-axis so main models are visible (high-RMSE models e.g. Universal Kriging may be off-scale)
  rmse_vals_block <- sum_block$mean_rmse[is.finite(sum_block$mean_rmse)]
  rmse_cap_block <- if (length(rmse_vals_block) > 0) {
    q <- quantile(rmse_vals_block, 0.92, na.rm = TRUE)
    max(0.25, min(q, 500))
  } else 0.25
  p_block_rmse <- ggplot(sum_block, aes(x = block_km, y = mean_rmse, colour = model, group = model)) +
    geom_line(linewidth = 1) + geom_point(size = 3) +
    scale_colour_manual(values = model_colours, na.value = "gray75") +
    coord_cartesian(ylim = c(0, rmse_cap_block)) +
    labs(x = "Block size (km)", y = "Mean RMSE", title = if (show_titles) "Performance over spatial block sizes" else NULL,
         subtitle = if (show_titles) paste0("RMSE axis capped at ", round(rmse_cap_block, 2), "; some models may be off-scale.") else NULL) +
    scale_x_discrete(labels = paste0(km_levels, " km")) +
    theme_minimal() + theme(legend.position = "bottom", plot.title = element_text(size = 10, hjust = 0.5), plot.subtitle = element_text(size = 7))
  p_block_combined <- p_block_r2 + p_block_rmse + plot_layout(ncol = 1)
  ggsave(file.path(out_dir, "cv_performance_over_block_sizes.png"), p_block_combined, width = 9, height = 7, dpi = dpi)
  cat("Saved", file.path(out_dir, "cv_performance_over_block_sizes.png"), "\n")
}

# Model selection for spatial prediction with uncertainty: highlight models that provide prediction intervals
key_methods <- c("random_split", "spatial_block_1000m", "spatial_block_5000m")
summary_unc <- summary_by_method_model %>%
  filter(method %in% key_methods)
if ("predictor_set" %in% names(summary_unc)) summary_unc <- summary_unc %>% filter(predictor_set == "all")
if (nrow(summary_unc) > 0) {
  summary_unc <- summary_unc %>%
    mutate(method_label = vapply(method, method_to_display, character(1)))
  # Order methods: random_split first, then spatial blocks by size
  method_order_unc <- order_methods_by_size(key_methods)
  method_order_labels <- vapply(method_order_unc, method_to_display, character(1))
  summary_unc$method_label <- factor(summary_unc$method_label, levels = method_order_labels)
  # Order models by performance (highest mean R² first = leftmost in each dodged group)
  model_order_unc <- summary_unc %>%
    group_by(model) %>%
    summarise(m = mean(mean_r2, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(m)) %>%
    pull(model)
  summary_unc$model <- factor(summary_unc$model, levels = model_order_unc)
  p_unc_simple <- ggplot(summary_unc, aes(x = method_label, y = mean_r2, fill = model)) +
    geom_col(position = position_dodge(0.85), width = 0.75) +
    geom_errorbar(aes(ymin = pmax(0, mean_r2 - sd_r2), ymax = pmin(1, mean_r2 + sd_r2)),
                  position = position_dodge(0.85), width = 0.2, linewidth = 0.4) +
    scale_fill_manual(values = model_colours, na.value = "gray75", name = "Model") +
    labs(
      x = "CV method",
      y = "Mean R² (± SD)",
      title = if (show_titles) "Model selection for spatial prediction with uncertainty" else NULL,
      subtitle = if (show_titles) "GAM and GPR provide prediction intervals; others give point predictions only. Spatial block CV tests transfer to nearby, unsampled areas." else NULL
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "bottom", plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(size = 8))
  ggsave(file.path(out_dir, "cv_uncertainty_capable_models.png"), p_unc_simple, width = 8, height = 5, dpi = dpi)
  cat("Saved", file.path(out_dir, "cv_uncertainty_capable_models.png"), "\n")
}

# All vs pruned predictors: compare R² by model and predictor set for key CV methods
predictor_set_labels <- c(
  all = "All predictors",
  pruned = "Pruned (GAM)",
  pruned_xgboost = "Pruned (XGBoost)",
  pruned_gpr = "Pruned (GPR, model permutation)"
)
if ("predictor_set" %in% names(summary_by_method_model) && length(unique(summary_by_method_model$predictor_set)) >= 2) {
  key_methods_ps <- c("random_split", "spatial_block_1000m", "spatial_block_5000m")
  summary_ps <- summary_by_method_model %>% filter(method %in% key_methods_ps)
  if (nrow(summary_ps) > 0) {
    summary_ps <- summary_ps %>%
      mutate(
        method_label = vapply(method, method_to_display, character(1)),
        predictor_set_label = ifelse(predictor_set %in% names(predictor_set_labels),
          predictor_set_labels[predictor_set], paste0("Pruned (", predictor_set, ")")
        )
      )
    # Order methods: random_split first, then spatial blocks by size
    method_order_ps <- order_methods_by_size(key_methods_ps)
    method_order_ps_labels <- vapply(method_order_ps, method_to_display, character(1))
    summary_ps$method_label <- factor(summary_ps$method_label, levels = method_order_ps_labels)
    # Order models by performance (highest mean R² first = leftmost in each dodged group)
    model_order_ps <- summary_ps %>%
      group_by(model) %>%
      summarise(m = mean(mean_r2, na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(m)) %>%
      pull(model)
    summary_ps$model <- factor(summary_ps$model, levels = model_order_ps)
    n_ps_facets <- length(unique(summary_ps$predictor_set_label))
    p_all_vs_pruned <- ggplot(summary_ps, aes(x = method_label, y = mean_r2, fill = model)) +
      geom_col(position = position_dodge(0.85), width = 0.75) +
      geom_errorbar(aes(ymin = pmax(0, mean_r2 - sd_r2), ymax = pmin(1, mean_r2 + sd_r2)),
                    position = position_dodge(0.85), width = 0.2, linewidth = 0.4) +
      scale_fill_manual(values = model_colours, na.value = "gray75", name = "Model") +
      facet_wrap(~ predictor_set_label, ncol = min(3, n_ps_facets)) +
      labs(
        x = "CV method",
        y = "Mean R² (± SD)",
        title = if (show_titles) "All vs pruned predictor sets (GAM, XGBoost, GPR)" else NULL,
        subtitle = if (show_titles) "Same folds; pruned sets reduce overfitting risk." else NULL
      ) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "bottom", plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(size = 8))
    ggsave(file.path(out_dir, "cv_all_vs_pruned_predictors.png"), p_all_vs_pruned, width = 4 + 3 * n_ps_facets, height = 5.5, dpi = dpi)
    cat("Saved", file.path(out_dir, "cv_all_vs_pruned_predictors.png"), "\n")
  }
}

# Note: cv_spatial_blocks_and_points.png requires sf objects (core_sf, cv_strategies) that are not saved.
# That plot can only be generated when running the full cv_pipeline.R script.

cat("\n========================================\n")
cat("All plots generated successfully!\n")
cat("========================================\n\n")
cat("Note: cv_spatial_blocks_and_points.png requires spatial data structures and is only generated\n")
cat("      when running the full cv_pipeline.R script.\n\n")
