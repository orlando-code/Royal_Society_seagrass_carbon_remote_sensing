# Predict seagrass carbon density at points and compare with measured values
#
# First install the renv package and run `renv::restore()` to install the
# necessary packages.

sys.source("modelling/R/init_repo.R", envir = .GlobalEnv)
project_root <- seagrass_init_repo(
  packages = c("dplyr", "readxl", "readr", "ggplot2"),
  source_files = c(
    "modelling/R/helpers.R",
    "modelling/R/ml.R",
    "modelling/R/extract_covariates_from_rasters.R"
  ),
  include_helpers = FALSE,
  require_core_inputs = FALSE,
  check_renv = FALSE
)
sys.source("review/env_training_comparison.R", envir = .GlobalEnv)

# -----------------------------------------------------------------------------
# Config
# -----------------------------------------------------------------------------
dat_fp <- file.path("data", "review", "PercOC_only_unique_reduced_df.xlsx")
model_fp <- file.path("data", "review", "GPR_final.rds")
output_fp <- file.path("output", "review", "PercOC_only_unique_reduced_df_predictions.csv")
plot_output_dir <- file.path("output", "review", "additional_samples_comparison")
env_comparison_output_dir <- file.path(plot_output_dir, "env_training_comparison")
train_data_fp <- file.path(project_root, "data", "all_extracted_new.rds")
measured_carbon_density_column <- "median_carbon_density_100cm_calc"

# make sure output directories exist
dir.create(dirname(output_fp), showWarnings = FALSE, recursive = TRUE)
dir.create(plot_output_dir, showWarnings = FALSE, recursive = TRUE)
# -----------------------------------------------------------------------------
# Load sample coordinates from file 
# -----------------------------------------------------------------------------
if (grepl("\\.xlsx$", dat_fp, ignore.case = TRUE)) {
  df_raw <- readxl::read_excel(dat_fp, sheet = 1)
} else if (grepl("\\.csv$", dat_fp, ignore.case = TRUE)) {
  df_raw <- readr::read_csv(dat_fp, show_col_types = FALSE)
} else {
  stop("File extension not recognized. Please provide a .csv or .xlsx file.")
}
required_columns <- c("longitude", "latitude", "seagrass_species")
missing_cols <- setdiff(required_columns, names(df_raw))
if (length(missing_cols) > 0) {
  stop(
    "Input data is missing required columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

points <- df_raw %>%
  dplyr::select(dplyr::all_of(required_columns))

# -----------------------------------------------------------------------------
# Load model and extract environmental covariates at point locations
# -----------------------------------------------------------------------------
model <- readRDS(model_fp)
predictor_vars <- model$predictor_vars

raster_covars <- setdiff(
  predictor_vars,
  c("seagrass_species", "longitude", "latitude")
)
missing_rasters <- setdiff(tolower(raster_covars), raster_covariates)
if (length(missing_rasters) > 0) {
  stop(
    "Model requires raster covariate(s) not found under data/env_rasters directory: ",
    paste(missing_rasters, collapse = ", ")
  )
}

cat("\n\nExtracting environmental covariates at", nrow(points), "point(s)...\n")
pred_data <- extract_covariates_at_points(
  points = points,
  covariates = raster_covars,
  use_closest = TRUE
)
pred_data <- process_rs_covariates(pred_data)

na_covars <- raster_covars[
  vapply(raster_covars, function(v) any(is.na(pred_data[[v]]), na.rm = TRUE), logical(1))
]
if (length(na_covars) > 0) {
  warning(
    "Missing covariate values at some points: ",
    paste(na_covars, collapse = ", "),
    call. = FALSE
  )
}

# -----------------------------------------------------------------------------
# Compare environmental vectors to training data
# -----------------------------------------------------------------------------
train_data <- readRDS(train_data_fp)
env_cmp <- compare_prediction_env_to_training(
  training_data = train_data,
  prediction_data = pred_data,
  predictor_vars = predictor_vars,
  output_dir = env_comparison_output_dir,
  make_plots = TRUE
)

# -----------------------------------------------------------------------------
# Predict carbon density
# -----------------------------------------------------------------------------
cat("\n\nPredicting carbon density at", nrow(points), "point(s)...\n")
pred <- predict_model(model, pred_data, se = TRUE)

# calculate carbon density
results <- dplyr::bind_cols(
  df_raw,
  tibble::tibble(
    predicted_carbon_density = pred$mean,
    predicted_se = pred$se
  ),
  env_cmp$point_summary %>%
    dplyr::select(
      env_similarity_score,
      env_flag_low_similarity,
      env_flag_any_predictor_outside_range,
      env_n_predictors_outside_range,
      env_pct_predictors_outside_range,
      dplyr::any_of("env_flag_novel_species"),
      env_flag_extrapolation
    )
)

# -----------------------------------------------------------------------------
# Compare predictions with measured values
# -----------------------------------------------------------------------------
cat("\n\nComparing predictions with species-mean predictions and measured values...\n")
if (!measured_carbon_density_column %in% names(results)) {
  stop("Measured column not found: ", measured_carbon_density_column)
}

results <- results %>%
  mutate(
    measured_carbon_density = .data[[measured_carbon_density_column]],
    density_residual = measured_carbon_density - predicted_carbon_density,
  )


# compute species-mean predictions (current state-of-the-art)
species_means <- compute_train_species_means(
  train_data,
  response_var = "median_carbon_density_100cm",
  species_var = "seagrass_species"
)
global_mean <- mean(train_data$median_carbon_density_100cm, na.rm = TRUE)
results <- results %>%
  mutate(
    species_mean_predicted_carbon_density = lookup_species_means(
      species_means,
      seagrass_species,
      fallback = global_mean
    ),
    species_mean_density_residual = measured_carbon_density - species_mean_predicted_carbon_density,
  )

model_density_metrics <- calculate_metrics(
  results$measured_carbon_density,
  results$predicted_carbon_density
)
species_density_metrics <- calculate_metrics(
  results$measured_carbon_density,
  results$species_mean_predicted_carbon_density
)


# -----------------------------------------------------------------------------
# Plot functions
# -----------------------------------------------------------------------------

comparison_colors <- c(
  "Model prediction" = "#2166AC",
  "Species-mean prediction" = "#D6604D"
)

save_and_show_plot <- function(plot, filename, width = 7, height = 5, dpi = 150) {
  out_path <- file.path(plot_output_dir, filename)
  ggplot2::ggsave(out_path, plot, width = width, height = height, dpi = dpi)
  if (interactive()) {
    print(plot)
  }
  cat("Saved plot:", out_path, "\n")
  invisible(out_path)
}

format_metrics_label <- function(metrics, label) {
  sprintf(
    "%s: R\u00b2 = %.3f, RMSE = %.4f",
    label,
    metrics$r2,
    metrics$rmse
  )
}

metrics_annotation <- function(model_metrics, species_metrics) {
  annotate(
    "text",
    x = Inf,
    y = Inf,
    hjust = 1.05,
    vjust = 1.4,
    label = paste(
      format_metrics_label(model_metrics, "Model"),
      format_metrics_label(species_metrics, "Species mean"),
      sep = "\n"
    ),
    size = 3.2,
    lineheight = 0.95,
    colour = "grey20"
  )
}

add_comparison_scatter <- function(plot, y_model, y_species, x_measured) {
  plot +
    geom_point(aes(x = .data[[x_measured]], y = .data[[y_model]], color = "Model prediction")) +
    geom_point(aes(x = .data[[x_measured]], y = .data[[y_species]], color = "Species-mean prediction")) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
    scale_color_manual(name = NULL, values = comparison_colors)
}

add_comparison_residuals <- function(plot, x_model, x_species, y_model, y_species) {
  plot +
    geom_point(aes(x = .data[[x_model]], y = .data[[y_model]], color = "Model prediction")) +
    geom_point(aes(x = .data[[x_species]], y = .data[[y_species]], color = "Species-mean prediction")) +
    geom_abline(slope = 0, intercept = 0, linetype = "dashed", color = "grey40") +
    scale_color_manual(name = NULL, values = comparison_colors)
}

# -----------------------------------------------------------------------------
# Plot results and compare model predictions with species-mean predictions
# -----------------------------------------------------------------------------
cat("\n\nPlotting and writing results...\n")
dir.create(plot_output_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Sanity check: model predictions on original training data – shows a nice scatter along 1:1 line
# -----------------------------------------------------------------------------
cat("\nSanity check: predicting on original training data (all_extracted_new.rds)...\n")
training_data <- process_rs_covariates(train_data)
training_data <- ensure_categorical_factors(training_data, model$predictor_vars)

# Exclude training rows with species absent from the fitted model.
training_species_levels <- if (identical(model$encoding$type, "gpr_species_adjusted") &&
                              !is.null(model$species_means)) {
  names(model$species_means)
} else {
  model$encoding$levels$seagrass_species %||% character(0)
}
if (length(training_species_levels) > 0L && "seagrass_species" %in% names(training_data)) {
  training_eval <- training_data %>%
    dplyr::filter(as.character(.data$seagrass_species) %in% training_species_levels)
  n_dropped <- nrow(training_data) - nrow(training_eval)
  if (n_dropped > 0L) {
    cat(
      "  Excluded ", n_dropped, " training row(s) with species not encoded in the saved model.\n",
      sep = ""
    )
  }
} else {
  training_eval <- training_data
}
training_pred <- predict_model(model, training_eval, se = FALSE)
stopifnot(length(training_pred$mean) == nrow(training_eval))

training_results <- training_eval %>%
  mutate(
    measured_carbon_density = median_carbon_density_100cm,
    predicted_carbon_density = training_pred$mean,
    species_mean_predicted_carbon_density = lookup_species_means(
      species_means,
      seagrass_species,
      fallback = global_mean
    )
  ) %>%
  filter(
    is.finite(measured_carbon_density),
    is.finite(predicted_carbon_density),
    is.finite(species_mean_predicted_carbon_density)
  )

training_model_metrics <- calculate_metrics(
  training_results$measured_carbon_density,
  training_results$predicted_carbon_density
)
training_species_metrics <- calculate_metrics(
  training_results$measured_carbon_density,
  training_results$species_mean_predicted_carbon_density
)

save_and_show_plot(
  add_comparison_scatter(
    ggplot(training_results, aes()),
    "predicted_carbon_density",
    "species_mean_predicted_carbon_density",
    "measured_carbon_density"
  ) +
    metrics_annotation(training_model_metrics, training_species_metrics) +
    labs(
      title = "Sanity check: original training data",
      x = "Measured carbon density",
      y = "Predicted carbon density"
    ) +
    theme_minimal() +
    theme(plot.margin = ggplot2::margin(5.5, 12, 5.5, 5.5)),
  "sanity_check_training_data_predicted_vs_measured.png",
  width = 7.5,
  height = 5.5
)

# -----------------------------------------------------------------------------
# Additional samples: compare model predictions with species-mean predictions (the current state of the art baseline)
# -----------------------------------------------------------------------------
save_and_show_plot(
  add_comparison_scatter(
    ggplot(results, aes()),
    "predicted_carbon_density",
    "species_mean_predicted_carbon_density",
    "measured_carbon_density"
  ) +
    metrics_annotation(model_density_metrics, species_density_metrics) +
    labs(x = "Measured carbon density", y = "Predicted carbon density") +
    theme_minimal() +
    theme(plot.margin = ggplot2::margin(5.5, 12, 5.5, 5.5)),
  "predicted_vs_measured_carbon_density.png",
  width = 7.5,
  height = 5.5
)


save_and_show_plot(
  add_comparison_residuals(
    ggplot(results, aes()),
    "predicted_carbon_density",
    "species_mean_predicted_carbon_density",
    "density_residual",
    "species_mean_density_residual"
  ) +
    metrics_annotation(model_density_metrics, species_density_metrics) +
    labs(x = "Predicted carbon density", y = "Density residual (measured - predicted)") +
    theme_minimal() +
    theme(plot.margin = ggplot2::margin(5.5, 12, 5.5, 5.5)),
  "carbon_density_residuals.png",
  width = 7.5,
  height = 5.5
)

# -----------------------------------------------------------------------------
# Carbon density distributions: training vs additional samples
# -----------------------------------------------------------------------------
distribution_colors <- c(
  "Training" = "#4D4D4D",
  "Additional samples" = "#2166AC"
)

additional_species <- sort(unique(as.character(results$seagrass_species)))

distribution_plot_data <- dplyr::bind_rows(
  train_data %>%
    dplyr::transmute(
      seagrass_species = as.character(.data$seagrass_species),
      carbon_density = .data$median_carbon_density_100cm,
      dataset = "Training"
    ),
  results %>%
    dplyr::transmute(
      seagrass_species = as.character(.data$seagrass_species),
      carbon_density = .data$measured_carbon_density,
      dataset = "Additional samples"
    )
) %>%
  dplyr::filter(
    .data$seagrass_species %in% additional_species,
    is.finite(.data$carbon_density)
  )

distribution_mean_lines <- distribution_plot_data %>%
  dplyr::group_by(.data$seagrass_species, .data$dataset) %>%
  dplyr::summarise(
    mean_carbon_density = mean(.data$carbon_density),
    n = dplyr::n(),
    .groups = "drop"
  )

save_and_show_plot(
  ggplot(distribution_plot_data, aes(x = .data$carbon_density, fill = .data$dataset, colour = .data$dataset)) +
    geom_histogram(
      bins = 30,
      boundary = 0,
      position = "identity",
      alpha = 0.45
    ) +
    geom_vline(
      data = distribution_mean_lines,
      aes(xintercept = .data$mean_carbon_density, colour = .data$dataset),
      linetype = "dashed",
      linewidth = 0.8
    ) +
    facet_wrap(~ .data$seagrass_species, scales = "free") +
    scale_fill_manual(name = NULL, values = distribution_colors) +
    scale_colour_manual(name = NULL, values = distribution_colors) +
    labs(
      title = "Carbon density: training vs additional samples",
      subtitle = "Dashed vertical lines show the species mean in each dataset",
      x = "Carbon density (g C/cm\u00b3)",
      y = "Count"
    ) +
    theme_minimal() +
    theme(
      plot.margin = ggplot2::margin(5.5, 12, 5.5, 5.5),
      strip.text = element_text(face = "bold")
    ),
  "carbon_density_training_vs_additional_histogram.png",
  width = 10,
  height = 6.5
)


# -----------------------------------------------------------------------------
# Write results to file
# -----------------------------------------------------------------------------
dir.create("output", showWarnings = FALSE, recursive = TRUE)
readr::write_csv(results, output_fp)
cat("\n\nWrote predictions to", output_fp, "\n")

