### Investigate effect of different cross-validation methods on model performance
# Refresh working directory
rm(list = ls())

# TODO: 
# move functions to separate file ,/
# hyperparameter tuning for models ,/
# investigate two-model approach: one model for aggregate, the second capturing depth effects
# compare effect of different target parameterisations (e.g. mean, sum, median) on model performance
# covariate pruning/feature importance investigation


# Set working directory
setwd("/Users/rt582/Library/CloudStorage/OneDrive-UniversityofCambridge/cambridge/phd/Paper_Conferences/seagrass/Fw_ Following up on modeling conversation at the Royal Society Meeting/SG_Cstock_modeling")

# source helper functions
source("helpers.R")

load_packages(c("here", "mgcv", "tidyverse", "ggplot2", "visreg", "randomForest", "corrplot", "blockCV"))

# Open dataset prepared for modeling (this is the same dataset used by Katalin Blix for machine learning models and
# myself for non-machine learning models)
For_modeling_df_shallow_carbon_density <- read_rds("data/For_modeling_df_shallow_carbon_density.rds")
summary(For_modeling_df_shallow_carbon_density)
dat <- For_modeling_df_shallow_carbon_density

# SPATIAL AUTOCORRELATION INVESTIGATION
# plot spatial distribution of unique data points
lon_min <- min(dat$longitude, na.rm = TRUE)
lon_max <- max(dat$longitude, na.rm = TRUE)
lat_min <- min(dat$latitude, na.rm = TRUE)
lat_max <- max(dat$latitude, na.rm = TRUE)
# background map
world <- map_data("world")
# get sample count for each unique location. Notice that some samples from exactly the same location i.e. multiple cores with exactly the same lat/lon pair
location_counts <- dat %>%
  group_by(longitude, latitude) %>%
  summarise(n = n(), .groups = "drop")
# plot
ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group), fill = "#eeeeee", color = "#a5a5a5") +
  geom_point(data = location_counts, aes(x = longitude, y = latitude, color = n), size = 1) +
  scale_color_viridis_c(name = "Sample count") +
  coord_cartesian(xlim = c(lon_min, lon_max), ylim = c(lat_min, lat_max)) +
  theme_minimal() +
  labs(x = "Longitude", y = "Latitude", title = "Spatial distribution of unique data points (colored by sample count)")
# ggsave("spatial_distribution.png", width = 10, height = 10)

# ============================================
# CROSS-VALIDATION COMPARISON: DIFFERENT TRAIN-TEST SCHEMAS
# ============================================
# This script compares the effect of different cross-validation methods on
# model performance for predicting median carbon density from remote sensing variables.
#
# Models tested:
#   - Random Forest (RF)
#   - XGBoost (boosted trees)
#   - Support Vector Machine (SVM)
#
# CV methods compared:
#   1. Random split (naive - ignores spatial structure)
#   2. Spatial block CV (accounts for spatial autocorrelation)
#   3. Core-level CV (non-spatial, but respects core structure)
#
# Response: Median carbon density per core (aggregated from samples)
# Predictors: Remote sensing variables (constant within cores)
# ============================================

# ============================================
# PREPARE DATA FOR MODELING
# ============================================

# Aggregate to core level: calculate median carbon density per core
# Remote sensing variables are constant within cores, so we use first() or mean()
core_data <- dat %>%
  group_by(random_core_variable) %>%
  summarise(
    # Response variable: median carbon density per core
    median_carbon_density = median(carbon_density_g_c_cm3, na.rm = TRUE),
    
    # Remote sensing predictors (constant within cores)
    KD_closest = first(KD_closest),
    RRS443_closest = first(RRS443_closest),
    wave_height_VHM0_p95_m_closest = first(wave_height_VHM0_p95_m_closest),
    po4_mean_1.5m_mmol_m3_closest = first(po4_mean_1.5m_mmol_m3_closest),
    pH_mean_1.5m_closest = first(pH_mean_1.5m_closest),
    bottomT_p95_C_closest = first(bottomT_p95_C_closest),
    vo_p90_1.5m_m_s_closest = first(vo_p90_1.5m_m_s_closest),
    uo_mean_1.5m_m_s_closest = first(uo_mean_1.5m_m_s_closest),
    Surf_fgco2_p95_molC_m2_yr_closest = first(Surf_fgco2_p95_molC_m2_yr_closest),
    
    # Additional predictors
    seagrass_species = first(seagrass_species),
    longitude = first(longitude),
    latitude = first(latitude),
    
    # Metadata
    n_samples = n(),
    .groups = "drop"
  ) %>%
  filter(!is.na(median_carbon_density))  # Remove cores with missing response

cat("Core-level dataset prepared:\n")
cat("  Number of cores:", nrow(core_data), "\n")
cat("  Number of samples:", sum(core_data$n_samples), "\n")
cat("  Median samples per core:", median(core_data$n_samples), "\n\n")

# Prepare predictor variables (remote sensing only, no depth or species for now)
predictor_vars <- c("KD_closest", "RRS443_closest", "wave_height_VHM0_p95_m_closest",
                    "po4_mean_1.5m_mmol_m3_closest", "pH_mean_1.5m_closest",
                    "bottomT_p95_C_closest", "vo_p90_1.5m_m_s_closest",
                    "uo_mean_1.5m_m_s_closest", "Surf_fgco2_p95_molC_m2_yr_closest")

# Check for missing values
missing_check <- core_data %>%
  dplyr::select(dplyr::all_of(predictor_vars), median_carbon_density) %>%
  dplyr::summarise_all(~sum(is.na(.)))

cat("Missing values per variable:\n")
print(missing_check)
cat("\n")

# Remove rows with missing values in predictors or response
# Use a simpler approach: create a temporary data frame with just the needed columns
temp_data <- core_data %>%
  dplyr::select(dplyr::all_of(predictor_vars), median_carbon_density)

# Get indices of complete cases
complete_indices <- complete.cases(temp_data)

# Filter the original data using these indices
core_data_complete <- core_data[complete_indices, ]

cat("Complete cases:", nrow(core_data_complete), "cores\n\n")

# ============================================
# CROSS-VALIDATION METHODS
# ============================================

# Method 1: Random split (naive approach - ignores spatial structure)
cat("=== METHOD 1: RANDOM SPLIT ===\n")
set.seed(42)
n_cores <- nrow(core_data_complete)
random_folds <- sample(rep(1:5, length.out = n_cores))
names(random_folds) <- 1:n_cores

# Method 2: Spatial block CV (using blockCV package)
cat("=== METHOD 2: SPATIAL BLOCK CV ===\n")

# Create sf object from core data
core_data_complete_spatial <- sf::st_as_sf(core_data_complete, 
                                           coords = c("longitude", "latitude"), 
                                           crs = 4326)

# Run spatial block CV
# Note: cv_spatial can work without a raster, but providing one ensures consistent block sizes
spatial_cv <- cv_spatial(x = core_data_complete_spatial,
                         k = 5,  # Number of folds
                         size = 10000,  # Size of blocks in meters (adjust based on your data spread)
                         selection = "random",  # Random blocks-to-fold assignment
                         iteration = 50,  # Number of iterations to find evenly dispersed folds
                         progress = TRUE,
                         biomod2 = FALSE,
                         hexagon = TRUE,
                         plot=FALSE
                         )  # Set to FALSE unless you need biomod2 format

# Extract fold assignments
spatial_folds <- spatial_cv$folds_ids

cat("  Spatial CV folds created.\n")


ecv <- cv_cluster(x = core_data_complete_spatial,
                  k = 5, 
                  scale = TRUE)
cv_similarity(cv=ecv, r=core_data_complete_spatial)
# Plot spatial blocks and fold numbers on world map

folds <- spatial_cv$folds_ids
blocks <- spatial_cv$blocks

highlight_fold <- 4

# Add fold information to location_counts by matching with core_data_complete
# Create a data frame with fold assignments for each location
core_folds_df <- core_data_complete %>%
  dplyr::select(longitude, latitude) %>%
  dplyr::mutate(fold = folds) %>%
  dplyr::group_by(longitude, latitude) %>%
  dplyr::summarise(fold = dplyr::first(fold), .groups = "drop") %>%
  dplyr::mutate(
    point_color = ifelse(fold == highlight_fold, "red", "black"),
  )

ggplot() +
  # World map as background
  geom_polygon(data = world, aes(x = long, y = lat, group = group), 
               fill = "#eeeeee", color = "#a5a5a5") +
  # Data points colored by fold (red for highlight_fold, black for others)
  geom_point(data = core_folds_df, 
             aes(x = longitude, y = latitude), 
             color = core_folds_df$point_color) +
  # Spatial folds polygons
  geom_sf(data = blocks, 
          aes(fill = as.factor(folds)),
          color = "red", alpha = 0.12, linewidth = 0.3, show.legend = FALSE) +
  # Fold number text labels
  geom_sf_text(data = blocks, 
               aes(label = folds), 
               size = 3, 
               fun.geometry = sf::st_centroid) +
  scale_size_continuous(name = "Sample count", range = c(0.5, 3)) +
  coord_sf(xlim = c(lon_min, lon_max), ylim = c(lat_min, lat_max)) +
  labs(x = "Longitude", y = "Latitude", title = paste("Spatial blocks (folds) on world map - Test fold", highlight_fold, "highlighted in red"), 
       fill = "Fold") +
  theme_minimal()

# plot the boxplot for the ensemble (all the model results together) as a function of the size of the cell
# test a range of spatial block sizes and evaluate model performance

cat("\n========================================\n")
cat("TESTING RANGE OF SPATIAL BLOCK SIZES\n")
cat("========================================\n\n")

# define range of block sizes to test (in meters)
block_sizes <- c(5000, 10000, 20000, 50000, 1000000)

# storage for results across all block sizes
all_size_results <- list()

# hyperparameter tuning options (same as above)
USE_HYPERPARAMETER_TUNING <- FALSE
NESTED_TUNING <- FALSE

# loop through each block size
for (block_size in block_sizes) {
  cat("\n--- Testing block size:", block_size, "meters ---\n")
  
  # run spatial block cv with this size
  spatial_cv_size <- tryCatch({
    cv_spatial(x = core_data_complete_spatial,
               k = 5,
               size = block_size,
               selection = "random",
               iteration = 50,
               progress = FALSE,
               biomod2 = FALSE,
               hexagon = FALSE,
               plot = FALSE)
  }, error = function(e) {
    cat("  Error creating spatial CV for size", block_size, ":", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(spatial_cv_size)) {
    cat("  Skipping block size", block_size, "\n")
    next
  }
  
  # extract fold assignments
  spatial_folds_size <- spatial_cv_size$folds_ids
  
  # run cv with this fold assignment
  results_size <- tryCatch({
    run_cv(paste0("spatial_split_size_", block_size), 
           spatial_folds_size, 
           core_data_complete, 
           predictor_vars,
           tune_hyperparams = USE_HYPERPARAMETER_TUNING,
           nested_tuning = NESTED_TUNING,
           verbose = FALSE)
  }, error = function(e) {
    cat("  Error running CV for size", block_size, ":", e$message, "\n")
    return(NULL)
  })
  
  if (!is.null(results_size)) {
    # add block size information to results
    results_size$block_size <- block_size
    all_size_results[[as.character(block_size)]] <- results_size
    cat("  Completed block size", block_size, "\n")
  }
}

# combine all results
if (length(all_size_results) > 0) {
  size_results_combined <- dplyr::bind_rows(all_size_results)
  
  # create boxplot showing ensemble performance (all models together) as a function of block size
  # calculate ensemble metrics: average across all models for each fold and block size
  ensemble_results <- size_results_combined %>%
    group_by(block_size, fold) %>%
    summarise(
      ensemble_r2 = mean(r2, na.rm = TRUE),
      ensemble_rmse = mean(rmse, na.rm = TRUE),
      ensemble_mae = mean(mae, na.rm = TRUE),
      ensemble_bias = mean(bias, na.rm = TRUE),
      .groups = "drop"
    )
  
  # plot r2 as a function of block size
  p_size_r2 <- ggplot(ensemble_results, aes(x = as.factor(block_size), y = ensemble_r2)) +
    geom_boxplot(fill = "lightblue", alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.5, color = "darkblue") +
    labs(x = "Block size (meters)", 
         y = "Ensemble R² (average across all models)", 
         title = "Model performance as a function of spatial block size") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  print(p_size_r2)
  ggsave("spatial_block_size_r2.png", width = 10, height = 6)
  
  # plot rmse as a function of block size
  p_size_rmse <- ggplot(ensemble_results, aes(x = as.factor(block_size), y = ensemble_rmse)) +
    geom_boxplot(fill = "lightcoral", alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.5, color = "darkred") +
    labs(x = "Block size (meters)", 
         y = "Ensemble RMSE (average across all models)", 
         title = "Model RMSE as a function of spatial block size") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  print(p_size_rmse)
  ggsave("spatial_block_size_rmse.png", width = 10, height = 6)
  
  # combined plot with both metrics
  ensemble_long <- ensemble_results %>%
    pivot_longer(cols = c(ensemble_r2, ensemble_rmse), 
                 names_to = "metric", 
                 values_to = "value") %>%
    mutate(metric = ifelse(metric == "ensemble_r2", "R²", "RMSE"))
  
  p_size_combined <- ggplot(ensemble_long, aes(x = as.factor(block_size), y = value)) +
    geom_boxplot(aes(fill = metric), alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    facet_wrap(~metric, scales = "free_y") +
    labs(x = "Block size (meters)", 
         y = "Ensemble metric value (average across all models)", 
         title = "Model performance metrics as a function of spatial block size",
         fill = "Metric") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "none")
  
  print(p_size_combined)
  ggsave("spatial_block_size_combined.png", width = 12, height = 6)
  
  # summary statistics by block size
  cat("\n=== Summary by block size ===\n")
  size_summary <- ensemble_results %>%
    group_by(block_size) %>%
    summarise(
      mean_r2 = mean(ensemble_r2, na.rm = TRUE),
      sd_r2 = sd(ensemble_r2, na.rm = TRUE),
      mean_rmse = mean(ensemble_rmse, na.rm = TRUE),
      sd_rmse = sd(ensemble_rmse, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(block_size)
  
  print(size_summary)
  
  # save results
  write.csv(size_results_combined, "spatial_block_size_results.csv", row.names = FALSE)
  write.csv(ensemble_results, "spatial_block_size_ensemble.csv", row.names = FALSE)
  write.csv(size_summary, "spatial_block_size_summary.csv", row.names = FALSE)
  
  cat("\nResults saved to:\n")
  cat("  - spatial_block_size_results.csv\n")
  cat("  - spatial_block_size_ensemble.csv\n")
  cat("  - spatial_block_size_summary.csv\n")
  cat("  - spatial_block_size_r2.png\n")
  cat("  - spatial_block_size_rmse.png\n")
  cat("  - spatial_block_size_combined.png\n\n")
} else {
  cat("\nNo results collected. Check for errors above.\n")
}

cat("========================================\n")
cat("SPATIAL BLOCK SIZE ANALYSIS COMPLETE\n")
cat("========================================\n\n")


# # Method 3: Environmental-only CV (split by environmental variables, but not spatially)
# cat("=== METHOD 3: ENVIRONMENTAL-ONLY CV ===\n")
# set.seed(42)
# environmental_cv <- cv_env(x = core_data_complete[, predictor_vars],
#                            k = 5,
#                            selection = "random",
#                            iteration = 100)

# ============================================
# RUN CROSS-VALIDATION FOR EACH METHOD
# ============================================


# Run CV for each method
cat("\n========================================\n")
cat("RUNNING CROSS-VALIDATION\n")
cat("========================================\n\n")

# Hyperparameter tuning options
USE_HYPERPARAMETER_TUNING <- FALSE  # Set to TRUE to enable hyperparameter tuning
NESTED_TUNING <- FALSE  # Set to TRUE for nested CV (tune per fold, more robust but slower)
                      # Set to FALSE to tune once on first fold (faster but less robust)

source("helpers.R")

tune_hyperparams <- USE_HYPERPARAMETER_TUNING
nested_tuning <- NESTED_TUNING

# Random split
results_random <- run_cv("random_split", random_folds, core_data_complete, predictor_vars, 
                         tune_hyperparams = USE_HYPERPARAMETER_TUNING, 
                         nested_tuning = NESTED_TUNING)

# Spatial block CV
results_spatial <- run_cv("spatial_split", spatial_folds, core_data_complete, predictor_vars,
                          tune_hyperparams = USE_HYPERPARAMETER_TUNING,
                          nested_tuning = NESTED_TUNING)

# Core-level CV
# results_env <- run_cv("env_split", env_folds, core_data_complete, predictor_vars)

# Combine all results
all_results <- dplyr::bind_rows(results_random, results_spatial)


# ============================================
# SUMMARIZE RESULTS
# ============================================

cat("\n========================================\n")
cat("RESULTS SUMMARY\n")
cat("========================================\n\n")

# Summary by method and model
summary_results <- all_results %>%
  group_by(method, model) %>%
  summarise(
    mean_r2 = mean(r2, na.rm = TRUE),
    sd_r2 = sd(r2, na.rm = TRUE),
    mean_rmse = mean(rmse, na.rm = TRUE),
    sd_rmse = sd(rmse, na.rm = TRUE),
    mean_mae = mean(mae, na.rm = TRUE),
    sd_mae = sd(mae, na.rm = TRUE),
    mean_bias = mean(bias, na.rm = TRUE),
    sd_bias = sd(bias, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(method, model)

print(summary_results)

# ============================================
# VISUALIZE RESULTS
# ============================================

# Plot R² comparison
p1 <- ggplot(all_results, aes(x = method, y = r2, fill = model)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.3) +
  labs(x = "CV Method", y = "R²", title = "Model Performance by CV Method") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5), plot.title = element_text(hjust = 0.5)
  )

print(p1)
ggsave("cv_comparison_r2.png", width = 10, height = 6)

# Plot RMSE comparison
p2 <- ggplot(all_results, aes(x = method, y = rmse, fill = model)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.3) +
  labs(x = "CV Method", y = "RMSE", title = "RMSE by CV Method") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5), plot.title = element_text(hjust = 0.5)
  )
print(p2)
ggsave("cv_comparison_rmse.png", width = 10, height = 6)

# Detailed comparison plot
# Facet R² and RMSE; label left facet as "R²", right facet as "RMSE"
metric_names <- c(r2 = "R²", rmse = "RMSE")
p3 <- all_results %>%
  dplyr::select(dplyr::all_of(c("method", "model", "r2", "rmse", "mae", "bias"))) %>%
  pivot_longer(cols = c(r2, rmse), names_to = "metric", values_to = "value") %>%
  ggplot(aes(x = method, y = value, fill = model)) +
  geom_boxplot() +
  facet_wrap(
    ~metric,
    scales = "free_y",
    labeller = as_labeller(metric_names)
  ) +
  labs(x = "CV Method", y = NULL, title = "Model Performance Metrics") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5),
  )

print(p3)
ggsave("cv_comparison_all_metrics.png", width = 12, height = 6)

# ============================================
# SAVE RESULTS
# ============================================

write.csv(all_results, "cv_results_detailed.csv", row.names = FALSE)
write.csv(summary_results, "cv_results_summary.csv", row.names = FALSE)

cat("\nResults saved to:\n")
cat("  - cv_results_detailed.csv\n")
cat("  - cv_results_summary.csv\n")
cat("  - cv_comparison_r2.png\n")
cat("  - cv_comparison_rmse.png\n")
cat("  - cv_comparison_all_metrics.png\n\n")

cat("========================================\n")
cat("ANALYSIS COMPLETE\n")
cat("========================================\n")

