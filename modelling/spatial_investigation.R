### Investigate effect of different cross-validation methods on model performance

## ================================ SETUP ================================
# refresh andset working directory
rm(list = ls())
setwd(here::here())

# source helper functions
source("modelling/helpers.R")

load_packages(c("here", "mgcv", "tidyverse", "ggplot2", "visreg", "randomForest", "corrplot", "blockCV"))
set.seed(42)

## get data
For_modeling_df_shallow_carbon_density <- read_rds("data/For_modeling_df_shallow_carbon_density.rds")
summary(For_modeling_df_shallow_carbon_density)
dat <- For_modeling_df_shallow_carbon_density

## ================================ CROSS-VALIDATION COMPARISON: DIFFERENT TRAIN-TEST SCHEMAS ================================
# Models tested:
#   - Random Forest (RF)
#   - XGBoost (boosted trees)
#   - Support Vector Machine (SVM)
#   - Neural Network (NN)
#   - Gaussian Process Regression (GPR)

# CV methods compared:
#   1. Random split (naive - ignores spatial structure)
#   2. Spatial block CV (accounts for spatial autocorrelation)
#   3. Environmental-only CV (split by environmental variables, but not spatially)
#
# Response: Median carbon density per core (aggregated from samples)
# Predictors: Remote sensing variables (constant within cores)

## aggregate data to core level. Currently only returns a subset of the environmental variables for speed. TODO: add all environmental variables.
core_data <- dat %>%
  group_by(random_core_variable) %>%
  summarise(median_carbon_density = median(carbon_density_g_c_cm3, na.rm = TRUE),
  # remote sensing predictors
  KD_closest = first(KD_closest),
  RRS443_closest = first(RRS443_closest),
  wave_height_VHM0_p95_m_closest = first(wave_height_VHM0_p95_m_closest),
  po4_mean_1.5m_mmol_m3_closest = first(po4_mean_1.5m_mmol_m3_closest),
  pH_mean_1.5m_closest = first(pH_mean_1.5m_closest),
  bottomT_p95_C_closest = first(bottomT_p95_C_closest),
  vo_p90_1.5m_m_s_closest = first(vo_p90_1.5m_m_s_closest),
  uo_mean_1.5m_m_s_closest = first(uo_mean_1.5m_m_s_closest),
  Surf_fgco2_p95_molC_m2_yr_closest = first(Surf_fgco2_p95_molC_m2_yr_closest),
  
  # additional predictors
  seagrass_species = first(seagrass_species),
  longitude = first(longitude),
  latitude = first(latitude),
  
  # metadata
  n_samples = n(),
  .groups = "drop"
  ) %>%
  filter(!is.na(median_carbon_density))  # Remove cores with missing response

cat("Core-level dataset prepared:\n")
cat("  Number of cores:", nrow(core_data), "\n")
cat("  Number of samples:", sum(core_data$n_samples), "\n")
cat("  Median samples per core:", median(core_data$n_samples), "\n\n")
cat("  Number of environmental variables cores:", length(colnames(core_data)) - 2, "\n\n")

# prepare predictor variables (remote sensing only, no depth or species for now)
predictor_vars <- c("KD_closest", "RRS443_closest", "wave_height_VHM0_p95_m_closest",
                    "po4_mean_1.5m_mmol_m3_closest", "pH_mean_1.5m_closest",
                    "bottomT_p95_C_closest", "vo_p90_1.5m_m_s_closest",
                    "uo_mean_1.5m_m_s_closest", "Surf_fgco2_p95_molC_m2_yr_closest")

# clean dataframe of missing values (remove rows with missing values in predictors or response)
core_data_complete <- core_data %>%
  dplyr::select(dplyr::all_of(predictor_vars), median_carbon_density) %>%
  dplyr::filter(complete.cases(.))

cat("Complete cases:", nrow(core_data_complete), "cores\n\n")

## ================================ CROSS-VALIDATION METHODS ================================
n_folds <- 5
n_cores <- nrow(core_data_complete)
# pick which fold to highlight (this one used for testing, all else for training)
highlight_fold <- 1

# Method 1: Random split (naive approach - ignores spatial structure)
cat("=== METHOD 1: RANDOM SPLIT ===\n")
random_folds <- sample(rep(1:n_folds, length.out = n_cores))
names(random_folds) <- 1:n_cores

# assign train/test: randomly split, e.g., Fold 1 is test, rest are train
core_data$cv_set <- ifelse(random_folds == highlight_fold, "Test", "Train")
world <- map_data("world")
lon_min <- min(core_data$longitude, na.rm = TRUE)
lon_max <- max(core_data$longitude, na.rm = TRUE)
lat_min <- min(core_data$latitude, na.rm = TRUE)
lat_max <- max(core_data$latitude, na.rm = TRUE)

ggplot() +
  # world map background (added first so it stays in the back)
  geom_polygon(
    data = world,
    aes(x = long, y = lat, group = group),
    fill = "#eeeeee",
    color = "#a5a5a5",
    inherit.aes = FALSE
  ) +
  # train (black) and test (red)
  geom_point(
    data = core_data,
    aes(x = longitude, y = latitude, color = cv_set),
    size = 2, alpha = 0.8
  ) +
  scale_color_manual(values = c("Train" = "black", "Test" = "red"), name = NULL) +
  theme(legend.position = "bottom",
) +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "Random split folds: train (black) vs test (blue)"
  ) +
  coord_cartesian(xlim = c(lon_min, lon_max), ylim = c(lat_min, lat_max)) 
ggsave("figures/random_split_folds.png", width = 10, height = 6)



# Method 2: Spatial block CV (using blockCV package)
cat("=== METHOD 2: SPATIAL BLOCK CV ===\n")

## visualise spatial blocking schema
core_data_complete_spatial <- sf::st_as_sf(core_data_complete, 
                                           coords = c("longitude", "latitude"), 
                                           crs = 4326)
spatial_cv <- cv_spatial(x = core_data_complete_spatial,
                         k = n_folds,  # number of folds
                         size = 10000,  # size of blocks in meters (adjust based on your data spread)
                         selection = "random",  # random blocks-to-fold assignment
                         iteration = 50,  # number of iterations to find evenly dispersed folds
                         progress = FALSE,
                         biomod2 = FALSE,
                         hexagon = TRUE,
                         plot=FALSE)

spatial_folds <- spatial_cv$folds_ids
folds <- spatial_cv$folds_ids
blocks <- spatial_cv$blocks

# dataframe with fold assignment for each core
core_folds_df <- core_data_complete %>%
  dplyr::select(longitude, latitude) %>%
  dplyr::mutate(fold = folds) %>%
  dplyr::group_by(longitude, latitude) %>%
  dplyr::summarise(fold = dplyr::first(fold), .groups = "drop") %>%
  dplyr::mutate(
    point_color = ifelse(fold == highlight_fold, "red", "black"),
  )
# plot distribution and fold assignment
ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group), 
               fill = "#eeeeee", color = "#a5a5a5") +
  # red for highlight_fold, black for others
  geom_point(data = core_folds_df, 
             aes(x = longitude, y = latitude), 
             color = core_folds_df$point_color) +
  # spatial folds polygons
  geom_sf(data = blocks, 
          aes(fill = as.factor(folds)),
          color = "red", alpha = 0.12, linewidth = 0.3, show.legend = FALSE) +
  # fold number text labels
  geom_sf_text(data = blocks, 
               aes(label = folds), 
               size = 3, 
               fun.geometry = sf::st_centroid) +
  scale_size_continuous(name = "Sample count", range = c(0.5, 3)) +
  coord_sf(xlim = c(lon_min, lon_max), ylim = c(lat_min, lat_max)) +
  labs(x = "Longitude", y = "Latitude", title = paste("Spatial blocks (folds) on world map - Test fold", highlight_fold, "highlighted in red"), 
       fill = "Fold") +
  theme_minimal()


## ================================ TESTING RANGE OF SPATIAL BLOCK SIZES ================================

# define range of block sizes to test (in meters). TODO: be more rigorous by quantifying spatial similarity from rasters
block_sizes <- c(5000, 10000, 20000, 50000, 1000000)

all_size_results <- list()

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


# TODO: Method 3: Environmental-only CV (split by environmental variables, but not spatially)
# cat("=== METHOD 3: ENVIRONMENTAL-ONLY CV ===\n")
# ecv <- cv_cluster(x = core_data_complete_spatial,
#                   k = 5, 
#                   scale = TRUE)
# cv_similarity(cv=ecv, r=core_data_complete_spatial)
# set.seed(42)
# environmental_cv <- cv_env(x = core_data_complete[, predictor_vars],
#                            k = 5,
#                            selection = "random",
#                            iteration = 100)


## ================================ RUN CROSS-VALIDATION FOR EACH METHOD ================================

# hyperparameter tuning options
USE_HYPERPARAMETER_TUNING <- FALSE  # set to TRUE to enable hyperparameter tuning
NESTED_TUNING <- FALSE  # set to TRUE for nested CV (tune per fold, more robust but slower)
                      # FALSE means tuning once on first fold (faster but less robust)


# random split
results_random <- run_cv("random_split", random_folds, core_data_complete, predictor_vars, 
                         tune_hyperparams = USE_HYPERPARAMETER_TUNING, 
                         nested_tuning = NESTED_TUNING)
# spatial block CV
results_spatial <- run_cv("spatial_split", spatial_folds, core_data_complete, predictor_vars,
                          tune_hyperparams = USE_HYPERPARAMETER_TUNING,
                          nested_tuning = NESTED_TUNING)
# environmental CV
# results_env <- run_cv("env_split", env_folds, core_data_complete, predictor_vars)

all_results <- dplyr::bind_rows(results_random, results_spatial)

# ================================ SUMMARIZE AND VISUALIZE RESULTS ================================

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


# plot R2 and RMSE together
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

# save results
write.csv(all_results, "cv_results_detailed.csv", row.names = FALSE)
write.csv(summary_results, "cv_results_summary.csv", row.names = FALSE)




# ================================ DEPRECATED ================================
# # Plot R² comparison
# p1 <- ggplot(all_results, aes(x = method, y = r2, fill = model)) +
#   geom_boxplot() +
#   geom_jitter(width = 0.2, alpha = 0.3) +
#   labs(x = "CV Method", y = "R²", title = "Model Performance by CV Method") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 0, hjust = 0.5), plot.title = element_text(hjust = 0.5)
#   )

# print(p1)
# ggsave("cv_comparison_r2.png", width = 10, height = 6)

# # Plot RMSE comparison
# p2 <- ggplot(all_results, aes(x = method, y = rmse, fill = model)) +
#   geom_boxplot() +
#   geom_jitter(width = 0.2, alpha = 0.3) +
#   labs(x = "CV Method", y = "RMSE", title = "RMSE by CV Method") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 0, hjust = 0.5), plot.title = element_text(hjust = 0.5)
#   )
# print(p2)
# ggsave("cv_comparison_rmse.png", width = 10, height = 6)