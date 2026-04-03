# sanity.R — diagnostic: how many unique covariate vectors exist?
#
# If the coarse raster resolution (~4 km) maps many distinct (lat, lon)
# locations to the same pixel, the covariate vectors will be identical and
# location-grouped folds won't fully prevent leakage.

# TODO: measure variance/similarity of covariate vectors across train and test for different folds under different seeds
setwd(here::here())

library(ggplot2)
library(maps)
library(tidyverse)
source("modelling/R/extract_covariates_from_rasters.R")
source("modelling/config/pipeline_config.R")

dat <- readr::read_rds("data/all_extracted_new.rds")
source("sanity_helpers.R")
dat <- process_rs_covariates(dat)
source("modelling/R/helpers.R")

target_var <- "median_carbon_density_100cm"

library(ncdf4)
library(raster)
# load gam model
cfg <- get_pipeline_config()
cv_regime_name <- cfg$cv_regime_name
final_models_dir <- file.path("output", cv_regime_name, "final_models")
gam_model <- readRDS(file.path(final_models_dir, "GAM_final.rds"))


# plot maps of the predictor variables
pred_vars <- gam_model$predictor_vars
for (pred_var in pred_vars) {
  fp <- file.path("data/env_rasters", paste0(pred_var, ".nc"))
  raster_data <- raster::raster(fp) # read .nc file directly with raster
  raster_points <- raster::rasterToPoints(raster_data)
  raster_points <- as.data.frame(raster_points)
  colnames(raster_points) <- c("longitude", "latitude", "value")
  p <- ggplot(raster_points, aes(x = longitude, y = latitude, fill = value)) +
    geom_raster() +
    theme_minimal() +
    labs(title = pred_var)
  ggsave(
    file.path(final_models_dir, "GAM_final", paste0(pred_var, ".png")),
    plot = p
  )
}

for (model in c("GAM", "XGB", "GPR")) {
  obj <- readRDS(file.path(final_models_dir, paste0(model, "_final.rds")))
  print(toupper(model))
  print(obj$hyperparams)
}


# === Seagrass map ===
# open csv of seagrass
fp <- "data/seagrass_eov_poly_2025.csv"
seagrass_eov_poly_2025 <- read.csv(fp)

# === Convert model parameter choice to nice readable versions === 

model_list <- c("GPR", "GAM", "XGB", "LR")
if (!exists("robust_fold_seed_list", inherits = TRUE)) {
  source("modelling/config/pipeline_config.R")
  robust_fold_seed_list <- get_pipeline_config()$robust_fold_seed_list
}
model_covariates <- read.csv(sprintf("output/pixel_grouped/covariate_selection/robust_pixel_grouped/pruned_model_variables_shap_robust_pixel_grouped_seeds_%s.csv", paste(robust_fold_seed_list, collapse = "-")))

# use plot_config.R to map the model covariates to human readable names
source("modelling/R/plot_config.R")
# Goal: one row per selected variable with neat label and number of models sharing it.
library(dplyr)

# Write to csv
out_dir_labels <- sprintf("output/pixel_grouped/covariate_selection/robust_pixel_grouped/seeds_%s", paste(robust_fold_seed_list, collapse = "-"))
dir.create(out_dir_labels, recursive = TRUE, showWarnings = FALSE)

out_df <- model_covariates %>%
  dplyr::mutate(variable = as.character(variable), model = as.character(model)) %>%
  dplyr::distinct(model, variable) %>%
  dplyr::group_by(variable) %>%
  dplyr::summarise(
    n_models_shared = dplyr::n_distinct(model),
    models = paste(sort(unique(model)), collapse = ", "),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    label = unname(label_vars(variable))
  ) %>%
  dplyr::select(variable, label, n_models_shared, models) %>%
  dplyr::arrange(dplyr::desc(n_models_shared), variable)

write.csv(out_df, file.path(out_dir_labels, "model_covariates_shared_counts.csv"), row.names = FALSE)

# === Covariate uniqueness diagnostic ===

exclude_regions <- c("Black Sea")
if (!"region" %in% names(dat)) dat <- assign_region_from_latlon(dat)
dat <- dat[is.na(dat$region) | !dat$region %in% exclude_regions, ]

covars <- raster_covariates[raster_covariates %in% colnames(dat)]

cat("=== Covariate uniqueness diagnostic ===\n\n")
cat("Total rows:             ", nrow(dat), "\n")

# Unique (lat, lon) locations
loc_key <- paste(dat$longitude, dat$latitude)
n_unique_locs <- length(unique(loc_key))
cat("Unique (lon, lat) pairs:", n_unique_locs, "\n")

# Unique covariate vectors (raster-derived columns only)
cov_mat <- dat[, covars, drop = FALSE]
cov_key <- do.call(paste, c(cov_mat, sep = "|"))
n_unique_cov <- length(unique(cov_key))
cat("Unique covariate vectors:", n_unique_cov, "\n")
cat("Ratio (unique covs / unique locs):", round(n_unique_cov / n_unique_locs, 3), "\n\n")

# How many locations share a covariate vector with at least one other location?
loc_to_cov <- data.frame(loc = loc_key, cov = cov_key, stringsAsFactors = FALSE)
loc_to_cov <- loc_to_cov[!duplicated(loc_to_cov), ]  # one row per unique location
cov_counts <- table(loc_to_cov$cov)
n_locs_sharing <- sum(cov_counts[match(loc_to_cov$cov, names(cov_counts))] > 1)
cat("Unique locations that share covariates with another location:",
    n_locs_sharing, "/", n_unique_locs,
    sprintf("(%.1f%%)\n", 100 * n_locs_sharing / n_unique_locs))

    

# count number of unique combinations of remote sensing variables
dat_unique <- dat %>%
  select(all_of(raster_covariates)) %>%
  unique()
dim(dat_unique)

# Calculate the mean and standard deviation of the carbon density column (median_carbon_density_100cm)
# at each unique combination of remote sensing variables
dat_clustered <- dat %>%
  group_by(across(all_of(raster_covariates))) %>%
  summarise(
    mean_carbon_density = mean(median_carbon_density_100cm, na.rm = TRUE),
    sd_carbon_density = sd(median_carbon_density_100cm, na.rm = TRUE),
    n_observations = n(),
    .groups = "drop"
  )
dat_clustered$n_observations


# assign each row in dat to a unique id based on the unique combination of remote sensing variables
dat$id <- match(
  do.call(paste, c(dat[, raster_covariates, drop = FALSE], sep = "|")),
  do.call(paste, c(dat_clustered[, raster_covariates, drop = FALSE], sep = "|"))
)

# === R2 craziness ===

# load by_seed_detailed for largest seed sweep
by_seed_detailed <- read.csv("output/pixel_grouped/cv_pipeline/robust_pixel_grouped_random_evaluation_robustSeeds_42-43-44-45-46-47-48-49-50/by_seed_detailed.csv")

# plot a histogram of r2 values
ggplot(by_seed_detailed, aes(x = r2)) +
  geom_histogram(binwidth = 0.1) +
  theme_minimal() +
  labs(x = "Mean R2", y = "Count")


# average r2 values over seeds and folds (the poor performance)
mean_r2 <- by_seed_detailed %>%
  group_by(model) %>%
  summarise(mean_r2 = mean(r2, na.rm = TRUE))
ggplot(mean_r2, aes(x = model, y = mean_r2)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(x = "Model", y = "Mean R2")




# reconstruct the pooled r2 values
pooled_reconstructed <- by_seed_detailed %>%
  group_by(model, fold_seed) %>%
  summarise(
    N = sum(n_eval, na.rm = TRUE),
    SS_res = sum(ss_residual, na.rm = TRUE),
    S1 = sum(sum_y, na.rm = TRUE),
    S2 = sum(sum_y2, na.rm = TRUE),
    SS_tot = S2 - (S1^2 / N),
    pooled_r2_reconstructed = ifelse(N >= 2 & SS_tot > 0, 1 - SS_res / SS_tot, NA_real_),
    .groups = "drop"
  )
# average pooled_r2_reconstructed over seeds
pooled_reconstructed %>%
  group_by(model) %>%
  summarise(mean_pooled_r2_reconstructed = mean(pooled_r2_reconstructed, na.rm = TRUE))




# Compute per-model pooled R2 using (1 - sum(ss_residual) / sum(ss_total)), grouped by model
pooled_r2_by_model <- by_seed_detailed %>%
  group_by(model) %>%
  summarise(
    pooled_r2 = 1 - (sum(ss_residual, na.rm = TRUE) / sum(ss_total, na.rm = TRUE))
  )
print(pooled_r2_by_model)


# Compute per-model pooled R2 within each seed, then average across seeds
pooled_r2_by_seed <- by_seed_detailed %>%
  group_by(model, fold_seed) %>%
  summarise(
    pooled_r2 = 1 - (sum(ss_residual, na.rm = TRUE) / sum(ss_total, na.rm = TRUE)),
    .groups = "drop"
  )

# Now average the pooled R2 over the seeds for each model
pooled_r2_by_model_avg <- pooled_r2_by_seed %>%
  group_by(model) %>%
  summarise(
    mean_pooled_r2 = mean(pooled_r2, na.rm = TRUE),
    sd_pooled_r2   = sd(pooled_r2, na.rm = TRUE),
    n_seeds        = dplyr::n(),
    .groups = "drop"
  )

print(pooled_r2_by_model_avg)



### splitting regimes
split_dat_list <- list()
for (i in 1:100) {
  split_dat_list[[i]] <- cov_id_split(dat, seed=i)
}

# Used by the optional covariate-vector-level train/test diagnostics below.
# `dat$id` corresponds to the row index in `dat_clustered` (via `match`), so
# these can be used to subset `dat_clustered` directly.
train_cov_ids <- sort(unique(split_dat_list[[1]]$train$id))

ratio_list <- sapply(split_dat_list, function(x) x$ratio)
ggplot(data.frame(ratio = ratio_list), aes(x = ratio)) +
  geom_histogram(binwidth = 0.01) +
  theme_minimal() +
  labs(x = "Train/test split ratio", y = "Count")

# on a grid, plot the first n train/test splits spatially
plot_n_splits_on_grid(dat, split_dat_list, n = 9, world_map, min_lon, max_lon, min_lat, max_lat)



# k-fold cross-validation on unique covariate vectors (i.e., `dat$id` clusters)
set.seed(42)
unique_cov_ids <- sort(unique(na.omit(dat$id)))
k_folds <- 5L
fold_assignment <- data.frame(
  id = unique_cov_ids,
  fold = as.integer(sample(rep(seq_len(k_folds), length.out = length(unique_cov_ids))))
)
dat <- dat %>%
  left_join(fold_assignment, by = "id")

# visualise the fold assignments spatially. N.B
ggplot(dat, aes(x = longitude, y = latitude, color = as.factor(fold))) +
  geom_point() +
  theme_minimal() +
  labs(x = "Longitude", y = "Latitude", color = "Fold")

# for each fold as the test, plot the train/test split ratio
ratio_list <- list()
for (i in 1:k_folds) {
  test_ids <- dat$id[dat$fold == i]
  train_ids <- dat$id[dat$fold != i]
  test_dat <- dat[dat$id %in% test_ids, ]
  train_dat <- dat[dat$id %in% train_ids, ]
  ratio_list[[i]] <- length(train_ids) / (length(train_ids) + length(test_ids))
}

ggplot(data.frame(ratio = unlist(ratio_list)), aes(x = ratio)) +
  geom_histogram(binwidth = 0.01) +
  theme_minimal() +
  labs(x = "Train/test split ratio", y = "Count")

# TODO: visualise environmental clusters for each fold


# Simplest predictor set: all raster covariates plus species factor if present.
predictor_vars <- raster_covariates[raster_covariates %in% names(dat)]
predictor_vars <- prune_by_correlation(dat, predictor_vars, target_var,
  cor_threshold = 0.8)  # filter vars by correlation with target


# === compute maximum possible R2 ===

result <- estimate_r2_upper_bound(dat, predictor_vars, target_var)
result$r2_max

between_var <- var(result$stats$mean_y, na.rm = TRUE)

within_var <- result$noise_variance
total_var  <- result$total_variance

c(
  within = within_var,
  between = between_var,
  total = total_var
)

mean_within_var <- mean(result$stats$var_y[result$stats$n > 1], na.rm = TRUE)


ggplot(result$stats, aes(x = n, y = var_y)) +
  geom_point() +
  scale_y_log10() +
  labs(
    x = "Number of replicates per covariate vector",
    y = "Variance of carbon density",
    title = "Within-group variance (noise) across environments"
  )


# checking distances between covariate vectors

# scale predictor variables
scale_params <- compute_scale_params(dat[, predictor_vars], predictor_vars)
dat_scaled <- apply_scaling(dat[, predictor_vars], scale_params, predictor_vars)

# get the variance of each column of the dat_scaled dataframe (should all be 1)
var_list <- list()
for (i in 1:ncol(dat_scaled)) {
  var_list[[i]] <- var(dat_scaled[, i])
}
var_df <- data.frame(var = unlist(var_list), colnames(dat_scaled))
summary(var_df$var)



X <- dat_scaled[, predictor_vars]
dists <- dist(X)
summary(dists)

# plot the distances between covariate vectors
ggplot(data.frame(dists = dists), aes(x = dists)) +
  geom_histogram(binwidth = 0.01) +
  theme_minimal() +
  labs(x = "Distance between covariate vectors", y = "Count")

library(FNN)

nn <- get.knn(X, k = 2)
summary(nn$nn.dist[,2])

# plot distribution of the distances between the second nearest neighbour and the first nearest neighbour
ggplot(data.frame(dists = nn$nn.dist[,2]), aes(x = dists)) +
  geom_histogram(binwidth = 0.01) +
  theme_minimal() +
  labs(x = "Distance between second nearest neighbour and first nearest neighbour", y = "Count")


if ("seagrass_species" %in% names(dat)) predictor_vars <- unique(c(predictor_vars, "seagrass_species"))
predictor_vars <- setdiff(predictor_vars, target_var)
if (length(predictor_vars) < 1L) stop("No predictor variables found for linear regression.")

# Model suite over grouped vs row-wise random CV, and with vs without log target

# Suite configuration (used by all subsequent CV variants)
covariate_vars_suite <- setdiff(predictor_vars, "seagrass_species")
suite_species_var <- if ("seagrass_species" %in% names(dat)) "seagrass_species" else NULL

# === CV suite: grouped covariate-vector folds (leakage control) ===
cat("\n=== CV suite: grouped covariate-vector folds (original target) ===\n")
suite_results <- run_cv_simple_model_suite(
  dat = dat,
  target_var = target_var,
  covariate_vars = covariate_vars_suite,
  species_var = suite_species_var,
  fold_col = "fold",
  rf_ntree = 200L,
  xgb_nrounds = 150L,
  gpr_kernel = "matern52",
  inverse_response_for_metrics = FALSE
)
print(suite_results$summary)
if (!is.null(suite_results$per_fold) && "status" %in% names(suite_results$per_fold)) {
  print(dplyr::count(suite_results$per_fold, model, status))
}

cat("\n=== CV suite: grouped covariate-vector folds (log-transformed target) ===\n")
dat_log <- dat
dat_log <- transform_response(dat_log, response_var = target_var, log = TRUE)
suite_results_log <- run_cv_simple_model_suite(
  dat = dat_log,
  target_var = target_var,
  covariate_vars = covariate_vars_suite,
  species_var = suite_species_var,
  fold_col = "fold",
  rf_ntree = 200L,
  xgb_nrounds = 150L,
  gpr_kernel = "matern52",
  inverse_response_for_metrics = TRUE
)
print(suite_results_log$summary)
if (!is.null(suite_results_log$per_fold) && "status" %in% names(suite_results_log$per_fold)) {
  print(dplyr::count(suite_results_log$per_fold, model, status))
}

# Sanity check: row-wise random folds (identical covariate vectors can land in both train/test).
# This should typically produce higher apparent performance due to overlap.
set.seed(123)
dat_random <- dat
dat_random$fold_random <- as.integer(
  sample(rep(seq_len(k_folds), length.out = nrow(dat_random)))
)

cat("\n=== CV suite: row-wise random folds ===\n")
suite_results_random <- run_cv_simple_model_suite(
  dat = dat_random,
  target_var = target_var,
  covariate_vars = covariate_vars_suite,
  species_var = suite_species_var,
  fold_col = "fold_random",
  rf_ntree = 200L,
  xgb_nrounds = 150L,
  gpr_kernel = "matern52",
  inverse_response_for_metrics = FALSE
)
print(suite_results_random$summary)
if (!is.null(suite_results_random$per_fold) && "status" %in% names(suite_results_random$per_fold)) {
  print(dplyr::count(suite_results_random$per_fold, model, status))
}

cat("\n=== CV suite: log-transformed target (row-wise random folds) ===\n")
dat_random_log <- dat_random
dat_random_log <- transform_response(dat_random_log, response_var = target_var, log = TRUE)

suite_results_random_log <- run_cv_simple_model_suite(
  dat = dat_random_log,
  target_var = target_var,
  covariate_vars = covariate_vars_suite,
  species_var = suite_species_var,
  fold_col = "fold_random",
  rf_ntree = 200L,
  xgb_nrounds = 150L,
  gpr_kernel = "matern52",
  inverse_response_for_metrics = TRUE
)
print(suite_results_random_log$summary)
if (!is.null(suite_results_random_log$per_fold) && "status" %in% names(suite_results_random_log$per_fold)) {
  print(dplyr::count(suite_results_random_log$per_fold, model, status))
}


# run cv suite across multiple seeds to check the stability of the results
run_cv_model_suite_over_seed_list(dat_random_log, target_var, covariate_vars_suite, suite_species_var, seed_list = 1:10, inverse_response_for_metrics = TRUE)



### VISUALISATION

# For each unique id, plot the distribution of the carbon density variable for that id
# Only include IDs with at least 2 observations for meaningful boxplots
id_counts <- as.data.frame(table(dat$id))
valid_ids <- id_counts$Var1[id_counts$Freq >= 2]
filtered_dat <- dat[dat$id %in% valid_ids, ]

# Order ids (clusters) by median carbon density for plotting
id_medians <- aggregate(filtered_dat[[target_var]], by = list(filtered_dat$id), FUN = median, na.rm = TRUE)
ordered_ids <- id_medians$Group.1[order(id_medians$x)]
filtered_dat$id <- factor(filtered_dat$id, levels = ordered_ids)

# Get n observations per ID (for color)
obs_per_id <- aggregate(filtered_dat[[target_var]], by = list(id = filtered_dat$id), FUN = length)
colnames(obs_per_id)[2] <- "n_obs"

# Join n_obs onto filtered_dat
filtered_dat <- merge(filtered_dat, obs_per_id, by.x = "id", by.y = "id", all.x = TRUE, sort = FALSE)

# Create the plot: boxplots colored by n_obs
p <- ggplot(filtered_dat, aes(x = id, y = .data[[target_var]], fill = n_obs)) +
  geom_boxplot() +
  scale_fill_viridis_c(option = "plasma", name = "n_obs") +
  theme_minimal() +
  labs(x = "Covariate vector ID", y = "Measured Carbon Density") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
print(p)


# do a pca on the covariate vectors
pca <- prcomp(dat_unique[, raster_covariates], scale. = TRUE)
summary(pca)
# plot the cumulative variance with lines at 90, 95, and 99%
cumulative_variance <- cumsum(pca$sdev^2) / sum(pca$sdev^2)
thresholds <- c(0.9, 0.95, 0.99)
threshold_colors <- c("red", "blue", "green")

# Find the index of the first component reaching or exceeding each threshold
threshold_indices <- sapply(thresholds, function(t) which(cumulative_variance >= t)[1])

df <- data.frame(x = 1:ncol(pca$x), y = cumulative_variance)

p <- ggplot(df, aes(x = x, y = y)) +
  geom_line() +
  # Add horizontal lines
  geom_hline(yintercept = thresholds[1], color = threshold_colors[1], linetype = "dashed") +
  geom_hline(yintercept = thresholds[2], color = threshold_colors[2], linetype = "dashed") +
  geom_hline(yintercept = thresholds[3], color = threshold_colors[3], linetype = "dashed") +
  # Add vertical lines where index meets threshold
  geom_vline(xintercept = threshold_indices[1], color = threshold_colors[1], linetype = "dotted") +
  geom_vline(xintercept = threshold_indices[2], color = threshold_colors[2], linetype = "dotted") +
  geom_vline(xintercept = threshold_indices[3], color = threshold_colors[3], linetype = "dotted") +
  # Add text labels at intersection points
  annotate("text", 
           x = threshold_indices,
           y = thresholds + 0.03,
           label = threshold_indices,
           color = threshold_colors,
           angle = 90,
           vjust = 0,
           size = 3.5) +
  labs(x = "Principal Component",
       y = "Cumulative Variance Explained") +
  theme_minimal()

print(p)


# do a pca on the correlated covariates and plot the cumulative variance
pca <- prcomp(dat_unique[, predictor_vars], scale. = TRUE)
summary(pca)
cumulative_variance <- cumsum(pca$sdev^2) / sum(pca$sdev^2)
thresholds <- c(0.9, 0.95, 0.99)
threshold_colors <- c("red", "blue", "green")
threshold_indices <- sapply(thresholds, function(t) which(cumulative_variance >= t)[1])
df <- data.frame(x = 1:ncol(pca$x), y = cumulative_variance)
p <- ggplot(df, aes(x = x, y = y)) +
  geom_line() +
  geom_hline(yintercept = thresholds[1], color = threshold_colors[1], linetype = "dashed") +
  geom_hline(yintercept = thresholds[2], color = threshold_colors[2], linetype = "dashed") +
  geom_hline(yintercept = thresholds[3], color = threshold_colors[3], linetype = "dashed") +
  geom_vline(xintercept = threshold_indices[1], color = threshold_colors[1], linetype = "dotted") +
  geom_vline(xintercept = threshold_indices[2], color = threshold_colors[2], linetype = "dotted") +
  geom_vline(xintercept = threshold_indices[3], color = threshold_colors[3], linetype = "dotted") +
  annotate("text", 
           x = threshold_indices,
           y = thresholds + 0.03,
           label = threshold_indices,
           color = threshold_colors,
           angle = 90,
           vjust = 0,
           size = 3.5) +
  labs(x = "Principal Component",
       y = "Cumulative Variance Explained") +
  theme_minimal()
print(p)


# Plot the signed correlation (not absolute) of each PCA component with the target variable,
# using the PCA fitted on the correlation-pruned predictors stored in `pca`.
# using color to indicate the sign (blue for positive, red for negative).
# Compute correlation of each PC (with correct names) with the target variable
#
# `dat_unique` contains only raster covariates, so it does not include the
# target column. Align a unique target value to each unique covariate vector
# used for the PCA (via `dat_clustered`).
cov_key_unique <- do.call(
  paste,
  c(dat_unique[, raster_covariates, drop = FALSE], sep = "|")
)
cov_key_clustered <- do.call(
  paste,
  c(dat_clustered[, raster_covariates, drop = FALSE], sep = "|")
)
y_unique <- as.numeric(
  dat_clustered$mean_carbon_density[match(cov_key_unique, cov_key_clustered)]
)
stopifnot(length(y_unique) == nrow(pca$x))

cor_with_target <- vapply(
  seq_len(ncol(pca$x)),
  function(j) cor(pca$x[, j], y_unique, use = "complete.obs"),
  numeric(1)
)

pc_names <- colnames(pca$x)
if (is.null(pc_names)) pc_names <- paste0("PC", seq_len(ncol(pca$x)))

sign_color <- ifelse(cor_with_target >= 0, "Positive", "Negative")
df_pc_cor <- data.frame(
  PC = factor(pc_names, levels = pc_names),
  correlation = as.numeric(cor_with_target),
  sign = sign_color
)

p <- ggplot(df_pc_cor, aes(x = PC, y = correlation, color = sign)) +
  geom_point(size = 3) +
  scale_color_manual(
    name = "Correlation Sign",
    values = c("Positive" = "blue", "Negative" = "red")
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(
    x = "Principal Component",
    y = "Correlation with Target"
  ) +
  theme_minimal()
print(p)


cat("\nDone.\n")