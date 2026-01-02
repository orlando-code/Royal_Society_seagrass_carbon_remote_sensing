### Investigate effect of different cross-validation methods on model performance
# refresh working directory
rm(list = ls())
# Set working directory
# setwd(here::here())
setwd("/Users/rt582/Library/CloudStorage/OneDrive-UniversityofCambridge/cambridge/phd/Paper_Conferences/seagrass/Fw_ Following up on modeling conversation at the Royal Society Meeting/SG_Cstock_modeling")

# source helper functions
source("helpers.R")

load_packages(c("here", "mgcv", "tidyverse", "ggplot2", "visreg", "randomForest", "corrplot", "blockCV", "GGally"))

# TODO: 
# move functions to separate file ,/
# hyperparameter tuning for all models
# compare effect of different target parameterisations (e.g. mean, sum, median) on model performance
# investigate two-model approach: one model for aggregate, the second capturing depth effects





# DATA VISUALISATION
For_modeling_df_shallow_carbon_density <- read_rds("data/For_modeling_df_shallow_carbon_density.rds")
dat <- For_modeling_df_shallow_carbon_density

# plot the average carbon density as a function of depth. include shaded region as plus or minus one standard deviation.
# fit exponential decay model separately to avoid stat_smooth issues
fit_exponential_carbon_decay <- function(dat) {
  # Clean data
  dat_clean <- dat %>%
    filter(!is.na(carbon_density_g_c_cm3) & !is.na(sediment_mean_depth_cm) &
           carbon_density_g_c_cm3 > 0 & sediment_mean_depth_cm >= 0)
  
  if (nrow(dat_clean) < 3) {
    warning("Insufficient data points after cleaning for NLS fitting.")
    return(NULL)
  }
  
  # Fit nls model: carbon_density = a * exp(b * depth)
  nls_model <- tryCatch({
    nls(carbon_density_g_c_cm3 ~ a * exp(b * sediment_mean_depth_cm),
        data = dat_clean,
        start = list(a = max(dat_clean$carbon_density_g_c_cm3, na.rm = TRUE), 
                     b = -0.01))
  }, error = function(e) {
    cat("nls fitting failed, trying alternative starting values\n")
    nls(carbon_density_g_c_cm3 ~ a * exp(b * sediment_mean_depth_cm),
        data = dat_clean,
        start = list(a = mean(dat_clean$carbon_density_g_c_cm3, na.rm = TRUE), 
                     b = -0.001),
        control = nls.control(maxiter = 500, warnOnly = TRUE))
  })
  
  # Create prediction data frame
  depth_seq <- seq(min(dat_clean$sediment_mean_depth_cm, na.rm = TRUE),
                   max(dat_clean$sediment_mean_depth_cm, na.rm = TRUE),
                   length.out = 200)
  pred_df <- data.frame(sediment_mean_depth_cm = depth_seq)
  pred_df$carbon_density_fit <- predict(nls_model, newdata = pred_df)
  
  # Calculate approximate confidence intervals based on residual standard error
  residuals_sd <- sd(residuals(nls_model), na.rm = TRUE)
  pred_df$carbon_density_se <- residuals_sd
  pred_df$carbon_density_lower <- pred_df$carbon_density_fit - 1.96 * pred_df$carbon_density_se
  pred_df$carbon_density_upper <- pred_df$carbon_density_fit + 1.96 * pred_df$carbon_density_se

  return(list(
    nls_model = nls_model,
    dat_clean = dat_clean,
    pred_df = pred_df,
    residuals_sd = residuals_sd
  ))
}

results <- fit_exponential_carbon_decay(dat)
nls_model <- results$nls_model
dat_clean <- results$dat_clean
pred_df <- results$pred_df
residuals_sd <- results$residuals_sd

# plot with fitted curve and confidence interval, with a line for each species
# points are colored by seagrass species
# fit a line to each species and display its formula in the legend

# Fit model for each species and generate predictions
library(dplyr)
library(purrr)
library(tidyr)

species_list <- unique(dat_clean$seagrass_species)
species_colors <- scales::hue_pal()(length(species_list))
names(species_colors) <- species_list

# Fit model and generate prediction data frame per species
fit_species <- function(spec) {
  dat_sub <- dat_clean %>% filter(seagrass_species == spec)
  if (nrow(dat_sub) < 3) return(NULL)
  fit <- tryCatch({
    nls(carbon_density_g_c_cm3 ~ a * exp(b * sediment_mean_depth_cm),
        data = dat_sub,
        start = list(
          a = max(dat_sub$carbon_density_g_c_cm3, na.rm = TRUE),
          b = -0.01
        ))
  }, error = function(e) {
    nls(carbon_density_g_c_cm3 ~ a * exp(b * sediment_mean_depth_cm),
        data = dat_sub,
        start = list(
          a = mean(dat_sub$carbon_density_g_c_cm3, na.rm = TRUE),
          b = -0.001
        ),
        control = nls.control(maxiter = 500, warnOnly = TRUE))
  })
  # predictions
  depth_seq <- seq(
    min(dat_sub$sediment_mean_depth_cm, na.rm = TRUE),
    max(dat_sub$sediment_mean_depth_cm, na.rm = TRUE),
    length.out = 200
  )
  pred_df <- data.frame(
    seagrass_species = spec,
    sediment_mean_depth_cm = depth_seq
  )
  pred_df$carbon_density_fit <- predict(fit, newdata = pred_df)
  # approx CI
  residuals_sd <- sd(residuals(fit), na.rm = TRUE)
  pred_df$carbon_density_se <- residuals_sd
  pred_df$carbon_density_lower <- pred_df$carbon_density_fit - 1.96 * residuals_sd
  pred_df$carbon_density_upper <- pred_df$carbon_density_fit + 1.96 * residuals_sd
  list(fit = fit, pred_df = pred_df)
}

species_results <- set_names(species_list) %>%
  map(fit_species)

# Get all non-null fits
species_results <- species_results[!map_lgl(species_results, is.null)]

# Extract predictions, and coefficient formulas for each
preds_df <- map_dfr(species_results, "pred_df")
fits_by_species <- map(species_results, "fit")

# Prepare observed data and assign color by species
dat_clean$seagrass_species <- factor(dat_clean$seagrass_species, levels = names(species_colors))

# Generate fit formulas for legend
fit_labels <- map2_chr(
  fits_by_species,
  names(fits_by_species),
  ~{
    if (is.null(.x)) return(.y)
    a_val <- round(coef(.x)[["a"]], 3)
    b_val <- round(coef(.x)[["b"]], 4)
    paste0(.y, ": ", a_val, " * exp(", b_val, " * depth)")
  }
)
names(fit_labels) <- names(fits_by_species)

# Combine for ggplot
# Make the legend spread over three rows & increase width to prevent cutoff
carbon_density_plot <- ggplot() +
  geom_point(
    data = dat_clean,
    aes(
      x = carbon_density_g_c_cm3,
      y = sediment_mean_depth_cm,
      color = seagrass_species
    ),
    alpha = 0.5,
    show.legend = TRUE
  ) +
  geom_ribbon(
    data = preds_df,
    aes(
      y = sediment_mean_depth_cm,
      xmin = carbon_density_lower,
      xmax = carbon_density_upper,
      fill = seagrass_species
    ),
    alpha = 0.18,
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
  geom_line(
    data = preds_df,
    aes(
      x = carbon_density_fit,
      y = sediment_mean_depth_cm,
      color = seagrass_species
    ),
    linewidth = 1,
    inherit.aes = FALSE,
    show.legend = TRUE
  ) +
  scale_y_reverse() +
  scale_color_manual(
    name = "Seagrass species fits",
    values = species_colors,
    labels = fit_labels
  ) +
  scale_fill_manual(
    name = "Seagrass species",
    values = species_colors
  ) +
  guides(
    color = guide_legend(
      order = 1,
      override.aes = list(alpha = 1, linewidth = 2),
      nrow = 3 # force the legend to span three rows
    ),
    fill = FALSE
  ) +
  labs(
    y = "Depth (cm)",
    x = expression("Carbon density (gC cm"^{-3}*")"),
    title = "Carbon density vs Depth by Seagrass Species"
  ) +
  # scale_x_continuous(limits = c(0, 0.1)) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.key.width = unit(2.5, "cm")  # widen color swatches if needed
  )

# Save with increased width to further ensure legend is not cut off
ggsave("carbon_density_vs_depth_by_species.png", plot = carbon_density_plot, width = 15, height = 6, dpi = 300)


# systematic investigation: does the depth-carbon density relationship vary by species or location?
# approach: fit models with and without interactions, compare statistically, visualize

# check available grouping variables and their levels
cat("=== data structure for grouping variables ===\n")
cat("seagrass_species levels:", length(unique(dat_clean$seagrass_species)), "\n")
print(table(dat_clean$seagrass_species))
cat("\nregion levels:", length(unique(dat_clean$Region)), "\n")
print(table(dat_clean$Region))




# test each grouping variable
results_species <- compare_depth_relationships("seagrass_species", "species")
results_region <- compare_depth_relationships("Region", "region")

# create plots
cat("\n=== creating plots ===\n")
plot_species <- plot_depth_by_group(results_species, "species", xlim = c(0, 0.1)) # some outliers in 'unspecified' stretch axis and obscure the fit
plot_region <- plot_depth_by_group(results_region, "region", xlim = c(0, 0.1)) # similarly, outliers from Black Sea stretch axis and obscure the fit

# summary table of coefficients by group
cat("\n=== summary: exponential decay coefficients by group ===\n")
if (!is.null(results_species)) {
  cat("\nspecies:\n")
  print(results_species$group_coefs %>% 
        arrange(desc(a)) %>% 
        mutate(a = round(a, 4), b = round(b, 6)))
}
if (!is.null(results_region)) {
  cat("\nregion:\n")
  print(results_region$group_coefs %>% 
        arrange(desc(a)) %>% 
        mutate(a = round(a, 4), b = round(b, 6)))
}
  
# prepare data subset with environmental variables and carbon density
# carbon density needs to be in the data frame for the custom function to access it
env_vars <- c(
  "KD_closest", "RRS443_closest", "wave_height_VHM0_p95_m_closest", 
  "po4_mean_1.5m_mmol_m3_closest", "pH_mean_1.5m_closest", "bottomT_p95_C_closest", 
  "vo_p90_1.5m_m_s_closest", "uo_mean_1.5m_m_s_closest", "Surf_fgco2_p95_molC_m2_yr_closest"
)

dat_pairs <- dat_clean %>%
  dplyr::select(dplyr::all_of(c(env_vars, "carbon_density_g_c_cm3")))

# custom function for scatter plots with continuous color using turbo colormap
# turbo colormap is from viridisLite package (dependency of ggplot2)
scatter_with_color <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_point(aes(color = carbon_density_g_c_cm3), alpha = 0.5, ...) +
    scale_color_gradientn(
      colors = viridisLite::turbo(256),
      name = "Carbon density\n(gC cm^-3)"
    )
}

# plot pairplots of all environmental variables
ggpairs(
  dat_pairs, 
  columns = env_vars,
  upper = list(continuous = scatter_with_color),
  diag = list(continuous = "densityDiag")
) + 
  theme_minimal()
