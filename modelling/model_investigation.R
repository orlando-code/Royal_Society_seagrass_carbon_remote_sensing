### Investigate effect of different cross-validation methods on model performance

## TODO:
# move functions to separate file ,/
# hyperparameter tuning for all models
# compare effect of different target parameterisations (e.g. mean, sum, median) on model performance
# investigate two-model approach: one model for aggregate, the second capturing depth effects
# covariate pruning/feature importance investigation
# investigate GAMs

## ================================ SETUP ================================
# refresh and set working directory
rm(list = ls())
setwd(here::here())

# source helper functions
source("modelling/helpers.R")

load_packages(c("here", "mgcv", "tidyverse", "ggplot2", "visreg", "randomForest", "corrplot", "blockCV", "GGally", "cowplot", "gridExtra"))

## get data
For_modeling_df_shallow_carbon_density <- read_rds("data/For_modeling_df_shallow_carbon_density.rds")
dat <- For_modeling_df_shallow_carbon_density


## ================================ INVESTIGATE SPATIAL DISTRIBUTION OF DATA ================================
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


## ================================ INVESTIGATE CARBON DENSITY VS DEPTH RELATIONSHIP ================================
results <- fit_exponential_carbon_decay(dat)
nls_model <- results$nls_model
dat_clean <- results$dat_clean
pred_df <- results$pred_df
residuals_sd <- results$residuals_sd

## plot with fitted curve and confidence interval, with a line for each species
species_list <- unique(dat_clean$seagrass_species)
species_colors <- scales::hue_pal()(length(species_list))
names(species_colors) <- species_list
fits_by_species <- map(species_list, ~fit_species(dat_clean, .x))
# extract prediction dataframe and combine into a single dataframe
preds_df <- fits_by_species %>%
  compact() %>%                # remove NULLs
  map("pred_df") %>%
  bind_rows()

dat_clean$seagrass_species <- factor(dat_clean$seagrass_species, levels = names(species_colors))
# create fit labels for legend
fit_labels <- map2_chr(
  fits_by_species,
  species_list,
  ~{
    if (is.null(.x)) return(.y)
    a_val <- round(coef(.x$fit)["a"], 3)
    b_val <- round(coef(.x$fit)["b"], 4)
    paste0(.y, ": ", a_val, " * exp(", b_val, " * depth)")
  }
)
names(fit_labels) <- species_list

## plot
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
    fill = FALSE,
    scale="none"
  ) +
  labs(
    y = "Depth (cm)",
    x = expression("Carbon density (gC cm"^{-3}*")"),
    title = "Carbon density vs Depth by Seagrass Species"
  ) +
  # scale_x_continuous(limits = c(0, 0.1)) # outliers obscure the fit +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.key.width = unit(2.5, "cm")  # widen color swatches if needed
  )
ggsave("figures/carbon_density_vs_depth_by_species.png", plot = carbon_density_plot, width = 15, height = 6, dpi = 300)


## ================================ HOW DOES THE DEPTH-CARBON DENSITY RELATIONSHIP VARY BY SPECIES OR REGION? ================================
# approach: fit models with and without species/region interactions, compare statistically, and visualize

## explore metadata
cat("seagrass_species levels:", length(unique(dat_clean$seagrass_species)), "\n")
print(table(dat_clean$seagrass_species))
cat("\nregion levels:", length(unique(dat_clean$Region)), "\n")
print(table(dat_clean$Region))

## test each grouping variable: compares AIC for exponential decay models with and without species/region interactions
results_species <- compare_depth_relationships(dat_clean, "seagrass_species", "species")
results_region <- compare_depth_relationships(dat_clean, "Region", "region")

## plot
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
  

## ================================ INVESTIGATE RELATIONSHIPS BETWEEN ENVIRONMENTAL VARIABLES AND CARBON DENSITY ================================
env_vars <- c(
  "KD_closest", "RRS443_closest", "wave_height_VHM0_p95_m_closest", 
  "po4_mean_1.5m_mmol_m3_closest", "pH_mean_1.5m_closest", "bottomT_p95_C_closest", 
  "vo_p90_1.5m_m_s_closest", "uo_mean_1.5m_m_s_closest", "Surf_fgco2_p95_molC_m2_yr_closest"
)

plot_env_pairs(
  data = dat_clean,
  env_vars = env_vars,
  color_var = "carbon_density_g_c_cm3",
  output_file = "figures/environmental_variables_pairs.png",
  width = 15,
  height = 15,
  dpi = 300
)
