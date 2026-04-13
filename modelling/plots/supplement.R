#!/usr/bin/env Rscript
# =============================================================================
# Supplement figures for the paper. Run after data and (optionally) covariate
# pruning exist. All outputs save to output/supplement/.
# Optional env for Figure 10 (covariate NA map): NA_MAP_N_LON, NA_MAP_N_LAT (defaults 240×200).
#
# Usage: setwd(project_root); source("modelling/plots/supplement.R")
#        Or run after run_paper.R (globals set).
# =============================================================================

if (!exists("seagrass_init_repo", mode = "function", inherits = TRUE)) {
  init_path <- file.path("modelling", "R", "init_repo.R")
  if (!file.exists(init_path)) {
    ff <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
    if (!length(ff)) stop("Run from repo root or with: Rscript /path/to/this_script.R", call. = FALSE)
    script_path <- normalizePath(sub("^--file=", "", ff[[1]]), winslash = "/", mustWork = FALSE)
    init_path <- normalizePath(file.path(dirname(script_path), "..", "R", "init_repo.R"), winslash = "/", mustWork = FALSE)
  }
  if (!file.exists(init_path)) stop("Missing bootstrap helper: modelling/R/init_repo.R", call. = FALSE)
  sys.source(init_path, envir = .GlobalEnv)
}
project_root <- seagrass_init_repo(
  include_helpers = FALSE,
  require_core_inputs = FALSE,
  check_renv = FALSE
)
source("modelling/R/helpers.R")
source("modelling/plots/plot_config.R")
source("modelling/R/assign_region_from_latlon.R")
source("modelling/pipeline_config.R")
source("modelling/R/extract_covariates_from_rasters.R")
source("modelling/R/covariate_na_map.R")
load_packages(c("here", "ggplot2", "dplyr", "tidyr", "maps", "sf", "patchwork", "FNN"))

cfg <- get_pipeline_config()
apply_pipeline_defaults(
  cfg,
  c(
    "dpi", "target_var", "exclude_regions", "use_correlation_filter",
    "correlation_filter_threshold", "n_lons", "n_lats",
    "robust_fold_seed_list", "cv_regime_name", "model_list"
  ),
  envir = .GlobalEnv
)
if (!"dplyr" %in% loadedNamespaces()) suppressPackageStartupMessages(library(dplyr))
if (exists("dpi", envir = .GlobalEnv)) dpi <- get("dpi", envir = .GlobalEnv) else dpi <- 150

# Output directory (theme: theme_paper() from modelling/plots/plot_config.R, sourced above)
OUT_DIR <- "output/supplement"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Data (load once; used by multiple figures)
# -----------------------------------------------------------------------------
target_var <- get("target_var", envir = .GlobalEnv)
exclude_regions <- get("exclude_regions", envir = .GlobalEnv)
use_correlation_filter <- isTRUE(get("use_correlation_filter", envir = .GlobalEnv))
cor_threshold <- get("correlation_filter_threshold", envir = .GlobalEnv)

dat <- readr::read_rds("data/all_extracted_new.rds")
dat <- process_rs_covariates(dat)
if (!"region" %in% names(dat)) dat <- assign_region_from_latlon(dat)
if (length(exclude_regions) > 0L) {
  dat <- dat[is.na(dat$region) | !dat$region %in% exclude_regions, ]
}

# Predictor candidates (raster/env only, no coords/target/species/region)
all_pred_vars <- setdiff(
  colnames(dat),
  c("latitude", "longitude", "number_id_final_version", "seagrass_species", "region", target_var)
)
all_pred_vars <- all_pred_vars[vapply(dat[all_pred_vars], is.numeric, logical(1))]

n_lons <- as.integer(get("n_lons", envir = .GlobalEnv))
n_lats <- as.integer(get("n_lats", envir = .GlobalEnv))
lon_min <- min(dat$longitude, na.rm = TRUE)
lon_max <- max(dat$longitude, na.rm = TRUE)
lat_min <- min(dat$latitude, na.rm = TRUE)
lat_max <- max(dat$latitude, na.rm = TRUE)

# Use first pruned covariate as reference raster (or find one that exists)
reference_raster <- all_pred_vars[1]
cat("Using", reference_raster, "as reference raster for filtering ocean cells\n\n")

cat("Creating prediction grid from rasters (", length(all_pred_vars), " covariates)...\n", sep = "")
prediction_grid <- create_prediction_grid_from_rasters(
  lon_range = c(lon_min, lon_max),
  lat_range = c(lat_min, lat_max),
  n_lon = n_lons,
  n_lat = n_lats,
  covariates = all_pred_vars,
  reference_raster = reference_raster,
  cache_path = sprintf("output/cache/prediction_grid_cache_lons%s_lats%s_all_pred_vars.rds", paste(n_lons, collapse = "_"), paste(n_lats, collapse = "_")),
  use_cache = TRUE,
  show_progress = TRUE,
  fill_coastal = TRUE
)

# Correlation filter: kept vars and full cor-with-target for plotting
if (use_correlation_filter && length(all_pred_vars) >= 2L) {
  kept_vars <- prune_by_correlation(dat, all_pred_vars, target_var, cor_threshold = cor_threshold)
  removed_vars <- setdiff(all_pred_vars, kept_vars)
} else {
  kept_vars <- all_pred_vars
  removed_vars <- character(0)
}

# Correlations with target (for bar chart)
dat_complete <- dat[complete.cases(dat[, c(target_var, all_pred_vars), drop = FALSE]), ]
cor_with_target <- if (nrow(dat_complete) >= 3L && length(all_pred_vars) > 0L) {
  setNames(
    vapply(all_pred_vars, function(v) cor(dat_complete[[v]], dat_complete[[target_var]], use = "complete.obs"), numeric(1)),
    all_pred_vars
  )
} else {
  setNames(rep(NA_real_, length(all_pred_vars)), all_pred_vars)
}
cor_df <- data.frame(
  variable = names(cor_with_target),
  correlation = as.numeric(cor_with_target),
  kept = names(cor_with_target) %in% kept_vars,
  stringsAsFactors = FALSE
)
cor_df <- cor_df[order(-abs(cor_df$correlation)), ]
cor_df$variable <- factor(cor_df$variable, levels = rev(cor_df$variable))

world <- map_data("world")

# -----------------------------------------------------------------------------
# Figure 1: Region outlines (MEOW-based)
# -----------------------------------------------------------------------------
meow_shp <- "data/MEOW/meow_ecos.shp"
if (file.exists(meow_shp)) {
  ecos <- sf::st_read(meow_shp, quiet = TRUE)
  ecoregion_groups <- tibble::tribble(
    ~region,                     ~ecoregion_list,
    "Mediterranean Sea",         c("Western Mediterranean", "Alboran Sea", "Levantine Sea", "Ionian Sea", "Aegean Sea"),
    "North European Atlantic",  c("Celtic Seas", "North Sea"),
    "South European Atlantic",  c("South European Atlantic Shelf"),
    "Baltic Sea",                "Baltic Sea",
    "Black Sea",                 "Black Sea"
  )
  region_shapes <- ecoregion_groups %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      geometry = sf::st_union(
        sf::st_geometry(ecos[ecos$ECOREGION %in% ecoregion_list, ])
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(region, geometry) %>%
    sf::st_as_sf()
  bbox <- sf::st_bbox(region_shapes)
  region_labels <- region_shapes %>%
    dplyr::mutate(
      label_point = suppressWarnings(sf::st_point_on_surface(geometry)),
      region = ifelse(region == "North European Atlantic", "North European\nAtlantic", region),
      region = ifelse(region == "South European Atlantic", "South European\nAtlantic", region)
    )

  p_region <- ggplot() +
    geom_polygon(data = world, aes(x = long, y = lat, group = group), fill = "#eeeeee", colour = "#a5a5a5") +
    geom_sf(data = region_shapes, fill = REGION_COLOURS[region_shapes$region], colour = "black", alpha = 0.5, linewidth = 0.2) +
    geom_sf_text(data = region_labels, aes(label = region, geometry = label_point), size = 4, fontface = "bold", colour = "black") +
    coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]), ylim = c(bbox["ymin"], bbox["ymax"]), expand = FALSE) +
    theme_paper() +
    theme(legend.position = "none") +
    labs(x = "Longitude", y = "Latitude")
  ggsave(file.path(OUT_DIR, "region_shapes.png"), p_region, width = 10, height = 10, dpi = dpi)
  cat("Saved", file.path(OUT_DIR, "region_shapes.png"), "\n")
} else {
  warning("Skipping region_shapes figure; missing optional file: ", meow_shp)
}

# -----------------------------------------------------------------------------
# Figure 2: Target vs log-transformed target (two histograms side by side)
# -----------------------------------------------------------------------------
y_lab <- if (target_var %in% names(VAR_LABELS)) VAR_LABELS[target_var] else gsub("_", " ", target_var)
dat_hist <- dat %>%
  dplyr::filter(!is.na(.data[[target_var]]) & .data[[target_var]] > 0)
p_raw <- ggplot(dat_hist, aes(x = .data[[target_var]])) +
  geom_histogram(bins = 40, fill = "#3366CC", colour = "grey30", alpha = 0.85) +
  labs(x = y_lab, y = "Count", title = "Raw") +
  theme_paper()
p_log <- ggplot(dat_hist, aes(x = log(pmax(.data[[target_var]], .Machine$double.eps)))) +
  geom_histogram(bins = 40, fill = "#CC6633", colour = "grey30", alpha = 0.85) +
  labs(x = paste("Log(", y_lab, ")"), y = "Count", title = "Log-transformed") +
  theme_paper()
p_target_hist <- p_raw + p_log +
  patchwork::plot_annotation(title = "Distribution of target variable", theme = theme(plot.title = element_text(face = "bold", hjust = 0.5)))
ggsave(file.path(OUT_DIR, "target_and_log_target_histogram.png"), p_target_hist, width = 10, height = 5, dpi = dpi)
cat("Saved", file.path(OUT_DIR, "target_and_log_target_histogram.png"), "\n")

# -----------------------------------------------------------------------------
# Figure 3: Maps of points by (a) species, (b) region, (c) count per location
# -----------------------------------------------------------------------------
lon_lim <- range(dat$longitude, na.rm = TRUE)
lat_lim <- range(dat$latitude, na.rm = TRUE)
base_map <- geom_polygon(data = world, aes(x = long, y = lat, group = group), fill = "#eeeeee", colour = "#a5a5a5", linewidth = 0.2)

# (a) Species
if ("seagrass_species" %in% names(dat)) {
  dat_spp <- dat %>% dplyr::filter(!is.na(seagrass_species))
  spp_cols <- SPECIES_COLOURS[names(SPECIES_COLOURS) %in% unique(dat_spp$seagrass_species)]
  p_species <- ggplot(dat_spp, aes(x = longitude, y = latitude, colour = seagrass_species)) +
    base_map +
    geom_point(size = 2, alpha = 0.8, stroke = 1) +
    scale_colour_manual(values = spp_cols, name = "Species", na.value = "grey70") +
    coord_cartesian(xlim = lon_lim, ylim = lat_lim) +
    labs(x = "Longitude", y = "Latitude") +
    theme_paper()
  ggsave(file.path(OUT_DIR, "map_points_by_species.png"), p_species, width = 9, height = 6, dpi = dpi)
  cat("Saved", file.path(OUT_DIR, "map_points_by_species.png"), "\n")
}

# (b) Region 
dat_reg <- dat %>% dplyr::filter(!is.na(region))
p_region_map <- ggplot(dat_reg, aes(x = longitude, y = latitude, colour = region)) +
  base_map +
  geom_point(size = 2, alpha = 0.8, stroke = 1) +
  scale_colour_manual(values = REGION_COLOURS, name = "Region", na.value = "grey70") +
  coord_cartesian(xlim = lon_lim, ylim = lat_lim) +
  labs(x = "Longitude", y = "Latitude") +
  theme_paper()
ggsave(file.path(OUT_DIR, "map_points_by_region.png"), p_region_map, width = 9, height = 6, dpi = dpi)
cat("Saved", file.path(OUT_DIR, "map_points_by_region.png"), "\n")

# (c) Count per unique lat/lon (only locations that appear in the data; table() would create a full grid of zeros)
location_counts <- dat %>%
  dplyr::count(longitude, latitude, name = "n")
p_count <- ggplot(location_counts, aes(x = longitude, y = latitude, fill = n)) +
  base_map +
  geom_point(aes(size = n), shape = 21, colour = "black", stroke = 0.3, alpha = 0.85) +
  scale_fill_viridis_c(option = "plasma", name = "Samples", begin = 0.15, end = 0.95, direction = -1) +
  scale_size_continuous(
    range = c(1, 6), 
    name = "", 
    breaks = c(1, 5, 10, 14),
    labels = c("1", "5", "10", "14"),
    guide = guide_legend(nrow = 2)
  ) +
  coord_cartesian(xlim = lon_lim, ylim = lat_lim) +
  labs(x = "Longitude", y = "Latitude") +
  theme_paper() +
  theme(legend.key.width = unit(2, "cm"))
ggsave(file.path(OUT_DIR, "map_points_by_count.png"), p_count, width = 9, height = 6, dpi = dpi)
cat("Saved", file.path(OUT_DIR, "map_points_by_count.png"), "\n")

# -----------------------------------------------------------------------------
# Figure 4: Correlation matrix of raster variables (post-correlation removal)
# -----------------------------------------------------------------------------
if (length(kept_vars) >= 2L) {
  dat_cor <- dat_complete[, kept_vars, drop = FALSE]
  C <- cor(dat_cor, use = "complete.obs")
  # Upper triangle only for a clean matrix plot
  C_upper <- C
  C_upper[lower.tri(C_upper, diag = TRUE)] <- NA
  C_long <- as.data.frame(as.table(C_upper))
  C_df <- C_long[!is.na(C_long$Freq), ]
  names(C_df) <- c("Var1", "Var2", "value")
  C_df <- C_df %>%
    dplyr::mutate(Var1 = factor(Var1, levels = rev(unique(as.character(Var1)))),
           Var2 = factor(Var2, levels = unique(as.character(Var2))))
  p_cormat <- ggplot(C_df, aes(Var2, Var1, fill = value)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0, limits = c(-1, 1), name = "Correlation") +
    labs(x = NULL, y = NULL) +
    theme_paper() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
      axis.text.y = element_text(size = 7),
      legend.position = "right"
    ) +
    coord_fixed()
  ggsave(file.path(OUT_DIR, "correlation_matrix_post_filter.png"), p_cormat, width = 10, height = 9, dpi = dpi)
  cat("Saved", file.path(OUT_DIR, "correlation_matrix_post_filter.png"), "\n")
}

# -----------------------------------------------------------------------------
# Figure 5: Bar chart of correlation with target (kept vs removed by filter)
# -----------------------------------------------------------------------------
cor_df_plot <- cor_df %>%
  dplyr::mutate(
    label = ifelse(as.character(variable) %in% names(VAR_LABELS), VAR_LABELS[as.character(variable)], as.character(variable)),
    status = ifelse(kept, "Retained", "Removed by correlation filter")
  )
y_labels <- label_vars(levels(cor_df_plot$variable))
y_labels <- ifelse(nchar(y_labels) > 50, paste0(substr(y_labels, 1, 47), "..."), y_labels)
p_cor_bar <- ggplot(cor_df_plot, aes(x = correlation, y = variable, fill = status)) +
  geom_col(width = 0.75) +
  geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.3) +
  scale_fill_manual(values = c("Retained" = "#2166ac", "Removed by correlation filter" = "#bdbdbd"), name = NULL) +
  scale_y_discrete(labels = y_labels) +
  labs(x = "Correlation with target", y = NULL,) +
  theme_paper() +
  theme(
    axis.text.y = element_text(size = 8),
    legend.position = "bottom"
  )
ggsave(file.path(OUT_DIR, "correlation_with_target_bar.png"), p_cor_bar, width = 8, height = max(5, length(all_pred_vars) * 0.15), dpi = dpi)
cat("Saved", file.path(OUT_DIR, "correlation_with_target_bar.png"), "\n")

# -----------------------------------------------------------------------------
# Figure 6: Environmental similarity (raster map + histogram). N.B. this uses all the predictor variables: would be more informative to use only the pruned variables for each model.
# -----------------------------------------------------------------------------

cat("\n==============================\n")
# for each model, get the pruned variables and calculate/save the environmental similarity scores
robust_fold_seed_list <- get("robust_fold_seed_list", envir = .GlobalEnv)
cv_regime_name_supp <- get("cv_regime_name", envir = .GlobalEnv)
seeds_str_supp <- paste(robust_fold_seed_list, collapse = "-")
robust_pruned_path <- file.path(
  "output", cv_regime_name_supp, "covariate_selection", "robust_pixel_grouped",
  paste0("pruned_model_variables_shap_robust_pixel_grouped_seeds_", seeds_str_supp, ".csv")
)
# Fall back to non-robust SHAP vars if robust file doesn't exist
if (!file.exists(robust_pruned_path)) {
  robust_pruned_path <- file.path("output", cv_regime_name_supp, "covariate_selection", "pruned_model_variables_shap.csv")
}
model_list <- get("model_list", envir = .GlobalEnv)
ess_similarity_percentile <- as.numeric(get0("ess_similarity_percentile", envir = .GlobalEnv, ifnotfound = 0.10))
ess_similarity_percentile <- max(0, min(1, ess_similarity_percentile))
ess_percent_label <- paste0(round(ess_similarity_percentile * 100), "th")

for (model in model_list) {
  cat("\nComputing environmental similarity for", model, "\n\n")
  all_models_pruned_variables <- read.csv(robust_pruned_path, stringsAsFactors = FALSE)
  pruned_variables <- all_models_pruned_variables[all_models_pruned_variables$model == model, "variable"]
  similarity <- compute_environmental_similarity(
    dat, prediction_grid, pruned_variables,
    method = "euclidean", verbose = TRUE
  )
  core_similarity_scores <- if (!is.null(similarity$data_similarity_reference_scores)) {
    similarity$data_similarity_reference_scores
  } else {
    similarity$data_similarity_scores
  }
  data_similarity_score_cutoff <- quantile(
    core_similarity_scores,
    ess_similarity_percentile,
    na.rm = TRUE
  )
  outside_ess <- is.finite(similarity$similarity_scores) &
    (similarity$similarity_scores < data_similarity_score_cutoff)
  # Mask cells outside ESS threshold for concise downstream use.
  similarity_scores_masked <- similarity$similarity_scores
  similarity_scores_masked[outside_ess] <- NaN

  p <- plot_applicability_domain(
    prediction_grid,
    similarity_scores = similarity_scores_masked,
    flag_outside = NULL,
    # flag_outside = outside_ess,
    plot_data_points = TRUE,
    point_size=1.5,
    dat = dat,
    data_similarity_scores = similarity$data_similarity_scores,
    world = world,
    xlim = EUROPE_LON_RANGE,
    ylim = EUROPE_LAT_RANGE,
    use_raster = TRUE
  )
  ggsave(file.path(OUT_DIR, paste0("applicability_domain_", model, ".png")), p, width = 8, height = 6)
  cat("Saved", file.path(OUT_DIR, paste0("applicability_domain_", model, ".png")), "\n")


  # calculate the 5th percentile similarity score in data
  # get indices of gebco cells between -70 and 0 depth
  in_range_indices <- which(prediction_grid$gebco_2025_n61.0_s34.0_w_10.0_e35.0 >= -70 & prediction_grid$gebco_2025_n61.0_s34.0_w_10.0_e35.0 <= 0)

  # filter similarity scores to only include the in_range_indices
  in_range_similarity_scores <- similarity_scores_masked[in_range_indices]


  # plot histogram of data similarity scores
  p <- ggplot(data.frame(data_similarity_scores = core_similarity_scores), aes(x = data_similarity_scores)) +
    geom_histogram(binwidth = 0.01, fill = "steelblue", color = "black", alpha = 0.7) +
    geom_vline(
      xintercept = data_similarity_score_cutoff,
      linetype = "dashed",
      color = "red",
      linewidth = 1
    ) +
    annotate(
      "text",
      x = data_similarity_score_cutoff,
      y = Inf,
      label = paste0(ess_percent_label, " percentile of core data similarity scores"),
      color = "red",
      angle = 90,
      vjust = -0.5,
      hjust = 1.1,
      size = 4
    ) +
    labs(
      x = "Euclidean Distance Similarity Score",
      y = "Count",
    ) +
    theme_paper()
  ggsave(file.path(OUT_DIR, paste0("data_similarity_scores_histogram_", model, ".png")), p, width = 8, height = 6)
  cat("Saved", file.path(OUT_DIR, paste0("data_similarity_scores_histogram_", model, ".png")), "\n")

  # plot a histogram of the similarity scores
  p <- ggplot(data.frame(similarity_scores = in_range_similarity_scores), aes(x = similarity_scores)) +
    geom_histogram(binwidth = 0.01, fill = "steelblue", color = "black", alpha = 0.7) +
    geom_vline(xintercept = data_similarity_score_cutoff, linetype = "dashed", color = "red", linewidth = 1) +
    labs(x = "Similarity Score", y = "Count") +
      annotate(
      "text",
      x = data_similarity_score_cutoff,
      y = Inf,
      label = paste0(ess_percent_label, " percentile of core data similarity scores"),
      color = "red",
      angle = 90,
      vjust = -0.5,
      hjust = 1.1,
      size = 4
    ) +
    labs(
      x = "Euclidean Distance Similarity Score",
      y = "Count",
    ) +
    theme_paper()
  ggsave(file.path(OUT_DIR, paste0("similarity_scores_histogram_", model, ".png")), p, width = 8, height = 6)
  cat("Saved", file.path(OUT_DIR, paste0("similarity_scores_histogram_", model, ".png")), "\n")
}
# TODO: this could be used to mask prediction grid so that only cells with a high-enough similarity score are used for prediction

# -----------------------------------------------------------------------------
# Figure 7: Scatter plot and correlation of surface pCO2 (95th percentile) and surface pCO2 (mean)
# -----------------------------------------------------------------------------

p <- ggplot(dat, aes(x = spco2_monthly_mean_pa, y = surf_spco2_mean_uatm)) +
  geom_point() +
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  labs(x = expression("Surface pCO"[2]*" (mean, Pa)"), 
       y = expression("Surface pCO"[2]*" (mean, µatm)")) +
  annotate(
    "text",
    x = Inf, y = Inf, 
    label = paste0("Correlation: ", round(cor(dat$spco2_monthly_p95_pa, dat$surf_spco2_p95_uatm, use = "complete.obs"), 2)),
    hjust = 1, vjust = 1, color = "red"
  ) +
  theme_paper()

ggsave(file.path(OUT_DIR, "scatter_spco2_monthly_mean_pa_vs_surf_spco2_mean_uatm.png"), p, width = 5, height = 5, dpi = 300)
cat("Saved", file.path(OUT_DIR, "scatter_spco2_monthly_p95_pa_vs_surf_spco2_p95_uatm.png"), "\n")

# -----------------------------------------------------------------------------
# Figure 8 – distribution of median carbon density for non-unique environmental covariates
# -----------------------------------------------------------------------------

# assign unique id to each environmental covariate vector
dat_clustered <- dat %>%
  group_by(across(all_of(raster_covariates))) %>%
  summarise(
    mean_carbon_density = mean(median_carbon_density_100cm, na.rm = TRUE),
    sd_carbon_density = sd(median_carbon_density_100cm, na.rm = TRUE),
    n_observations = n(),
    .groups = "drop"
  )
dat$id <- match(
  do.call(paste, c(dat[, raster_covariates, drop = FALSE], sep = "|")),
  do.call(paste, c(dat_clustered[, raster_covariates, drop = FALSE], sep = "|"))
)

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
  scale_fill_viridis_c(option = "plasma", direction = -1, name = "# of observations") +
  theme_paper() +
  # set the x tick labels as 1 to len(filtered_dat$id)
  scale_x_discrete(labels = as.character(1:length(filtered_dat$id))) +
  labs(y = expression("Distribution of median carbon density (gC/cm"^3*")"), x = "Unique covariate set index")
ggsave(file.path(OUT_DIR, "boxplot_carbon_density_by_id.png"), p, width = 10, height = 5, dpi = dpi)
cat("Saved", file.path(OUT_DIR, "boxplot_carbon_density_by_id.png"), "\n")


# -----------------------------------------------------------------------------
# Figure 9 – stacked histogram of median carbon density, shaded by region
# -----------------------------------------------------------------------------

# assign region to dat via assign_region_from_latlon

plot_dat <- dat %>%
  dplyr::filter(!is.na(.data[[target_var]]), !is.na(region))

region_medians <- aggregate(plot_dat[[target_var]], by = list(region = plot_dat$region), FUN = median, na.rm = TRUE)
ordered_regions <- region_medians$region[order(region_medians$x)]
plot_dat$region <- factor(plot_dat$region, levels = ordered_regions)

region_cols <- REGION_COLOURS[names(REGION_COLOURS) %in% levels(plot_dat$region)]

p <- ggplot(plot_dat, aes(x = .data[[target_var]], fill = region)) +
  geom_histogram(position = "stack", binwidth = 0.01, boundary = 0, color = NA, alpha = 0.5) +
  scale_fill_manual(values = region_cols, name = "Region", drop = FALSE) +
  theme_paper() +
  labs(
    x = expression("Median carbon density (gC/cm"^3*")"),
    y = "Sample count"
  )
ggsave(file.path(OUT_DIR, "histogram_carbon_density_by_region.png"), p, width = 10, height = 5, dpi = dpi)
cat("Saved", file.path(OUT_DIR, "histogram_carbon_density_by_region.png"), "\n")


# -----------------------------------------------------------------------------
# Figure 10: Map of how many base raster covariates are NA per cell (Europe)
# Raw extraction (no nearest-neighbour fill); see build_covariate_na_count_grid().
# -----------------------------------------------------------------------------
na_map_lon_range <- EUROPE_LON_RANGE
na_map_lat_range <- EUROPE_LAT_RANGE
na_map_n_lon <- as.integer(Sys.getenv("NA_MAP_N_LON", unset = "240"))
na_map_n_lat <- as.integer(Sys.getenv("NA_MAP_N_LAT", unset = "200"))
na_map_cache_rds <- "output/cache/covariate_na_count_grid_europe.rds"
cov_na_grid <- build_covariate_na_count_grid(
  lon_range = na_map_lon_range,
  lat_range = na_map_lat_range,
  n_lon = na_map_n_lon,
  n_lat = na_map_n_lat,
  covariates = raster_covariates,
  use_closest = FALSE,
  method = "nearest",
  drop_all_na = FALSE,
  cache_path = na_map_cache_rds,
  use_cache = TRUE,
  show_progress = TRUE
)

cov_na_plot_df <- cov_na_grid |>
  dplyr::arrange(latitude, longitude) |>
  dplyr::mutate(
    n_na_ocean = dplyr::if_else(n_na >= n_covariates, NA_integer_, n_na)
  )

oc_na <- cov_na_plot_df$n_na < cov_na_plot_df$n_covariates
mx_na <- if (any(oc_na, na.rm = TRUE)) {
  as.integer(max(cov_na_plot_df$n_na[oc_na], na.rm = TRUE))
} else {
  0L
}
cov_na_plot_df$n_na_lab <- factor(cov_na_plot_df$n_na_ocean, levels = seq.int(0L, mx_na))

p_cov_na <- ggplot(cov_na_plot_df, aes(x = longitude, y = latitude, fill = n_na_ocean)) +
  geom_raster(interpolate = FALSE) +
  geom_polygon(
    data = world,
    aes(x = long, y = lat, group = group),
    inherit.aes = FALSE,
    fill = "#e8e8e8",
    colour = "#b0b0b0",
    linewidth = 0.15
  ) +
  coord_cartesian(xlim = na_map_lon_range, ylim = na_map_lat_range, expand = FALSE) +
  scale_fill_viridis_c(
    option = "plasma",
    name = "Number of NA values:",
    begin = 0.4,
    end = 0.95,
    direction = -1,
    na.value = "transparent"
  ) +
  geom_point(
    data = dat,
    aes(x = longitude, y = latitude),
    inherit.aes = FALSE,
    color = "black",
    size = 1.5,
    alpha = 0.5
  ) +
  labs(x = "Longitude", y = "Latitude") +
  theme_paper() +
  guides(
    fill = guide_colourbar(
      barheight = 0.5,
      barwidth = 18,
      title.position = "left",
      title.theme = element_text(margin = ggplot2::margin(r = 14))
    )
  )

ggsave(file.path(OUT_DIR, "covariate_na_count_europe.png"), p_cov_na, width = 10, height = 8, dpi = dpi)
cat("Saved", file.path(OUT_DIR, "covariate_na_count_europe.png"), "\n")

cat("\nSupplement figures complete. Outputs in", OUT_DIR, "\n")


