#!/usr/bin/env Rscript

# Plot approximate outlines of the broad Regions used for GPR analyses,
# with labels placed inside each region.
#
# Usage:
#   source("modelling/plot_region_outlines.R")
#   p <- plot_region_outlines()
#   ggsave("figures/interpolation/region_outlines.png", p, width = 10, height = 6)

rm(list = ls())
setwd(here::here())

source("modelling/R/helpers.R")
source("modelling/R/plot_config.R")
load_packages(c("here", "ggplot2", "dplyr", "maps"))

world <- map_data("world")

# load meow ecoregions shapefile
ecos <- sf::st_read("data/MEOW/meow_ecos.shp")
# ecos <- sf::st_read("data/ICES_ecoregions/ICES_ecoregions_20171207_erase_ESRI.shp")
unique(ecos$ECOREGION)


# create ecoregion clusters
med_sea <- c("Western Mediterranean", "Alboran Sea", "Levantine Sea", "Ionian Sea", "Aegean Sea")
north_european_atlantic <- c("Celtic Seas", "North Sea")
south_european_atlantic <- c("South European Atlantic Shelf")

# Group relevant MEOW ecoregion shapes into larger target regions and output as a dataframe

# Define lists of source eco-regions for grouping
med_sea <- c("Western Mediterranean", "Alboran Sea", "Levantine Sea", "Ionian Sea", "Aegean Sea")
north_european_atlantic <- c("Celtic Seas", "North Sea")
south_european_atlantic <- c("South European Atlantic Shelf")
baltic_sea <- c("Baltic Sea")
black_sea <- c("Black Sea")

# Create a data frame for grouped ecoregions
ecoregion_groups <- tibble::tribble(
  ~region,                     ~ecoregion_list,
  "Mediterranean Sea",         med_sea,
  "North European Atlantic",   north_european_atlantic,
  "South European Atlantic",   south_european_atlantic,
  "Baltic Sea",                baltic_sea,
  "Black Sea",                 black_sea
)

# For each grouped region, union the relevant MEOW polygons
region_shapes <- ecoregion_groups %>%
  rowwise() %>%
  mutate(
    geometry = sf::st_union(
      ecos %>%
        filter(ECOREGION %in% ecoregion_list) %>%
        sf::st_geometry()
    )
  ) %>%
  ungroup() %>%
  dplyr::select(region, geometry)


library(sf)
# Cast to sf data frame
region_shapes <- sf::st_sf(region_shapes)
world <- ggplot2::map_data("world")

# Convert to sf object
region_shapes <- sf::st_as_sf(region_shapes)

bbox <- st_bbox(region_shapes)
xmin <- bbox$xmin
xmax <- bbox$xmax
ymin <- bbox$ymin
ymax <- bbox$ymax

region_labels <- region_shapes %>%
  mutate(label_point = sf::st_point_on_surface(geometry))

# put a newline before Atlantic for each region label
region_labels <- region_labels %>%
  mutate(
    region = ifelse(region == "North European Atlantic", "North European\nAtlantic", region),
    region = ifelse(region == "South European Atlantic", "South European\nAtlantic", region)
  )



ggplot() +
  geom_polygon(
    data = world,
    aes(x = long, y = lat, group = group),
    fill = "#eeeeee",
    colour = "#a5a5a5"
  ) +
  geom_sf(
    data = region_shapes,
    fill = REGION_COLOURS[region_shapes$region],
    colour = "black",
    alpha = 0.5,
    linewidth = 0.2
  ) +
  geom_sf_text(
    data = region_labels,
    aes(label = region, geometry = label_point),
    size = 4,
    fontface = "bold",
    colour = "black"
  ) +
  coord_sf(
    xlim = c(xmin, xmax),
    ylim = c(ymin, ymax),
    expand = FALSE
  ) +
  theme_minimal() +
  theme(legend.position = "none") + labs(x = "Longitude", y = "Latitude")
# save
ggsave("figures/supplement/region_shapes.png", width = 10, height = 10)



# #' Return a ggplot of region outlines and labels.
# #'
# #' The bounding boxes are chosen to be consistent with those used in
# #' assign_region_from_latlon.R and can be refined later.
# plot_region_outlines <- function() {
#   regions_df <- data.frame(
#     region = c(
#       "Mediterranean Sea",
#       "Baltic Sea",
#       "Black Sea",
#       "North European Atlantic",
#       "South European Atlantic"
#     ),
#     xmin = c(-6, 10, 27, -20, -20),
#     xmax = c(36, 30, 42, 10, 10),
#     ymin = c(30, 53, 40, 45, 30),
#     ymax = c(46, 66, 47, 66, 45),
#     stringsAsFactors = FALSE
#   )

#   regions_df <- regions_df %>%
#     mutate(
#       xmid = (xmin + xmax) / 2,
#       ymid = (ymin + ymax) / 2
#     )

#   p <- ggplot() +
#     geom_polygon(
#       data = world,
#       aes(x = long, y = lat, group = group),
#       fill = "grey95",
#       color = "grey80",
#       linewidth = 0.2
#     ) +
#     geom_rect(
#       data = regions_df,
#       aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
#       fill = NA,
#       color = REGION_COLOURS[regions_df$region],
#       linewidth = 0.6
#     ) +
#     geom_text(
#       data = regions_df,
#       aes(x = xmid, y = ymid, label = region),
#       size = 3.5,
#       fontface = "bold",
#       color = "black"
#     ) +
#     coord_quickmap(
#       xlim = c(-25, 45),
#       ylim = c(25, 70),
#       expand = FALSE
#     ) +
#     theme_minimal() +
#     theme(
#       panel.background = element_rect(fill = "white", color = NA),
#       panel.grid = element_line(color = "grey90"),
#       axis.title = element_text(size = 10),
#       axis.text = element_text(size = 8)
#     ) +
#     labs(
#       x = "Longitude",
#       y = "Latitude",
#       title = "Approximate Region Outlines for GPR Analyses"
#     )

#   p
# }

# if (interactive()) {
#   print(plot_region_outlines())
# }

