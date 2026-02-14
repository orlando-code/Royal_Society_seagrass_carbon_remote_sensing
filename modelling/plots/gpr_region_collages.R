#!/usr/bin/env Rscript

# Helpers to create region-based collages of GPR predictions and SE.
#
# Usage (from an interactive session after running gpr_predictions.R):
#   source("modelling/gpr_region_collages.R")
#   collage <- create_gpr_region_collage(gpr_result$prediction_grid)
#   ggsave("figures/interpolation/gpr_region_collage.png", collage, width = 14, height = 10)
#

# rm(list = ls())
setwd(here::here())

source("modelling/R/helpers.R")
source("modelling/R/assign_region_from_latlon.R")
load_packages(c("here", "tidyverse", "ggplot2", "patchwork", "maps"))

world <- map_data("world")

#' Create a collage of regional GPR mean/SE maps.
#'
#' @param pred_grid Prediction grid data frame with columns longitude, latitude,
#'        gpr_mean, gpr_se, and region (inferred or observed).
#' @param regions Optional character vector of region names to include. If NULL,
#'        uses all non-NA unique regions in pred_grid.
#' @param sample_every Optional integer to thin points for plotting (e.g. 1 = all,
#'        5 = every 5th point). Defaults to 1.
#' @return A patchwork object combining per-region mean/SE panels.
create_gpr_region_collage <- function(pred_grid,
                                      regions = NULL,
                                      sample_every = 1L) {
  if (exists("show_titles", envir = .GlobalEnv)) show_titles <- get("show_titles", envir = .GlobalEnv) else show_titles <- TRUE
  dat <- pred_grid
  if (!"region" %in% names(dat)) {
    if (all(c("longitude", "latitude") %in% names(dat))) {
      dat <- assign_region_from_latlon(dat)
    } else {
      stop("create_gpr_region_collage: pred_grid must have 'region' or lon/lat to infer it.")
    }
  }
  dat <- dat[!is.na(dat$region) &
    !is.na(dat$gpr_mean) &
    !is.na(dat$gpr_se) &
    !is.na(dat$longitude) &
    !is.na(dat$latitude), , drop = FALSE]
  if (nrow(dat) == 0) stop("create_gpr_region_collage: no valid prediction points with region and GPR outputs.")

  if (is.null(regions)) {
    regions <- sort(unique(dat$region))
  }
  if (length(regions) == 0) stop("create_gpr_region_collage: no regions to plot.")

  region_plots <- list()
  for (r in regions) {
    df <- dat[dat$region == r, , drop = FALSE]
    if (nrow(df) == 0) next

    # Optional thinning to keep plots light
    if (!is.null(sample_every) && sample_every > 1L) {
      idx <- seq(1L, nrow(df), by = sample_every)
      df <- df[idx, , drop = FALSE]
    }

    lon_range <- range(df$longitude, na.rm = TRUE)
    lat_range <- range(df$latitude, na.rm = TRUE)
    lon_buffer <- max(0.5, 0.1 * diff(lon_range))
    lat_buffer <- max(0.5, 0.1 * diff(lat_range))

    p_mean <- ggplot() +
      geom_polygon(
        data = world,
        aes(x = long, y = lat, group = group),
        fill = "grey95",
        color = "grey60",
        linewidth = 0.2
      ) +
      geom_point(
        data = df,
        aes(x = longitude, y = latitude, color = gpr_mean),
        size = 0.2,
        alpha = 0.9
      ) +
      scale_color_gradientn(
        colors = viridisLite::turbo(256),
        name = "GPR mean"
      ) +
      coord_cartesian(
        xlim = c(lon_range[1] - lon_buffer, lon_range[2] + lon_buffer),
        ylim = c(lat_range[1] - lat_buffer, lat_range[2] + lat_buffer),
        expand = FALSE
      ) +
      theme_minimal() +
      theme(
        legend.position = "bottom",
        legend.key.width = unit(2, "cm"),
        panel.background = element_rect(fill = "white", color = NA),
        panel.grid = element_line(color = "grey90")
      ) +
      labs(
        x = "Longitude",
        y = "Latitude",
        title = if (show_titles) paste("Mean:", r) else NULL
      )

    p_se <- ggplot() +
      geom_polygon(
        data = world,
        aes(x = long, y = lat, group = group),
        fill = "grey95",
        color = "grey60",
        linewidth = 0.2
      ) +
      geom_point(
        data = df,
        aes(x = longitude, y = latitude, color = gpr_se),
        size = 0.2,
        alpha = 0.9
      ) +
      scale_color_gradientn(
        colors = viridisLite::viridis(256),
        name = "SE"
      ) +
      coord_cartesian(
        xlim = c(lon_range[1] - lon_buffer, lon_range[2] + lon_buffer),
        ylim = c(lat_range[1] - lat_buffer, lat_range[2] + lat_buffer),
        expand = FALSE
      ) +
      theme_minimal() +
      theme(
        legend.position = "bottom",
        legend.key.width = unit(2, "cm"),
        panel.background = element_rect(fill = "white", color = NA),
        panel.grid = element_line(color = "grey90")
      ) +
      labs(
        x = "Longitude",
        y = "Latitude",
        title = if (show_titles) paste("SE:", r) else NULL
      )

    region_plots[[r]] <- p_mean + p_se + plot_layout(ncol = 2)
  }

  if (length(region_plots) == 0) stop("create_gpr_region_collage: no regions had data to plot.")

  # Combine all region panels vertically by default
  collage <- wrap_plots(region_plots, ncol = 1)
  collage
}
