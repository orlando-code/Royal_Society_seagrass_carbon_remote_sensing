# GPR spatial prediction maps (shared by gpr_predictions.R and run_final_gpr.R)
#
# Creates mean and standard-error maps using geom_raster (dense, continuous-looking)
# with land mask on top and optional observation points.
#
# Usage:
#   source("modelling/plots/plot_gpr_spatial_maps.R")
#   plots <- plot_gpr_spatial_maps(grid, lon_range, lat_range, obs_data = dat, ...)

#' Plot GPR spatial predictions as raster maps.
#'
#' @param grid Data frame with longitude, latitude, gpr_mean, gpr_se
#' @param lon_range Length-2 vector c(min, max) for x-axis
#' @param lat_range Length-2 vector c(min, max) for y-axis
#' @param obs_data Optional data frame with longitude, latitude and value column (for observation points)
#' @param value_col Name of value column in obs_data (default: median_carbon_density_100cm)
#' @param mean_title Title for mean map
#' @param se_title Title for SE map
#' @param mean_name Legend label for mean
#' @param se_name Legend label for SE
#' @param width Plot width (default 8)
#' @param height Plot height (default 6)
#' @param dpi Resolution (default 150)
#' @return List with p_mean, p_se (ggplot objects)
plot_gpr_spatial_maps <- function(grid,
                                  lon_range,
                                  lat_range,
                                  obs_data = NULL,
                                  value_col = "median_carbon_density_100cm",
                                  mean_title = "GPR Mean Predictions",
                                  se_title = "GPR Prediction Uncertainty (Standard Error)",
                                  mean_name = "Predicted median\ncarbon density",
                                  se_name = "Standard error on\npredicted carbon density",
                                  width = 8,
                                  height = 6,
                                  dpi = 150) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plot_gpr_spatial_maps.")
  }
  if (!requireNamespace("maps", quietly = TRUE)) {
    stop("maps package is required for world coastlines.")
  }
  if (nrow(grid) == 0L || !all(c("longitude", "latitude", "gpr_mean", "gpr_se") %in% names(grid))) {
    return(list(p_mean = NULL, p_se = NULL))
  }

  world <- ggplot2::map_data("world")

  p_mean <- ggplot2::ggplot() +
    ggplot2::geom_raster(
      data = grid,
      ggplot2::aes(x = longitude, y = latitude, fill = gpr_mean)
    ) +
    ggplot2::scale_fill_viridis_c(option = "turbo", name = mean_name, na.value = "transparent") +
    ggplot2::geom_polygon(
      data = world,
      ggplot2::aes(x = long, y = lat, group = group),
      fill = "#eeeeee",
      colour = "#a5a5a5"
    )
  if (!is.null(obs_data) && nrow(obs_data) > 0L && all(c("longitude", "latitude", value_col) %in% names(obs_data))) {
    p_mean <- p_mean +
      ggplot2::geom_point(
        data = obs_data,
        ggplot2::aes(x = longitude, y = latitude, fill = .data[[value_col]]),
        shape = 21,
        size = 1.3,
        colour = "white",
        stroke = 0.5,
        show.legend = FALSE
      )
  }
  p_mean <- p_mean +
    ggplot2::coord_sf(xlim = lon_range, ylim = lat_range) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom", legend.key.width = ggplot2::unit(2, "cm")) +
    ggplot2::labs(x = "Longitude", y = "Latitude", title = mean_title)

  p_se <- ggplot2::ggplot() +
    ggplot2::geom_raster(
      data = grid,
      ggplot2::aes(x = longitude, y = latitude, fill = gpr_se)
    ) +
    ggplot2::scale_fill_viridis_c(option = "turbo", name = se_name, na.value = "transparent") +
    ggplot2::geom_polygon(
      data = world,
      ggplot2::aes(x = long, y = lat, group = group),
      fill = "#eeeeee",
      colour = "#a5a5a5"
    )
  if (!is.null(obs_data) && nrow(obs_data) > 0L && all(c("longitude", "latitude") %in% names(obs_data))) {
    p_se <- p_se +
      ggplot2::geom_point(
        data = obs_data,
        ggplot2::aes(x = longitude, y = latitude),
        shape = 21,
        size = 1.5,
        fill = "white",
        colour = "black",
        stroke = 1.2,
        show.legend = FALSE
      )
  }
  p_se <- p_se +
    ggplot2::coord_sf(xlim = lon_range, ylim = lat_range) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom", legend.key.width = ggplot2::unit(2, "cm")) +
    ggplot2::labs(x = "Longitude", y = "Latitude", title = se_title)

  list(p_mean = p_mean, p_se = p_se)
}
