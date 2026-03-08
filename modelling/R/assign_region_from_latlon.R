# Assign Region categories using MEOW ecoregions (Marine Ecoregions of the World).
#
# Maps lon/lat points to the project's broad Regions by spatially intersecting
# points with grouped MEOW ecoregion polygons (same grouping as supplement map).

meow_region_shapes <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)

    if (!requireNamespace("sf", quietly = TRUE)) stop("Package 'sf' is required.")
    if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' is required.")

    shp_path <- "data/MEOW/meow_ecos.shp"
    if (!file.exists(shp_path) && requireNamespace("here", quietly = TRUE)) {
      shp_path <- here::here("data/MEOW/meow_ecos.shp")
    }
    if (!file.exists(shp_path)) stop("MEOW shapefile not found at: ", shp_path)

    ecos <- sf::st_read(shp_path, quiet = TRUE)
    if (!"ECOREGION" %in% names(ecos)) stop("MEOW shapefile missing ECOREGION column: ", shp_path)

    groups <- list(
      "Mediterranean Sea" = c("Western Mediterranean", "Alboran Sea", "Levantine Sea", "Ionian Sea", "Aegean Sea"),
      "North European Atlantic" = c("Celtic Seas", "North Sea"),
      "South European Atlantic" = c("South European Atlantic Shelf"),
      "Baltic Sea" = c("Baltic Sea"),
      "Black Sea" = c("Black Sea")
    )

    region_shapes <- dplyr::bind_rows(lapply(names(groups), function(r) {
      g <- groups[[r]]
      geom <- sf::st_union(sf::st_geometry(ecos[ecos$ECOREGION %in% g, ]))
      sf::st_sf(region = r, geometry = geom)
    }))

    cache <<- region_shapes
    region_shapes
  }
})

#' Assign Region based on lon/lat using MEOW ecoregions.
#'
#' @param dat Data frame with numeric columns longitude and latitude.
#' @return Same data frame with a character column region (if not present) or
#'         with missing values in region filled where possible.
assign_region_from_latlon <- function(dat) {
  if (!all(c("longitude", "latitude") %in% names(dat))) {
    stop("assign_region_from_latlon: dat must have 'longitude' and 'latitude' columns.")
  }
  if (!requireNamespace("sf", quietly = TRUE)) stop("Package 'sf' is required.")

  dat <- as.data.frame(dat)
  dat$longitude <- as.numeric(dat$longitude)
  dat$latitude <- as.numeric(dat$latitude)
  if (!"region" %in% names(dat)) dat$region <- NA_character_

  # only assign where missing
  need <- is.na(dat$region) & is.finite(dat$longitude) & is.finite(dat$latitude)
  if (!any(need)) return(dat)

  pts <- sf::st_as_sf(dat[need, , drop = FALSE], coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)
  shapes <- meow_region_shapes()

  hits <- sf::st_intersects(pts, shapes)
  assigned <- vapply(hits, function(ix) if (length(ix) > 0) shapes$region[ix[1]] else NA_character_, character(1))
  dat$region[need] <- assigned

  dat
}