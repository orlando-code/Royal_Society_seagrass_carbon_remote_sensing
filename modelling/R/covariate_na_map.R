# Count per grid cell how many base raster covariates are NA (no coastal fill).
# Depends on modelling/R/extract_covariates_from_rasters.R (RASTER_CONFIG, raster_covariates).

#' Build a lon/lat grid and count NAs across covariate rasters at each cell
#'
#' Uses [extract_covariates_at_points()] with \code{use_closest = FALSE} so values
#' reflect true missingness in each NetCDF layer (no nearest-neighbour fill).
#'
#' @param lon_range,lat_range Longitude/latitude ranges (length-2).
#' @param n_lon,n_lat Grid dimensions (\code{length.out} for \code{seq}).
#' @param covariates Character vector of covariate names; default all
#'   \code{raster_covariates} (typically 42 layers under \code{data/env_rasters}).
#' @param use_closest Must stay \code{FALSE} for an NA diagnostic map.
#' @param method Passed to \code{extract_from_nc} (e.g. \code{"nearest"}).
#' @param drop_all_na If \code{TRUE}, drop rows where every covariate is NA (e.g. land).
#'   Keep \code{FALSE} when plotting with \code{geom_raster} so the lon/lat grid stays regular.
#' @param cache_path If set, save/load RDS with list(grid, params).
#' @param use_cache Whether to read cache when params match.
#' @param show_progress Progress bar over covariates.
#' @return \code{data.frame} with \code{longitude}, \code{latitude}, \code{n_na},
#'   \code{n_covariates}, and one column per covariate.
build_covariate_na_count_grid <- function(lon_range,
                                          lat_range,
                                          n_lon = 200L,
                                          n_lat = 160L,
                                          covariates = NULL,
                                          use_closest = FALSE,
                                          method = "nearest",
                                          drop_all_na = FALSE,
                                          cache_path = NULL,
                                          use_cache = TRUE,
                                          show_progress = TRUE) {
  if (!exists("RASTER_CONFIG") || !exists("raster_covariates")) {
    stop("Source modelling/R/extract_covariates_from_rasters.R first.")
  }
  if (is.null(covariates)) {
    covariates <- raster_covariates
  }
  covariates <- unique(tolower(covariates))
  miss <- setdiff(covariates, names(RASTER_CONFIG))
  if (length(miss) > 0L) {
    stop("Unknown covariate(s): ", paste(miss, collapse = ", "))
  }

  lon_range <- range(lon_range)
  lat_range <- range(lat_range)
  params <- list(
    lon_range = lon_range,
    lat_range = lat_range,
    n_lon = as.integer(n_lon)[1L],
    n_lat = as.integer(n_lat)[1L],
    covariates = covariates,
    use_closest = use_closest,
    method = method,
    version = 1L
  )

  if (is.character(cache_path) && nzchar(cache_path) && isTRUE(use_cache)) {
    cache_file <- if (requireNamespace("here", quietly = TRUE)) {
      here::here(cache_path)
    } else {
      cache_path
    }
    if (file.exists(cache_file)) {
      cached <- tryCatch(readRDS(cache_file), error = function(e) NULL)
      if (is.list(cached) && is.data.frame(cached$grid) && identical(cached$params, params)) {
        message("Loaded NA-count grid from cache: ", cache_file)
        return(cached$grid)
      }
      message("Cache miss or params differ; recomputing NA-count grid.")
    }
  }

  lon_seq <- seq(lon_range[1L], lon_range[2L], length.out = params$n_lon)
  lat_seq <- seq(lat_range[1L], lat_range[2L], length.out = params$n_lat)
  grid <- expand.grid(longitude = lon_seq, latitude = lat_seq)

  n_cov <- length(covariates)
  pb <- NULL
  progress_callback <- NULL
  if (isTRUE(show_progress) && n_cov > 0L) {
    pb <- utils::txtProgressBar(min = 0, max = n_cov, style = 3, width = 50)
    progress_callback <- function(i, n, name) utils::setTxtProgressBar(pb, i)
  }

  grid <- extract_covariates_at_points(
    grid,
    covariates = covariates,
    use_closest = use_closest,
    method = method,
    progress_callback = progress_callback
  )
  if (!is.null(pb)) close(pb)

  cov_cols <- covariates[covariates %in% names(grid)]
  mat <- as.matrix(grid[, cov_cols, drop = FALSE])
  grid$n_na <- rowSums(is.na(mat))
  grid$n_covariates <- length(cov_cols)

  if (isTRUE(drop_all_na)) {
    grid <- grid[grid$n_na < grid$n_covariates, , drop = FALSE]
  }

  if (is.character(cache_path) && nzchar(cache_path)) {
    cache_file <- if (requireNamespace("here", quietly = TRUE)) {
      here::here(cache_path)
    } else {
      cache_path
    }
    dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE)
    saveRDS(list(grid = grid, params = params), cache_file)
    message("Saved NA-count grid cache: ", cache_file)
  }

  grid
}
