# Efficient extraction of environmental covariates from NetCDF raster files
# and creation of prediction grids from raster data
#
# This replaces the inefficient point-by-point extraction in 
# data/VRE/Environmental_covariates_matching_VRE.R
#
# Key improvements:
# - Vectorized extraction (much faster)
# - Creates prediction grids directly from raster data (ocean cells only)
# - No predictions on land (raster data naturally excludes land)
#
# Usage:
#   source("modelling/extract_covariates_from_rasters.R")
#   # Then use extract_covariates_at_points() or create_prediction_grid_from_rasters()
#
# Required packages: here, RNetCDF, dplyr, terra

if (!requireNamespace("here", quietly = TRUE)) {
  stop("Package 'here' is required. Install with: install.packages('here')")
}
if (!requireNamespace("RNetCDF", quietly = TRUE)) {
  stop("Package 'RNetCDF' is required. Install with: install.packages('RNetCDF')")
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
  stop("Package 'dplyr' is required. Install with: install.packages('dplyr')")
}
if (!requireNamespace("terra", quietly = TRUE)) {
  stop("Package 'terra' is required. Install with: install.packages('terra')")
}

library(here)
library(dplyr)
library(RNetCDF)



# Base path for final rasters (used by default)
RASTER_DIR <- "data/env_rasters"

#' Inspect a single NetCDF file and return the first 2D data variable and its lat/lon dimension names
inspect_nc_file <- function(nc_path) {
  nc <- RNetCDF::open.nc(nc_path)
  on.exit(RNetCDF::close.nc(nc), add = TRUE)
  nvars <- RNetCDF::file.inq.nc(nc)$nvars
  for (j in seq_len(nvars) - 1L) {
    v <- RNetCDF::var.inq.nc(nc, j)
    if (length(v$dimids) != 2L) next
    dim_names <- vapply(v$dimids, function(d) RNetCDF::dim.inq.nc(nc, d)$name, character(1L))
    if (!any(grepl("lat", dim_names, ignore.case = TRUE)) ||
        !any(grepl("lon", dim_names, ignore.case = TRUE))) next
    lat_var <- dim_names[grepl("lat", dim_names, ignore.case = TRUE)][1L]
    lon_var <- dim_names[grepl("lon", dim_names, ignore.case = TRUE)][1L]
    return(list(varname = v$name, lat_var = lat_var, lon_var = lon_var))
  }
  NULL
}

#' Build COVARIATE_CONFIG from all NetCDF rasters in a directory (auto-discover variable and dim names)
build_covariate_config_from_dir <- function(raster_dir = RASTER_DIR) {
  base <- here::here(raster_dir)
  if (!dir.exists(base)) stop("Raster directory not found: ", base)
  files <- list.files(base, pattern = "[.]nc$", full.names = FALSE)
  if (length(files) == 0L) return(list())
  out <- list()
  n_files <- length(files)
  cat("Building covariate configuration from", n_files, "raster(s) in directory:", raster_dir, "\n")
  for (i in seq_along(files)) {
    f <- files[i]
    cat("  ", i, "/", n_files, " ", f, "\n", sep = "")
    path <- file.path(base, f)
    info <- tryCatch(inspect_nc_file(path), error = function(e) {
      warning("Skipping ", f, ": ", conditionMessage(e))
      NULL
    })
    if (is.null(info)) next
    name <- gsub("[.]nc$", "", f, ignore.case = TRUE)
    name <- gsub("-", "_", name, fixed = TRUE)
    out[[name]] <- list(
      file = file.path(raster_dir, f),
      varname = info$varname,
      lat_var = info$lat_var,
      lon_var = info$lon_var
    )
  }
  # Normalize configuration names to lower-case so that everything downstream
  # (extracted data columns, requested covariates, etc.) can use lower-case
  # consistently.
  if (length(out) > 0L) {
    names(out) <- tolower(names(out))
  }
  cat("  ", length(out), " raster(s) configured.\n", sep = "")
  out
}

# out <- build_covariate_config_from_dir("/Users/rt582/Library/CloudStorage/OneDrive-UniversityofCambridge/cambridge/phd/Paper_Conferences/seagrass/data/env_rasters")

# # Auto-build config from all rasters in env_rasters (one entry per .nc file)
COVARIATE_CONFIG <- build_covariate_config_from_dir(RASTER_DIR)



#' Find nearest grid index for a vector of coordinates (vectorized)
#' Handles grid_vals that are unsorted or contain NAs (uses only non-NA, sorts for findInterval).
nearest_grid_index <- function(coords, grid_vals) {
  grid_vals <- as.numeric(grid_vals)
  n <- length(grid_vals)
  if (n == 1L) return(rep(1L, length(coords)))
  # findInterval requires sorted non-decreasing and no NAs
  valid <- !is.na(grid_vals)
  if (!any(valid)) stop("grid_vals are all NA")
  valid_idx <- which(valid)
  grid_clean <- grid_vals[valid_idx]
  ord <- order(grid_clean)
  sorted <- grid_clean[ord]
  orig_idx <- valid_idx[ord]
  n_sorted <- length(sorted)
  if (n_sorted == 1L) return(rep(orig_idx[1L], length(coords)))
  i <- findInterval(coords, sorted, left.open = FALSE)
  i <- pmax(1L, pmin(i, n_sorted - 1L))
  i_hi <- i + 1L
  i_lo_orig <- orig_idx[i]
  i_hi_orig <- orig_idx[i_hi]
  d_lo <- abs(coords - grid_vals[i_lo_orig])
  d_hi <- abs(coords - grid_vals[i_hi_orig])
  ifelse(d_lo <= d_hi, i_lo_orig, i_hi_orig)
}

#' Extract values from a NetCDF raster at given point locations
#'
#' @param nc_file Path to NetCDF file (relative to project root)
#' @param varname Name of variable in NetCDF file
#' @param points Data frame with columns 'latitude' and 'longitude'
#' @param lat_var Name of latitude variable in NetCDF file (ignored when using terra; kept for backwards compatibility)
#' @param lon_var Name of longitude variable in NetCDF file (ignored when using terra; kept for backwards compatibility)
#' @param use_closest If TRUE, attempt to replace NA values with the nearest non-NA raster value (using terra's search_radius)
#' @param method Extraction method: "nearest" (default; terra's "simple") or "bilinear" for bilinear interpolation
#' @return Vector of extracted values (same length as nrow(points))
extract_from_nc <- function(nc_file, varname, lat_var, lon_var, points, use_closest = TRUE, method = "nearest") {
  nc_path <- here::here(nc_file)
  if (!file.exists(nc_path)) {
    stop("NetCDF file not found: ", nc_path)
  }

  # Use terra for raster I/O and interpolation
  r <- terra::rast(nc_path)

  # If a specific variable name is provided and exists as a layer, select it.
  if (!is.null(varname) && nzchar(varname) && varname %in% names(r)) {
    r <- r[[varname]]
  } else if (terra::nlyr(r) > 1) {
    # Fall back to first layer but warn once if varname is not found.
    warning("Variable '", varname, "' not found in ", nc_file, 
            "; using first layer in file. Check COVARIATE_CONFIG if this is unintended.")
    r <- r[[1]]
  }

  # Prepare coordinate matrix (lon, lat)
  if (!all(c("longitude", "latitude") %in% names(points))) {
    stop("points must have 'latitude' and 'longitude' columns")
  }
  pts <- cbind(points$longitude, points$latitude)

  terra_method <- if (identical(tolower(method), "bilinear")) "bilinear" else "simple"

  # Primary extraction (matrix interface; no ID argument here)
  vals_mat <- terra::extract(r, pts, method = terra_method)
  extracted <- as.numeric(vals_mat[, 1])

  # Optionally fill NA values with nearest non-NA cell using FNN on terra cell centers
  if (use_closest && any(is.na(extracted)) && requireNamespace("FNN", quietly = TRUE)) {
    na_idx <- which(is.na(extracted))

    # Get all non-NA cell centers and values from the raster
    vals_all <- terra::values(r, mat = FALSE)
    valid_cells <- which(!is.na(vals_all))
    if (length(valid_cells) > 0) {
      valid_xy <- terra::xyFromCell(r, valid_cells)
      valid_vals <- vals_all[valid_cells]

      # Nearest-neighbour in geographic space from NA points to valid cells
      nn <- FNN::get.knnx(
        data = valid_xy,
        query = pts[na_idx, , drop = FALSE],
        k = 1
      )$nn.index[, 1]

      extracted[na_idx] <- valid_vals[nn]
    }
  }

  extracted
}


#' Extract all covariates at point locations
#'
#' @param points Data frame with columns 'latitude' and 'longitude'
#' @param covariates Character vector of covariate names to extract (default: all)
#' @param use_closest If TRUE, use closest valid value when exact match is NA
#' @param method Extraction method: "nearest" (default) or "bilinear" for bilinear interpolation
#' @param progress_callback Optional function(i, n, name) called after each covariate (e.g. for progress bar)
#' @return Data frame with original points plus extracted covariate columns
extract_covariates_at_points <- function(points, 
                                         covariates = names(COVARIATE_CONFIG),
                                         use_closest = TRUE,
                                         method = "nearest",
                                         progress_callback = NULL) {
  if (!all(c("latitude", "longitude") %in% names(points))) {
    stop("points must have 'latitude' and 'longitude' columns")
  }
  
  # Work in lower-case consistently and filter to known covariates
  covariates <- tolower(covariates)
  covariates <- covariates[covariates %in% names(COVARIATE_CONFIG)]
  n_cov <- length(covariates)
  
  result <- points
  
  for (i in seq_along(covariates)) {
    covar_name <- covariates[i]
    config <- COVARIATE_CONFIG[[covar_name]]
    
    if (is.null(progress_callback)) {
      cat("Extracting", covar_name, "...\n")
    }
    
    result[[covar_name]] <- extract_from_nc(
      nc_file = config$file,
      varname = config$varname,
      lat_var = config$lat_var,
      lon_var = config$lon_var,
      points = points,
      use_closest = use_closest,
      method = method
    )
    
    if (is.function(progress_callback)) {
      progress_callback(i, n_cov, covar_name)
    }
  }
  
  # Note: Negative values are kept as-is (not converted to NA)
  # If you need to filter negative values, do it downstream in your analysis
  # Previous behavior converted negative remote-sensing values to NA, but this
  # was removed per user request to preserve original raster values

  result
}


#' Create an arbitrary prediction grid and extract covariates from rasters
#'
#' Creates a regular grid over specified lon/lat ranges and extracts covariates
#' from raster data. Only grid cells where at least one raster has valid data
#' are included (naturally excludes land areas). Optionally caches the result
#' to disk so subsequent calls with the same parameters load instantly.
#'
#' @param lon_range Length-2 vector (min, max) longitude
#' @param lat_range Length-2 vector (min, max) latitude
#' @param n_lon Number of longitude points
#' @param n_lat Number of latitude points
#' @param covariates Character vector of covariate names to include in grid
#' @param reference_raster Name of covariate to use for filtering valid cells
#'   (cells with NA in this raster will be excluded)
#' @param cache_path If non-NULL, path to save/load the grid as RDS. When
#'   \code{use_cache} is TRUE and this file exists with matching parameters,
#'   the grid is loaded from cache instead of recomputing.
#' @param use_cache If TRUE and \code{cache_path} is set and the cache file
#'   exists with matching parameters, load from cache. If FALSE, always
#'   recompute (and overwrite cache if \code{cache_path} is set).
#' @param show_progress If TRUE, show a progress bar while extracting covariates
#' @param method Extraction method: "nearest" (default) or "bilinear" for bilinear interpolation
#' @param fill_coastal If TRUE (default), use nearest non-NA values to fill coastal/ocean NAs
#' @return Data frame with columns: longitude, latitude, and all requested covariates
create_prediction_grid_from_rasters <- function(lon_range,
                                                lat_range,
                                                covariates,
                                                n_lon = 200,
                                                n_lat = 200,
                                                reference_raster = "KD_closest",
                                                cache_path = NULL,
                                                use_cache = TRUE,
                                                show_progress = TRUE,
                                                method = "nearest",
                                                fill_coastal = TRUE) {
  # Work in lower-case consistently and align with COVARIATE_CONFIG
  covariates <- tolower(covariates)
  covariates <- covariates[covariates %in% names(COVARIATE_CONFIG)]

  reference_raster <- tolower(reference_raster)
  # Resolve legacy reference_raster name (e.g. "kd_closest") to actual config name if needed
  if (!reference_raster %in% covariates) {
    legacy <- list(kd_closest = "kd_490|kd490", bathy = "bathy")
    for (alias in names(legacy)) {
      if (reference_raster == alias) {
        pat <- legacy[[alias]]
        match <- covariates[grepl(pat, covariates, ignore.case = TRUE)][1L]
        if (!is.na(match)) {
          reference_raster <- match
          if (!reference_raster %in% covariates) covariates <- c(covariates, reference_raster)
        }
        break
      }
    }
  }
  params <- list(
    lon_range = lon_range,
    lat_range = lat_range,
    n_lon = n_lon,
    n_lat = n_lat,
    covariates = covariates,
    reference_raster = reference_raster
  )

  # Try load from cache
  if (is.character(cache_path) && length(cache_path) == 1L && nzchar(cache_path) && use_cache) {
    cache_file <- here::here(cache_path)
    if (file.exists(cache_file)) {
      cached <- tryCatch(readRDS(cache_file), error = function(e) NULL)
      if (is.list(cached) &&
          identical(names(cached), c("grid", "params")) &&
          is.data.frame(cached$grid) &&
          is.list(cached$params)) {
        if (isTRUE(all.equal(cached$params, params))) {
          cat("Loaded prediction grid from cache:", cache_file, "(", nrow(cached$grid), "cells)\n")
          return(cached$grid)
        }
        warning("Cached grid parameters differ from request; recomputing and updating cache.")
      } else {
        warning("Cache file invalid or outdated; recomputing and updating cache.")
      }
    }
  }

  # Create regular grid
  lon_seq <- seq(lon_range[1], lon_range[2], length.out = n_lon)
  lat_seq <- seq(lat_range[1], lat_range[2], length.out = n_lat)

  grid <- expand.grid(longitude = lon_seq, latitude = lat_seq)
  cat("Created", nrow(grid), "grid cells\n")

  # Progress bar: one tick per covariate
  n_cov <- length(covariates)

  progress_callback <- NULL
  pb <- NULL
  if (show_progress && n_cov > 0) {
    pb <- utils::txtProgressBar(min = 0, max = n_cov, style = 3, width = 50)
    progress_callback <- function(i, n, name) {
      utils::setTxtProgressBar(pb, i)
    }
  }

  # Extract all covariates.
  # If fill_coastal = TRUE, use_closest = TRUE so coarse fields are extended toward coasts
  # using nearest non-NA values (implemented via FNN on terra cell centers).
  grid <- extract_covariates_at_points(
    grid,
    covariates = covariates,
    use_closest = fill_coastal,
    method = method,
    progress_callback = progress_callback
  )

  if (!is.null(pb)) {
    close(pb)
    cat("\n")
  }

  # Filter to cells where reference raster has valid data (excludes land)
  if (reference_raster %in% names(grid)) {
    valid_mask <- !is.na(grid[[reference_raster]])
    dropped_mask <- !valid_mask
    num_dropped <- sum(dropped_mask)
    grid <- grid[valid_mask, , drop = FALSE]
    if (num_dropped > 0) {
      cat(sprintf(
        "Filtered out %d grid cells with NA in '%s'. Example dropped cell coordinates (up to 5):\n", 
        num_dropped, reference_raster
      ))
      example_rows <- which(dropped_mask)[seq_len(min(num_dropped, 5))]
      print(grid[FALSE, c("longitude", "latitude")]) # empty placeholder for format
      # Because unfiltered grid is gone, reload example coordinates from a pre-subset copy
      # (rebuild mask for original)
      all_coords <- expand.grid(longitude = unique(grid$longitude), latitude = unique(grid$latitude))
      suppressWarnings({
        # Try to reproduce original grid (for display of dropped points)
        if (exists("lon_seq") && exists("lat_seq")) {
          full_grid <- expand.grid(longitude = lon_seq, latitude = lat_seq)
          dropped_coords <- full_grid[which(dropped_mask)[seq_len(min(num_dropped, 5))], , drop = FALSE]
          print(dropped_coords)
        }
      })
    }
    cat("Filtered to", nrow(grid), "cells with valid", reference_raster, "data (ocean only)\n")
  } else {
    warning("Reference raster '", reference_raster, "' not found in grid. Keeping all cells.")
  }

  # Save to cache if requested
  if (is.character(cache_path) && length(cache_path) == 1L && nzchar(cache_path)) {
    cache_file <- here::here(cache_path)
    cache_dir <- dirname(cache_file)
    if (!dir.exists(cache_dir)) {
      dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    }
    saveRDS(list(grid = grid, params = params), cache_file)
    cat("Saved prediction grid to cache:", cache_file, "\n")
  }

  grid
}


#' Verify extracted covariates match existing data
#'
#' Compares extracted values with existing data to ensure correctness.
#'
#' @param existing_data Data frame with latitude, longitude, and covariate columns
#' @param covariates Character vector of covariate names to verify
#' @param tolerance Numeric tolerance for comparison (default 1e-6)
#' @param max_mismatches Maximum number of mismatches to report (default 100)
#' @return List with verification results
verify_extraction <- function(existing_data,
                              covariates = names(COVARIATE_CONFIG),
                              tolerance = 1e-6,
                              max_mismatches = 100) {
  if (!all(c("latitude", "longitude") %in% names(existing_data))) {
    stop("existing_data must have 'latitude' and 'longitude' columns")
  }
  
  # Extract covariates for existing points
  cat("Re-extracting covariates for", nrow(existing_data), "points...\n")
  re_extracted <- extract_covariates_at_points(
    existing_data[, c("latitude", "longitude")],
    covariates = covariates,
    use_closest = TRUE
  )
  
  # Compare
  results <- list()
  mismatches <- list()
  
  for (covar in covariates) {
    if (!covar %in% names(existing_data)) {
      warning("Covariate '", covar, "' not found in existing_data. Skipping.")
      next
    }
    
    existing_vals <- existing_data[[covar]]
    new_vals <- re_extracted[[covar]]
    
    # Handle NAs
    both_na <- is.na(existing_vals) & is.na(new_vals)
    both_valid <- !is.na(existing_vals) & !is.na(new_vals)
    existing_na_new_val <- is.na(existing_vals) & !is.na(new_vals)
    existing_val_new_na <- !is.na(existing_vals) & is.na(new_vals)
    
    # Compare non-NA values
    if (sum(both_valid) > 0) {
      diffs <- abs(existing_vals[both_valid] - new_vals[both_valid])
      matches <- diffs <= tolerance
      n_match <- sum(matches)
      n_diff <- sum(!matches)
      
      results[[covar]] <- list(
        n_total = length(existing_vals),
        n_both_na = sum(both_na),
        n_both_valid = sum(both_valid),
        n_match = n_match,
        n_mismatch = n_diff,
        pct_match = 100 * n_match / sum(both_valid),
        max_diff = if (n_diff > 0) max(diffs[!matches], na.rm = TRUE) else 0,
        mean_diff = if (n_diff > 0) mean(diffs[!matches], na.rm = TRUE) else 0
      )
      
      # Store mismatches
      if (n_diff > 0) {
        mismatch_idx <- which(both_valid)[!matches]
        n_report <- min(n_diff, max_mismatches)
        mismatches[[covar]] <- data.frame(
          row = mismatch_idx[1:n_report],
          latitude = existing_data$latitude[mismatch_idx[1:n_report]],
          longitude = existing_data$longitude[mismatch_idx[1:n_report]],
          existing = existing_vals[mismatch_idx[1:n_report]],
          re_extracted = new_vals[mismatch_idx[1:n_report]],
          diff = diffs[!matches][1:n_report]
        )
      }
    } else {
      results[[covar]] <- list(
        n_total = length(existing_vals),
        n_both_na = sum(both_na),
        n_both_valid = 0,
        n_match = 0,
        n_mismatch = 0,
        pct_match = NA,
        max_diff = NA,
        mean_diff = NA
      )
    }
    
    # Report NA mismatches
    if (sum(existing_na_new_val) > 0 || sum(existing_val_new_na) > 0) {
      results[[covar]]$n_na_mismatch <- sum(existing_na_new_val) + sum(existing_val_new_na)
      results[[covar]]$n_existing_na_new_val <- sum(existing_na_new_val)
      results[[covar]]$n_existing_val_new_na <- sum(existing_val_new_na)
    }
  }
  
  list(results = results, mismatches = mismatches)
}
