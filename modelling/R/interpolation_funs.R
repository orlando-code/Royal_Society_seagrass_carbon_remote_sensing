# Interpolation helper functions for ordinary kriging and related workflows.
# Use with any data subset (e.g. dat, or dat %>% filter(Region == "Mediterranean Sea")).


#' Aggregate to one row per location (required for kriging; duplicate locations cause singular matrix).
#'
#' @param dat Data frame with coordinate and value columns.
#' @param coords Character vector of coordinate column names (e.g. c("longitude", "latitude")).
#' @param value_var Name of the variable to aggregate (e.g. "carbon_density_g_c_cm3").
#' @param fun Function to aggregate (default median).
#' @param extra_vars Optional character vector of extra columns to keep; aggregated by median per location.
#' @return Data frame with one row per unique (coords) and aggregated value; NA values dropped.
aggregate_to_unique_locations <- function(dat,
                                          coords = c("longitude", "latitude"),
                                          value_var = "carbon_density_g_c_cm3",
                                          fun = median,
                                          extra_vars = NULL) {
  extra_vars <- intersect(if (is.null(extra_vars)) character(0) else extra_vars, names(dat))
  if (length(extra_vars) == 0) {
    out <- dat %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(coords))) %>%
      dplyr::summarise(
        !!value_var := fun(.data[[value_var]], na.rm = TRUE),
        .groups = "drop"
      )
  } else {
    out <- dat %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(coords))) %>%
      dplyr::summarise(
        !!value_var := fun(.data[[value_var]], na.rm = TRUE),
        dplyr::across(dplyr::all_of(extra_vars), ~ if (is.numeric(.)) median(., na.rm = TRUE) else .[1]),
        .groups = "drop"
      )
  }
  out %>% dplyr::filter(!is.na(.data[[value_var]]))
}


#' Compute empirical variogram.
#'
#' @param dat_locations Data frame with coordinates and one value column.
#' @param formula Formula (e.g. carbon_density_g_c_cm3 ~ 1).
#' @param coords Formula for locations (e.g. ~ longitude + latitude).
#' @return gstat variogram object.
empirical_variogram <- function(dat_locations, formula, coords = ~ longitude + latitude) {
  gstat::variogram(formula, data = dat_locations, locations = coords)
}


#' Check if a fitted vgm is valid (finite psill/range, positive structure ranges, no try-error).
is_valid_variogram_fit <- function(fit) {
  basic_ok <- !inherits(fit, "try-error") &&
    is.data.frame(fit) &&
    all(c("psill", "range") %in% names(fit)) &&
    all(!is.na(fit$psill)) && all(!is.na(fit$range)) &&
    all(is.finite(fit$psill)) && all(is.finite(fit$range)) &&
    nrow(fit) >= 2L
  if (!basic_ok) return(FALSE)

  # gstat / geoR can occasionally return negative ranges for some trial fits,
  # which then cause downstream errors like "variogram range can never be negative".
  # Treat any non-nugget structure with range <= 0 as invalid.
  structure_rows <- fit$model != "Nug"
  if (any(structure_rows)) {
    all(fit$range[structure_rows] > 0)
  } else {
    TRUE
  }
}


#' Build default trial variogram models from empirical variogram summaries.
default_trial_models <- function(vario) {
  sill_guess <- max(vario$gamma, na.rm = TRUE)
  nugget_guess <- min(vario$gamma, na.rm = TRUE)
  range_guess <- max(vario$dist, na.rm = TRUE) / 3
  range_med <- median(vario$dist, na.rm = TRUE)
  list(
    gstat::vgm(psill = sill_guess * 0.5, "Exp", range = range_guess, nugget = nugget_guess),
    gstat::vgm(psill = sill_guess, "Exp", range = range_med, nugget = nugget_guess),
    gstat::vgm(psill = sill_guess * 0.5, "Sph", range = range_guess, nugget = nugget_guess),
    gstat::vgm(psill = sill_guess, "Sph", range = range_med, nugget = nugget_guess),
    gstat::vgm(psill = sill_guess * 0.5, "Gau", range = range_guess, nugget = nugget_guess),
    gstat::vgm(psill = sill_guess, "Gau", range = range_med, nugget = nugget_guess),
    gstat::vgm(psill = sill_guess * 0.5, "Mat", range = range_guess, nugget = nugget_guess, kappa = 1.5),
    gstat::vgm(psill = sill_guess, "Mat", range = range_med, nugget = nugget_guess, kappa = 0.5),
    gstat::vgm("Exp"),
    gstat::vgm("Sph"),
    gstat::vgm("Gau")
  )
}


#' Fit multiple variogram models and select best by Cressie WSS.
#'
#' @param vario Empirical variogram from gstat::variogram.
#' @param trial_models List of vgm() initial models; if NULL, uses default_trial_models(vario).
#' @param verbose If TRUE, print comparison table and best model.
#' @return List with vario_fit (best vgm or NA), comparison_df, best_meta, vario.
fit_variogram_best <- function(vario, trial_models = NULL, verbose = TRUE) {
  if (is.null(trial_models)) trial_models <- default_trial_models(vario)
  fitted_list <- list()
  fit_metadata <- list()
  for (mod in trial_models) {
    fit <- try(gstat::fit.variogram(vario, model = mod, fit.method = 7), silent = TRUE)
    if (!is_valid_variogram_fit(fit)) next
    gamma_fit <- gstat::variogramLine(fit, dist_vector = vario$dist)$gamma
    gamma_emp <- vario$gamma
    np_bin <- vario$np
    gamma_fit_safe <- pmax(gamma_fit, 1e-6)
    wss <- sum((np_bin / gamma_fit_safe^2) * (gamma_emp - gamma_fit)^2, na.rm = TRUE)
    rss <- sum((gamma_emp - gamma_fit)^2, na.rm = TRUE)
    model_label <- paste(fit$model[fit$model != "Nug"], collapse = "+")
    fitted_list <- append(fitted_list, list(fit))
    fit_metadata <- append(fit_metadata, list(list(
      model = model_label, wss = wss, rss = rss,
      nugget = fit$psill[1], sill = sum(fit$psill), range = max(fit$range)
    )))
  }
  vario_fit <- NA
  comparison_df <- NULL
  best_meta <- NULL
  if (length(fitted_list) > 0) {
    wss_all <- vapply(fit_metadata, function(m) m$wss, numeric(1))
    best_idx <- which.min(wss_all)
    vario_fit <- fitted_list[[best_idx]]
    best_meta <- fit_metadata[[best_idx]]
    comparison_df <- data.frame(
      model = vapply(fit_metadata, function(m) m$model, character(1)),
      WSS = vapply(fit_metadata, function(m) m$wss, numeric(1)),
      RSS = vapply(fit_metadata, function(m) m$rss, numeric(1))
    )
    comparison_df <- comparison_df[order(comparison_df$WSS), ]
    rownames(comparison_df) <- NULL
    if (verbose) {
      cat("Variogram model comparison (Cressie WSS, lower is better):\n")
      print(comparison_df)
      cat("\nBest model:", best_meta$model, "(WSS =", round(best_meta$wss, 6), ")\n")
    }
  } else if (verbose) {
    warning("Variogram fit did not converge for any model or starting values.")
  }
  list(vario_fit = vario_fit, comparison_df = comparison_df, best_meta = best_meta, vario = vario)
}


#' Build a prediction grid.
#'
#' @param lon_range Length-2 vector (min, max) longitude.
#' @param lat_range Length-2 vector (min, max) latitude.
#' @param n_lon Number of longitude points.
#' @param n_lat Number of latitude points.
#' @param coords Character vector of coordinate names (default c("longitude", "latitude")).
#' @return Data frame with columns coords[1], coords[2].
build_prediction_grid <- function(lon_range, lat_range, n_lon = 200, n_lat = 200,
                                  coords = c("longitude", "latitude")) {
  gr <- expand.grid(
    x = seq(lon_range[1], lon_range[2], length.out = n_lon),
    y = seq(lat_range[1], lat_range[2], length.out = n_lat)
  )
  names(gr) <- coords
  gr
}


#' Validate variogram and data before kriging (structure rows must have positive range).
validate_kriging_inputs <- function(vario_fit, prediction_grid, dat_locations,
                                    value_var, coords = c("longitude", "latitude")) {
  if (!identical(vario_fit, NA) && is.data.frame(vario_fit)) {
    has_na <- any(is.na(vario_fit$psill)) || any(is.na(vario_fit$range))
    structure_rows <- vario_fit$model != "Nug"
    bad_range <- any(vario_fit$range[structure_rows] <= 0, na.rm = TRUE)
    if (has_na || bad_range)
      stop("Invalid variogram: psill/range NA or structure range <= 0.")
  }
  if (any(is.na(prediction_grid[[coords[1]]])) || any(is.na(prediction_grid[[coords[2]]])))
    stop("prediction_grid contains NA coordinates.")
  if (any(is.na(dat_locations[[coords[1]]])) || any(is.na(dat_locations[[coords[2]]])) ||
      any(is.na(dat_locations[[value_var]])))
    stop("dat_locations contains NA in coordinates or value variable.")
  invisible(TRUE)
}


#' Run ordinary kriging and return predictions/variance as vectors.
#'
#' Uses local kriging (nmax) to avoid singular covariance matrices. Optionally
#' enforces a minimum nugget on the variogram for numerical stability.
#'
#' @param dat_locations Data frame with one row per location (coordinates + value_var).
#' @param vario_fit Fitted vgm from fit_variogram_best.
#' @param prediction_grid Data frame with coordinate columns.
#' @param value_var Name of response variable.
#' @param coords Character vector of coordinate names.
#' @param nmax Maximum number of nearest observations to use per prediction (default 25). Smaller values reduce singular matrix errors.
#' @param min_nugget_ratio Minimum nugget as fraction of sill (default 1e-4). Ensures positive definite covariance.
#' @return List with kriged (full krige result), pred (vector), var (vector), success (logical).
run_ordinary_kriging <- function(dat_locations, vario_fit, prediction_grid,
                                 value_var = "carbon_density_g_c_cm3",
                                 coords = c("longitude", "latitude"),
                                 nmax = 25,
                                 min_nugget_ratio = 1e-4) {
  if (identical(vario_fit, NA) || !is.data.frame(vario_fit)) {
    return(list(kriged = NULL, pred = NULL, var = NULL, success = FALSE))
  }
  # Enforce minimum nugget for numerical stability (avoids singular covariance)
  sill_total <- sum(vario_fit$psill)
  nugget_current <- vario_fit$psill[1]
  nugget_min <- max(min_nugget_ratio * sill_total, 1e-10)
  if (nugget_current < nugget_min) {
    vario_fit$psill[1] <- nugget_min
  }
  # gstat::krige requires unique coordinates in data (duplicate locations → singular covariance → NA predictions)
  n_unique <- nrow(unique(dat_locations[, coords, drop = FALSE]))
  if (n_unique < nrow(dat_locations)) {
    data_for_krige <- aggregate_to_unique_locations(dat_locations, coords = coords, value_var = value_var)
  } else {
    data_for_krige <- dat_locations
  }
  formula <- as.formula(paste(value_var, "~ 1"))
  locations <- as.formula(paste("~", paste(coords, collapse = " + ")))
  kriged <- tryCatch(
    gstat::krige(formula, locations = locations, data = data_for_krige,
                 newdata = prediction_grid, model = vario_fit, nmax = nmax),
    error = function(e) {
      warning("Kriging failed: ", conditionMessage(e))
      return(NULL)
    }
  )
  if (is.null(kriged)) return(list(kriged = NULL, pred = NULL, var = NULL, success = FALSE))
  if (is.data.frame(kriged)) {
    pred <- kriged$var1.pred
    var <- kriged$var1.var
  } else if (!is.null(methods::slotNames(kriged)) && "data" %in% methods::slotNames(kriged)) {
    pred <- kriged@data$var1.pred
    var <- kriged@data$var1.var
  } else {
    return(list(kriged = kriged, pred = NULL, var = NULL, success = FALSE))
  }
  list(kriged = kriged, pred = pred, var = var, success = TRUE)
}


#' Extract kriged predictions at observation locations via nearest grid point.
#'
#' @param kriged_pred Vector of predictions on prediction_grid (same order as prediction_grid rows).
#' @param kriged_var Vector of kriging variances.
#' @param prediction_grid Data frame with coordinate columns.
#' @param dat_locations Data frame with observation coordinates.
#' @param coords Character vector of coordinate names.
#' @param max_dist_deg Maximum distance (degrees) to accept; beyond this, prediction set to NA.
#' @return List with pred (vector), var (vector), dist_to_nn (vector).
extract_predictions_at_locations <- function(kriged_pred, kriged_var, prediction_grid,
                                             dat_locations, coords = c("longitude", "latitude"),
                                             max_dist_deg = 1) {
  if (!requireNamespace("FNN", quietly = TRUE)) stop("Package FNN required for extract_predictions_at_locations.")
  nn_idx <- FNN::get.knnx(
    data = as.matrix(prediction_grid[, coords]),
    query = as.matrix(dat_locations[, coords]),
    k = 1
  )$nn.index[, 1]
  pred <- kriged_pred[nn_idx]
  var <- kriged_var[nn_idx]
  dist_to_nn <- sqrt(
    (dat_locations[[coords[1]]] - prediction_grid[[coords[1]]][nn_idx])^2 +
    (dat_locations[[coords[2]]] - prediction_grid[[coords[2]]][nn_idx])^2
  )
  pred[dist_to_nn > max_dist_deg] <- NA_real_
  var[dist_to_nn > max_dist_deg] <- NA_real_
  list(pred = pred, var = var, dist_to_nn = dist_to_nn)
}


#' Full ordinary kriging pipeline: aggregate → variogram → fit → krige → optional plots.
#'
#' Run with different subsets by passing e.g. dat or dat %>% filter(Region == "Mediterranean Sea").
#'
#' @param dat Data frame with coordinate and value columns (can have duplicate locations).
#' @param value_var Name of variable to predict (default "carbon_density_g_c_cm3").
#' @param coords Character vector of coordinate names (default c("longitude", "latitude")).
#' @param prediction_grid Optional; if NULL, built from data extent with n_lon, n_lat.
#' @param aggregate_to_unique If TRUE, aggregate to one row per (coords) before kriging; if FALSE (default), use all rows (duplicate locations with same coords, different target values are allowed).
#' @param n_lon,n_lat Grid size when prediction_grid is NULL.
#' @param make_plots If TRUE, create spatial and scatter plots.
#' @param save_plots If TRUE, save plots to files.
#' @param plot_prefix Path prefix for saved plots (e.g. "figures/interpolation/").
#' @param verbose If TRUE, print progress and variogram comparison.
#' @param kriging_nmax Maximum number of nearest observations per prediction (default 25). Use to avoid "covariance matrix singular" and all-NA predictions.
#' @return List with dat_locations, vario, fit_result (vario_fit, comparison_df, best_meta),
#'   prediction_grid, kriged, kriged_pred, kriged_var, dat_locations_with_pred (copy with predicted, predicted_error),
#'   kriged_sum, plots (list of ggplot objects if make_plots).
ordinary_kriging_pipeline <- function(dat,
                                     value_var = "carbon_density_g_c_cm3",
                                     coords = c("longitude", "latitude"),
                                     prediction_grid = NULL,
                                     aggregate_to_unique = FALSE,
                                     n_lon = 100,
                                     n_lat = 100,
                                     make_plots = TRUE,
                                     save_plots = TRUE,
                                     plot_prefix = "figures/interpolation/",
                                     verbose = TRUE,
                                     kriging_nmax = 25) {
  if (aggregate_to_unique) {
    dat_locations <- aggregate_to_unique_locations(dat, coords = coords, value_var = value_var)
    if (verbose) cat("Aggregated to", nrow(dat_locations), "unique locations (from", nrow(dat), "samples).\n")
  } else {
    dat_locations <- dat
    if (verbose) cat("Using", nrow(dat_locations), "samples (duplicate coordinates allowed).\n")
  }

  # Empirical variogram
  formula <- as.formula(paste(value_var, "~ 1"))
  locations <- as.formula(paste("~", paste(coords, collapse = " + ")))
  vario <- empirical_variogram(dat_locations, formula, locations)

  # Fit best variogram model
  fit_result <- fit_variogram_best(vario, trial_models = NULL, verbose = verbose)
  vario_fit <- fit_result$vario_fit

  # Prediction grid
  if (is.null(prediction_grid)) {
    lon_range <- range(dat_locations[[coords[1]]], na.rm = TRUE)
    lat_range <- range(dat_locations[[coords[2]]], na.rm = TRUE)
    prediction_grid <- build_prediction_grid(lon_range, lat_range, n_lon, n_lat, coords)
  }

  # Store original grid and filter out rows with NaN/NA/Inf coordinates before kriging
  # (kriging requires valid coordinates; we'll restore NA predictions for filtered rows afterward)
  prediction_grid_original <- prediction_grid
  coord_valid <- !is.na(prediction_grid[[coords[1]]]) & !is.na(prediction_grid[[coords[2]]]) &
                 is.finite(prediction_grid[[coords[1]]]) & is.finite(prediction_grid[[coords[2]]])
  n_original <- nrow(prediction_grid)
  
  if (!all(coord_valid)) {
    n_removed <- sum(!coord_valid)
    if (verbose) cat("Filtering out", n_removed, "prediction grid cells with invalid coordinates (NaN/NA/Inf) before kriging.\n")
    prediction_grid <- prediction_grid[coord_valid, , drop = FALSE]
  }

  # Validate and run kriging on filtered grid (nmax avoids singular covariance matrices)
  if (!identical(vario_fit, NA) && is.data.frame(vario_fit)) {
    validate_kriging_inputs(vario_fit, prediction_grid, dat_locations, value_var, coords)
  }
  krig_result <- run_ordinary_kriging(dat_locations, vario_fit, prediction_grid, value_var, coords, nmax = kriging_nmax)
  kriged <- krig_result$kriged
  
  # Expand predictions back to original grid size, with NA for filtered-out rows
  if (!all(coord_valid)) {
    kriged_pred_full <- rep(NA_real_, n_original)
    kriged_var_full <- rep(NA_real_, n_original)
    kriged_pred_full[coord_valid] <- krig_result$pred
    kriged_var_full[coord_valid] <- krig_result$var
    kriged_pred <- kriged_pred_full
    kriged_var <- kriged_var_full
    # Restore original grid for output
    prediction_grid <- prediction_grid_original
  } else {
    kriged_pred <- krig_result$pred
    kriged_var <- krig_result$var
  }

  # Diagnostics
  if (verbose && !is.null(kriged_pred)) {
    pred_all_na <- all(is.na(kriged_pred))
    n_na <- sum(is.na(kriged_pred))
    n_tot <- length(kriged_pred)
    if (pred_all_na) warning("Kriging produced all NaN (", n_na, "/", n_tot, ").")
    else if (n_na > 0) cat("Kriging:", n_na, "/", n_tot, "predictions are NA.\n")
  }

  kriged_sum <- sum(kriged_pred, na.rm = TRUE)

  # Predictions at data locations
  dat_locations_with_pred <- dat_locations
  if (!is.null(kriged_pred) && !is.null(kriged_var)) {
    at_pts <- extract_predictions_at_locations(kriged_pred, kriged_var, prediction_grid, dat_locations, coords)
    dat_locations_with_pred$predicted <- at_pts$pred
    # Handle negative/NaN variance: clamp to 0 before sqrt (kriging variance can be negative due to numerical issues)
    var_safe <- pmax(at_pts$var, 0, na.rm = FALSE)
    var_safe[is.na(var_safe) | !is.finite(var_safe)] <- NA_real_
    dat_locations_with_pred$predicted_error <- sqrt(var_safe)
  } else {
    dat_locations_with_pred$predicted <- NA_real_
    dat_locations_with_pred$predicted_error <- NA_real_
  }

  # Plots (variogram + kriged map + scatter)
  plots <- list()
  if (make_plots) {
    vario <- fit_result$vario
    p_vario <- ggplot2::ggplot(vario, ggplot2::aes(x = dist, y = gamma)) +
      ggplot2::geom_point(ggplot2::aes(size = np)) +
      ggplot2::scale_size_continuous(name = "Number of point pairs", range = c(1, 10)) +
      ggplot2::labs(x = "Average bin distance (degrees)", y = "Semi-variance", title = "Variogram of carbon density")
    plots$variogram <- p_vario
    if (save_plots) ggplot2::ggsave(paste0(plot_prefix, "variogram.png"), p_vario, width = 10, height = 10)

    if (!identical(vario_fit, NA) && is.data.frame(vario_fit)) {
      fit_line <- gstat::variogramLine(vario_fit, dist_vector = seq(0, max(vario$dist, na.rm = TRUE), length.out = 200))
      p_fit <- ggplot2::ggplot(vario, ggplot2::aes(x = dist, y = gamma)) +
        ggplot2::geom_point(ggplot2::aes(size = np)) +
        ggplot2::scale_size_continuous(name = "Number of point pairs", range = c(1, 10)) +
        ggplot2::geom_line(data = fit_line, ggplot2::aes(x = dist, y = gamma), color = "red", linewidth = 1) +
        ggplot2::labs(x = "Average bin distance (degrees)", y = "Semi-variance",
                      title = paste0("Fitted variogram (best: ", fit_result$best_meta$model, ")"))
      plots$variogram_fit <- p_fit
      if (save_plots) ggplot2::ggsave(paste0(plot_prefix, "fitted_variogram.png"), p_fit, width = 10, height = 10)
    }
  }
  if (make_plots && !is.null(kriged_pred)) {
    # Ensure predictions are finite (replace NaN/Inf with NA for ggplot2)
    kriged_pred_plot <- kriged_pred
    kriged_pred_plot[!is.finite(kriged_pred_plot)] <- NA_real_
    
    # Create full data frame with all predictions (NaN/Inf converted to NA)
    df_map <- data.frame(
      prediction_grid[[coords[1]]],
      prediction_grid[[coords[2]]],
      kriged_pred_plot,
      stringsAsFactors = FALSE
    )
    names(df_map) <- c(coords[1], coords[2], value_var)
    world <- ggplot2::map_data("world")
    
    # Compute plot limits from prediction grid coordinates (always valid and finite)
    x_vals <- prediction_grid[[coords[1]]][is.finite(prediction_grid[[coords[1]]])]
    y_vals <- prediction_grid[[coords[2]]][is.finite(prediction_grid[[coords[2]]])]
    if (length(x_vals) > 0 && length(y_vals) > 0) {
      xlim <- range(x_vals, na.rm = TRUE)
      ylim <- range(y_vals, na.rm = TRUE)
    } else {
      # Fallback: use data locations
      x_vals <- dat_locations[[coords[1]]][is.finite(dat_locations[[coords[1]]])]
      y_vals <- dat_locations[[coords[2]]][is.finite(dat_locations[[coords[2]]])]
      xlim <- if (length(x_vals) > 0) range(x_vals, na.rm = TRUE) else c(-180, 180)
      ylim <- if (length(y_vals) > 0) range(y_vals, na.rm = TRUE) else c(-90, 90)
    }
    # Final safety check: ensure coordinate limits are finite and have positive width
    if (!all(is.finite(c(xlim, ylim)))) {f
      xlim <- c(-180, 180)
      ylim <- c(-90, 90)
    }
    if (diff(xlim) <= 0) xlim <- c(xlim[1], xlim[1] + 1)
    if (diff(ylim) <= 0) ylim <- c(ylim[1], ylim[1] + 1)
    # Fill scale limits: use range of observed data so scale is always finite (avoids viewport error when most predictions are NA)
    obs_vals <- dat_locations[[value_var]][is.finite(dat_locations[[value_var]])]
    fill_lim <- if (length(obs_vals) > 0 && all(is.finite(obs_vals))) range(obs_vals, na.rm = TRUE) else c(0, 1)
    if (!all(is.finite(fill_lim))) fill_lim <- c(0, 1)
    
    # Filter dat_locations_with_pred to valid coordinates and finite predictions for plotting
    locs_valid <- is.finite(dat_locations_with_pred[[coords[1]]]) & 
                  is.finite(dat_locations_with_pred[[coords[2]]]) &
                  is.finite(dat_locations_with_pred$predicted)
    dat_locs_plot <- if (all(locs_valid)) dat_locations_with_pred else dat_locations_with_pred[locs_valid, ]
    
    # Ensure predicted_error is finite for plotting
    if ("predicted_error" %in% names(dat_locs_plot)) {
      dat_locs_plot$predicted_error[!is.finite(dat_locs_plot$predicted_error)] <- NA_real_
    }
    
    p_spatial <- ggplot2::ggplot(df_map, ggplot2::aes(x = .data[[coords[1]]], y = .data[[coords[2]]], fill = .data[[value_var]])) +
      ggplot2::geom_raster() +
      ggplot2::geom_polygon(data = world, ggplot2::aes(x = long, y = lat, group = group), fill = "grey40", color = NA, alpha = 0.1) +
      ggplot2::geom_point(data = dat_locs_plot, ggplot2::aes(x = .data[[coords[1]]], y = .data[[coords[2]]], fill = .data[[value_var]]), shape = 21, color = "black", stroke = 0.5, size = 2) +
      ggplot2::scale_fill_viridis_c(option = "viridis", name = "Carbon density", na.value = "transparent", limits = fill_lim) +
      ggplot2::coord_cartesian(xlim = xlim, ylim = ylim) +
      ggplot2::labs(x = "Longitude", y = "Latitude", title = "Kriged predictions")
    plots$spatial <- p_spatial
    if (save_plots) ggplot2::ggsave(paste0(plot_prefix, "kriged_predictions_spatial.png"), p_spatial, width = 10, height = 10)

    # Uncertainty map (prediction SE): same structure as predictions (raster + land + points)
    # Handle negative/NaN variance: clamp to 0 before sqrt (kriging variance can be negative due to numerical issues)
    var_safe <- pmax(kriged_var, 0, na.rm = FALSE)
    var_safe[!is.finite(var_safe)] <- NA_real_
    prediction_se <- sqrt(var_safe)
    # Ensure SE is finite (replace NaN/Inf with NA)
    prediction_se[!is.finite(prediction_se)] <- NA_real_
    
    # Create data frame with all SE values (NaN/Inf converted to NA)
    df_unc <- data.frame(
      prediction_grid[[coords[1]]],
      prediction_grid[[coords[2]]],
      prediction_se = prediction_se,
      stringsAsFactors = FALSE
    )
    names(df_unc)[1:2] <- coords
    # Filter dat_locs_plot for uncertainty plot (need finite predicted_error)
    if ("predicted_error" %in% names(dat_locs_plot)) {
      unc_locs_valid <- is.finite(dat_locs_plot[[coords[1]]]) & 
                        is.finite(dat_locs_plot[[coords[2]]]) &
                        is.finite(dat_locs_plot$predicted_error)
      dat_locs_unc <- if (all(unc_locs_valid)) dat_locs_plot else dat_locs_plot[unc_locs_valid, ]
    } else {
      dat_locs_unc <- dat_locs_plot[is.finite(dat_locs_plot[[coords[1]]]) & is.finite(dat_locs_plot[[coords[2]]]), ]
    }
    
    # SE scale limits: finite range so viewport stays valid when most SE are NA
    se_vals <- prediction_se[is.finite(prediction_se)]
    se_lim <- if (length(se_vals) > 0) range(se_vals, na.rm = TRUE) else c(0, 1)
    if (!all(is.finite(se_lim))) se_lim <- c(0, 1)
    
    p_uncertainty <- ggplot2::ggplot(df_unc, ggplot2::aes(x = .data[[coords[1]]], y = .data[[coords[2]]], fill = prediction_se)) +
      ggplot2::geom_raster() +
      ggplot2::geom_polygon(data = world, ggplot2::aes(x = long, y = lat, group = group), fill = "grey40", color = NA, alpha = 0.1) +
      {if ("predicted_error" %in% names(dat_locs_unc) && nrow(dat_locs_unc) > 0) 
        ggplot2::geom_point(data = dat_locs_unc, ggplot2::aes(x = .data[[coords[1]]], y = .data[[coords[2]]], fill = predicted_error), shape = 21, color = "black", stroke = 0.5, size = 2)} +
      ggplot2::scale_fill_viridis_c(option = "plasma", name = "Prediction SE", na.value = "transparent", limits = se_lim) +
      ggplot2::coord_cartesian(xlim = xlim, ylim = ylim) +
      ggplot2::labs(x = "Longitude", y = "Latitude", title = "Kriged prediction uncertainty (SE)")
    plots$spatial_uncertainty <- p_uncertainty
    if (save_plots) ggplot2::ggsave(paste0(plot_prefix, "kriged_uncertainty_spatial.png"), p_uncertainty, width = 10, height = 10)

    # Scatter plot: use only rows with valid observed and predicted to avoid "Removed N rows" warning
    scatter_use <- stats::complete.cases(dat_locations_with_pred[[value_var]], dat_locations_with_pred$predicted)
    dat_scatter <- dat_locations_with_pred[scatter_use, ]
    n_omit <- nrow(dat_locations_with_pred) - nrow(dat_scatter)
    
    # Performance label for scatter (R², adj R², RMSE, MAE)
    obs <- dat_scatter[[value_var]]
    pred <- dat_scatter$predicted
    n_use <- length(obs)
    if (n_use > 0) {
      r2 <- stats::cor(obs, pred)^2
      adj_r2 <- 1 - (1 - r2) * (n_use - 1) / max(n_use - 2, 1)
      rmse <- sqrt(mean((obs - pred)^2))
      mae <- mean(abs(obs - pred))
      perf_label <- sprintf("R² = %.3f  adj R² = %.3f  RMSE = %.4f  MAE = %.4f  n = %d", r2, adj_r2, rmse, mae, n_use)
      if (n_omit > 0) perf_label <- paste0(perf_label, " (", n_omit, " pts without prediction omitted)")
    } else {
      perf_label <- if (n_omit > 0) paste0("No valid pairs (all ", n_omit, " predictions missing)") else "No valid pairs for metrics"
    }
    p_scatter <- ggplot2::ggplot(dat_scatter, ggplot2::aes(x = .data[[value_var]], y = predicted)) +
      ggplot2::geom_point(alpha = 0.7, color = "darkblue") +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = predicted - predicted_error, ymax = predicted + predicted_error), width = 0.001, alpha = 0.1) +
      ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
      ggplot2::annotate("label", x = -Inf, y = Inf, label = perf_label, hjust = -0.05, vjust = 1.1, size = 3.5, fill = "white", alpha = 0.9) +
      ggplot2::labs(x = "Observed", y = "Kriged predicted", title = "Observed vs Kriged (at data locations)") +
      ggplot2::theme_minimal()
    plots$scatter <- p_scatter
    if (save_plots) ggplot2::ggsave(paste0(plot_prefix, "kriged_predictions_scatter.png"), p_scatter, width = 8, height = 8)
  }

  list(
    dat_locations = dat_locations,
    vario = vario,
    fit_result = fit_result,
    prediction_grid = prediction_grid,
    kriged = kriged,
    kriged_pred = kriged_pred,
    kriged_var = kriged_var,
    dat_locations_with_pred = dat_locations_with_pred,
    kriged_sum = kriged_sum,
    plots = plots
  )
}


# ---- Prediction intervals and uncertainty quantification ----

#' Compute prediction intervals for GAM predictions.
#'
#' @param gam_fit Fitted GAM model from mgcv::gam.
#' @param newdata Data frame with predictor variables for prediction.
#' @param level Confidence level (default 0.95 for 95% intervals).
#' @param type Type of interval: "response" (on response scale) or "link" (on link scale).
#' @return List with fit (predictions), se.fit (standard errors), lower, upper (interval bounds).
compute_gam_prediction_intervals <- function(gam_fit, newdata, level = 0.95, type = "response") {
  if (!requireNamespace("mgcv", quietly = TRUE)) {
    stop("Package mgcv is required.")
  }
  
  # Get predictions with standard errors
  pred <- try(mgcv::predict.gam(gam_fit, newdata = newdata, type = type, se.fit = TRUE), silent = FALSE)
  if (inherits(pred, "try-error")) {
    return(list(fit = NULL, se.fit = NULL, lower = NULL, upper = NULL, error = attr(pred, "condition")$message))
  }
  
  fit <- as.numeric(pred$fit)
  se <- as.numeric(pred$se.fit)
  
  # Compute critical value (z-score for normal approximation)
  alpha <- 1 - level
  z_crit <- stats::qnorm(1 - alpha / 2)
  
  # Compute intervals
  lower <- fit - z_crit * se
  upper <- fit + z_crit * se
  
  # If using log link (Gamma family), ensure positive bounds
  if (inherits(gam_fit$family, "family") && gam_fit$family$link == "log") {
    lower <- pmax(lower, 0)  # Ensure non-negative
  }
  
  list(
    fit = fit,
    se.fit = se,
    lower = lower,
    upper = upper,
    level = level
  )
}

#' Compute prediction intervals for kriging predictions.
#'
#' @param kriged_pred Vector of kriged predictions.
#' @param kriged_var Vector of kriging variances.
#' @param level Confidence level (default 0.95 for 95% intervals).
#' @return List with pred, se, lower, upper (interval bounds).
compute_kriging_prediction_intervals <- function(kriged_pred, kriged_var, level = 0.95) {
  if (is.null(kriged_pred) || is.null(kriged_var)) {
    return(list(pred = NULL, se = NULL, lower = NULL, upper = NULL))
  }
  
  # Ensure variance is non-negative
  var_safe <- pmax(kriged_var, 0, na.rm = FALSE)
  se <- sqrt(var_safe)
  
  # Compute critical value
  alpha <- 1 - level
  z_crit <- stats::qnorm(1 - alpha / 2)
  
  # Compute intervals
  lower <- kriged_pred - z_crit * se
  upper <- kriged_pred + z_crit * se
  
  list(
    pred = kriged_pred,
    se = se,
    lower = lower,
    upper = upper,
    level = level
  )
}

#' Create spatial uncertainty map data frame with prediction intervals.
#'
#' @param prediction_grid Data frame with coordinate columns.
#' @param predictions Vector of predictions.
#' @param se Vector of standard errors (or NULL if using intervals directly).
#' @param lower Vector of lower interval bounds (optional, computed from se if NULL).
#' @param upper Vector of upper interval bounds (optional, computed from se if NULL).
#' @param coords Character vector of coordinate names.
#' @param level Confidence level (used if computing intervals from se).
#' @return Data frame with coordinates, predictions, se, lower, upper, and uncertainty metrics.
create_uncertainty_map_data <- function(prediction_grid, predictions, se = NULL, 
                                       lower = NULL, upper = NULL,
                                       coords = c("longitude", "latitude"),
                                       level = 0.95) {
  df <- data.frame(
    prediction_grid[, coords, drop = FALSE],
    prediction = predictions,
    stringsAsFactors = FALSE
  )
  
  # Compute intervals if not provided
  if (is.null(lower) || is.null(upper)) {
    if (!is.null(se)) {
      alpha <- 1 - level
      z_crit <- stats::qnorm(1 - alpha / 2)
      df$se <- se
      df$lower <- predictions - z_crit * se
      df$upper <- predictions + z_crit * se
    } else {
      df$se <- NA_real_
      df$lower <- NA_real_
      df$upper <- NA_real_
    }
  } else {
    df$lower <- lower
    df$upper <- upper
    if (!is.null(se)) {
      df$se <- se
    } else {
      # Approximate se from interval width
      df$se <- (upper - lower) / (2 * stats::qnorm(1 - (1 - level) / 2))
    }
  }
  
  # Compute additional uncertainty metrics
  df$interval_width <- df$upper - df$lower
  df$cv <- df$se / pmax(abs(df$prediction), 1e-10)  # Coefficient of variation
  
  df
}

# ---- Applicability domain analysis ----

#' Compute environmental similarity scores between prediction locations and training data.
#'
#' @param training_data Data frame with training observations and predictor variables.
#' @param prediction_data Data frame with prediction locations and predictor variables.
#' @param predictor_vars Character vector of predictor variable names.
#' @param method Similarity metric: "euclidean" (scaled), "mahalanobis", or "range_overlap".
#' @param verbose If TRUE, print diagnostic information.
#' @return List with similarity scores (vector), flag_outside_range (logical vector), and summary statistics.
compute_environmental_similarity <- function(training_data, prediction_data, predictor_vars,
                                            method = "euclidean", verbose = FALSE) {
  # Get common predictors
  pv <- intersect(predictor_vars, names(training_data))
  pv <- intersect(pv, names(prediction_data))
  
  if (length(pv) == 0) {
    warning("No common predictor variables found")
    return(list(similarity_scores = NULL, flag_outside_range = NULL, summary = NULL))
  }
  
  # Extract predictor matrices
  X_train <- as.matrix(training_data[, pv, drop = FALSE])
  X_pred <- as.matrix(prediction_data[, pv, drop = FALSE])
  
  # Check for non-finite values and replace with NA
  X_train[!is.finite(X_train)] <- NA_real_
  X_pred[!is.finite(X_pred)] <- NA_real_
  
  # For training data, we need complete cases to compute statistics
  train_complete <- stats::complete.cases(X_train)
  n_train_complete <- sum(train_complete)
  
  if (verbose) {
    cat("Environmental similarity computation:\n")
    cat("  Training data: ", nrow(training_data), " rows, ", n_train_complete, " complete cases\n", sep = "")
    cat("  Prediction data: ", nrow(prediction_data), " rows\n", sep = "")
    cat("  Predictor variables: ", length(pv), " (", paste(pv, collapse = ", "), ")\n", sep = "")
    
    # Check missing data patterns in prediction grid
    pred_missing <- colSums(is.na(X_pred))
    cat("  Missing values in prediction grid per variable:\n")
    for (i in seq_along(pv)) {
      cat("    ", pv[i], ": ", pred_missing[i], " (", round(100 * pred_missing[i] / nrow(X_pred), 1), "%)\n", sep = "")
    }
  }
  
  if (n_train_complete == 0) {
    warning("No complete cases in training data for similarity computation")
    return(list(similarity_scores = rep(NA_real_, nrow(prediction_data)), 
                flag_outside_range = rep(FALSE, nrow(prediction_data)), 
                summary = NULL))
  }
  
  # Use only complete training cases for computing statistics
  X_train_complete <- X_train[train_complete, , drop = FALSE]
  
  # Remove columns that are all NA in training data
  train_cols_valid <- colSums(!is.na(X_train_complete)) > 0
  if (sum(train_cols_valid) == 0) {
    warning("No valid predictor columns in training data")
    return(list(similarity_scores = rep(NA_real_, nrow(prediction_data)), 
                flag_outside_range = rep(FALSE, nrow(prediction_data)), 
                summary = NULL))
  }
  
  X_train_complete <- X_train_complete[, train_cols_valid, drop = FALSE]
  pv_valid <- pv[train_cols_valid]
  
  if (verbose) {
    cat("  Using ", length(pv_valid), " predictors with valid training data\n", sep = "")
  }
  
  # Compute similarity based on method
  if (method == "euclidean") {
    # Scale by training data standard deviations
    train_mean <- colMeans(X_train, na.rm = TRUE)
    train_sd <- apply(X_train, 2, sd, na.rm = TRUE)
    train_sd[train_sd == 0 | !is.finite(train_sd)] <- 1  # Avoid division by zero or Inf/NaN
    
    # Check for Inf/NaN in means
    train_mean[!is.finite(train_mean)] <- 0
    
    # Scale training dataset
    X_train_scaled <- scale(X_train_complete, center = train_mean, scale = train_sd)
    X_train_scaled[!is.finite(X_train_scaled)] <- 0
    
    # Compute minimum distance from each prediction point to training data
    # Handle missing values by computing distance only on available predictors per point
    similarity_scores <- apply(X_pred, 1, function(x) {
      # Find which predictors are available (finite) for this prediction point
      x_valid <- is.finite(x)
      
      if (sum(x_valid) == 0) {
        # No valid predictors for this point
        return(NA_real_)
      }
      
      # Use only available predictors for this point
      x_available <- x[x_valid]
      train_mean_available <- train_mean[x_valid]
      train_sd_available <- train_sd[x_valid]
      
      # Scale this point using available predictors only
      x_scaled <- (x_available - train_mean_available) / train_sd_available
      x_scaled[!is.finite(x_scaled)] <- 0
      
      # Compute distances to training data using only the available predictors
      train_scaled_available <- X_train_scaled[, x_valid, drop = FALSE]
      diff_mat <- train_scaled_available - matrix(x_scaled, nrow = nrow(train_scaled_available), 
                                                   ncol = length(x_scaled), byrow = TRUE)
      diff_mat[!is.finite(diff_mat)] <- 0
      dists <- sqrt(rowSums(diff_mat^2))
      dists <- dists[is.finite(dists)]
      
      if (length(dists) == 0) {
        return(NA_real_)
      }
      
      # Return minimum distance (will normalize later)
      min(dists, na.rm = TRUE)
    })
    
    # Remove NA values before computing max_dist
    valid_scores <- similarity_scores[is.finite(similarity_scores)]
    if (verbose) {
      cat("  Euclidean distances computed: ", length(valid_scores), " finite out of ", length(similarity_scores), "\n", sep = "")
    }
    if (length(valid_scores) == 0) {
      warning("All similarity scores are non-finite after distance computation")
      similarity_scores[] <- NA_real_
    } else {
      # Lower score = more similar (closer to training data)
      # Convert to similarity (higher = more similar)
      max_dist <- max(valid_scores, na.rm = TRUE)
      if (verbose) {
        cat("  Max distance: ", max_dist, "\n", sep = "")
      }
      if (!is.finite(max_dist) || max_dist <= 0) {
        # If max_dist is not valid, set all to NA or use a default
        warning("max_dist is not finite or <= 0, setting similarity scores to NA")
        similarity_scores[] <- NA_real_
      } else {
        # Only convert finite scores
        finite_idx <- is.finite(similarity_scores)
        similarity_scores[finite_idx] <- 1 - (similarity_scores[finite_idx] / (max_dist + 1e-10))
        # Ensure similarity scores are in [0, 1]
        similarity_scores[finite_idx] <- pmax(0, pmin(1, similarity_scores[finite_idx]))
        if (verbose) {
          cat("  Similarity scores range: ", min(similarity_scores[finite_idx], na.rm = TRUE), 
              " to ", max(similarity_scores[finite_idx], na.rm = TRUE), "\n", sep = "")
        }
      }
    }
    
  } else if (method == "mahalanobis") {
    # Compute Mahalanobis distance (using only valid predictors)
    train_cov <- cov(X_train_complete, use = "complete.obs")
    train_mean <- colMeans(X_train_complete, na.rm = TRUE)
    
    # Check for non-finite values
    train_mean[!is.finite(train_mean)] <- 0
    train_cov[!is.finite(train_cov)] <- 0
    
    # Check if covariance matrix is invertible
    if (det(train_cov) < 1e-10 || !is.finite(det(train_cov))) {
      warning("Covariance matrix near-singular or non-finite, using euclidean method instead")
      return(compute_environmental_similarity(training_data, prediction_data, predictor_vars, "euclidean", verbose))
    }
    
    # Compute Mahalanobis distance handling missing values per prediction point
    similarity_scores <- apply(X_pred, 1, function(x) {
      # Find which predictors are available
      x_valid <- is.finite(x)
      if (sum(x_valid) == 0) {
        return(NA_real_)
      }
      
      # Use only available predictors
      x_available <- x[x_valid]
      train_mean_available <- train_mean[x_valid]
      train_cov_available <- train_cov[x_valid, x_valid, drop = FALSE]
      
      diff <- x_available - train_mean_available
      if (any(!is.finite(diff))) {
        return(NA_real_)
      }
      tryCatch({
        dist_val <- sqrt(t(diff) %*% solve(train_cov_available) %*% diff)
        if (is.finite(dist_val)) dist_val else NA_real_
      }, error = function(e) NA_real_)
    })
    
    # Convert distance to similarity
    valid_scores <- similarity_scores[is.finite(similarity_scores)]
    if (length(valid_scores) == 0) {
      warning("All Mahalanobis distances are non-finite")
      similarity_scores[] <- NA_real_
    } else {
      max_dist <- max(valid_scores, na.rm = TRUE)
      if (!is.finite(max_dist) || max_dist <= 0) {
        warning("max_dist is not finite or <= 0 for Mahalanobis, setting similarity scores to NA")
        similarity_scores[] <- NA_real_
      } else {
        finite_idx <- is.finite(similarity_scores)
        similarity_scores[finite_idx] <- 1 - (similarity_scores[finite_idx] / (max_dist + 1e-10))
        similarity_scores[finite_idx] <- pmax(0, pmin(1, similarity_scores[finite_idx]))
      }
    }
    
  } else if (method == "range_overlap") {
    # Compute fraction of predictors within training range (using only valid predictors)
    train_ranges <- apply(X_train_complete, 2, range, na.rm = TRUE)
    
    # Check for non-finite ranges
    train_ranges[!is.finite(train_ranges)] <- 0
    
    similarity_scores <- apply(X_pred, 1, function(x) {
      # Find which predictors are available
      x_valid <- is.finite(x)
      if (sum(x_valid) == 0) {
        return(NA_real_)
      }
      
      # Use only available predictors
      x_available <- x[x_valid]
      train_ranges_available <- train_ranges[, x_valid, drop = FALSE]
      
      in_range <- (x_available >= train_ranges_available[1, ]) & (x_available <= train_ranges_available[2, ])
      mean(in_range, na.rm = TRUE)  # Fraction of available predictors in range
    })
    
    # Ensure scores are finite
    similarity_scores[!is.finite(similarity_scores)] <- NA_real_
  } else {
    stop("Unknown method: ", method)
  }
  
  # Ensure similarity_scores has correct length (should match nrow(prediction_data))
  if (length(similarity_scores) != nrow(prediction_data)) {
    warning("Length mismatch: similarity_scores (", length(similarity_scores), 
            ") != nrow(prediction_data) (", nrow(prediction_data), ")")
    # This shouldn't happen, but handle it gracefully
    if (length(similarity_scores) > nrow(prediction_data)) {
      similarity_scores <- similarity_scores[1:nrow(prediction_data)]
    } else {
      similarity_scores <- c(similarity_scores, rep(NA_real_, nrow(prediction_data) - length(similarity_scores)))
    }
  }
  
  # Flag predictions outside training range
  flag_outside_range <- rep(FALSE, nrow(prediction_data))
  finite_scores <- is.finite(similarity_scores)
  if (any(finite_scores)) {
    flag_outside_range[finite_scores] <- similarity_scores[finite_scores] < 0.5
  }
  
  # Compute summary statistics (only on finite scores)
  finite_scores <- is.finite(similarity_scores)
  n_finite <- sum(finite_scores)
  
  if (n_finite > 0) {
    summary_stats <- list(
      mean_similarity = mean(similarity_scores[finite_scores], na.rm = TRUE),
      min_similarity = min(similarity_scores[finite_scores], na.rm = TRUE),
      max_similarity = max(similarity_scores[finite_scores], na.rm = TRUE),
      n_outside_range = sum(flag_outside_range, na.rm = TRUE),
      pct_outside_range = 100 * mean(flag_outside_range, na.rm = TRUE),
      n_finite_scores = n_finite,
      n_total_pred = nrow(prediction_data),
      pct_finite = 100 * n_finite / nrow(prediction_data)
    )
  } else {
    summary_stats <- list(
      mean_similarity = NA_real_,
      min_similarity = NA_real_,
      max_similarity = NA_real_,
      n_outside_range = 0,
      pct_outside_range = 0,
      n_finite_scores = 0,
      n_total_pred = nrow(prediction_data),
      pct_finite = 0
    )
  }
  
  # similarity_scores already has length matching prediction_data
  similarity_full <- similarity_scores
  
  list(
    similarity_scores = similarity_full,
    flag_outside_range = flag_outside_range,
    summary = summary_stats,
    method = method
  )
}

#' Flag predictions where covariates are outside training range (values are more extreme than any in training data)
#'
#' @param training_data Data frame with training observations.
#' @param prediction_data Data frame with prediction locations.
#' @param predictor_vars Character vector of predictor variable names.
#' @param strict If TRUE, flag if ANY predictor is outside range; if FALSE, flag if >50% are outside.
#' @return Logical vector indicating which predictions are outside applicability domain.
flag_outside_applicability_domain <- function(training_data, prediction_data, predictor_vars,
                                             strict = FALSE) {
  pv <- intersect(predictor_vars, names(training_data))
  pv <- intersect(pv, names(prediction_data))
  
  if (length(pv) == 0) {
    return(rep(FALSE, nrow(prediction_data)))
  }
  
  # Compute training ranges
  train_ranges <- lapply(pv, function(v) {
    range(training_data[[v]], na.rm = TRUE)
  })
  names(train_ranges) <- pv
  
  # Check each prediction location
  flags <- logical(nrow(prediction_data))
  
  for (i in seq_len(nrow(prediction_data))) {
    outside_count <- 0
    for (v in pv) {
      pred_val <- prediction_data[[v]][i]
      if (!is.na(pred_val)) {
        if (pred_val < train_ranges[[v]][1] || pred_val > train_ranges[[v]][2]) {
          outside_count <- outside_count + 1
        }
      }
    }
    
    if (strict) {
      flags[i] <- outside_count > 0  # Any predictor outside range
    } else {
      flags[i] <- outside_count > length(pv) * 0.5  # >50% outside range
    }
  }
  
  flags
}

#' Create applicability domain map.
#'
#' @param prediction_grid Data frame with coordinates.
#' @param similarity_scores Vector of similarity scores (or NULL).
#' @param flag_outside Logical vector indicating outside-domain locations.
#' @param coords Character vector of coordinate names.
#' @param world Optional world map data.
#' @param use_raster If TRUE, use geom_raster(); if FALSE, use geom_point() (default FALSE to avoid NaN issues with irregular grids).
#' @param point_size Size of points when use_raster = FALSE (default 0.3).
#' @param xlim Optional x-axis limits.
#' @param ylim Optional y-axis limits.
#' @return ggplot object with applicability domain map.
plot_applicability_domain <- function(prediction_grid, similarity_scores = NULL, plot_data_points = FALSE, dat = NULL,
                                     flag_outside = NULL, 
                                     coords = c("longitude", "latitude"),
                                     world = NULL,
                                     use_raster = FALSE,
                                     point_size = 0.3,
                                     xlim = NULL,
                                     ylim = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package ggplot2 is required.")
  }
  
  # Prepare data
  df <- data.frame(
    prediction_grid[, coords, drop = FALSE],
    stringsAsFactors = FALSE
  )
  
  if (!is.null(similarity_scores)) {
    df$similarity <- similarity_scores
  }
  
  if (!is.null(flag_outside)) {
    df$outside_domain <- flag_outside
  }
  
  # Get world map if not provided
  if (is.null(world)) {
    if (requireNamespace("maps", quietly = TRUE)) {
      world <- ggplot2::map_data("world")
    }
  }
  
  # Filter out NA/NaN/Inf values before plotting
  if (!is.null(similarity_scores)) {
    df <- df[is.finite(df$similarity), , drop = FALSE]
  }
  if (!is.null(flag_outside)) {
    df <- df[!is.na(df$outside_domain), , drop = FALSE]
  }
  
  # Create base plot
  p <- ggplot2::ggplot()

  
  # Create plot with similarity scores or flags
  if (!is.null(similarity_scores) && nrow(df) > 0) {
    if (use_raster) {
      p <- p + ggplot2::geom_raster(
        data = df,
        ggplot2::aes(x = .data[[coords[1]]], y = .data[[coords[2]]], fill = similarity),
        na.rm = TRUE
      ) +
      ggplot2::scale_fill_viridis_c(
        option = "plasma",
        name = "Similarity",
        na.value = "transparent"
      )
    } else {
      p <- p + ggplot2::geom_point(
        data = df,
        ggplot2::aes(x = .data[[coords[1]]], y = .data[[coords[2]]], color = similarity),
        size = point_size,
        alpha = 0.7,
        na.rm = TRUE
      ) +
      ggplot2::scale_color_viridis_c(
        option = "plasma",
        name = "Similarity",
        na.value = "transparent"
      )
    }
    p <- p + ggplot2::labs(
      x = "Longitude",
      y = "Latitude",
      title = "Environmental similarity to training data"
    )
  } else if (!is.null(flag_outside) && nrow(df) > 0) {
    if (use_raster) {
      p <- p + ggplot2::geom_raster(
        data = df,
        ggplot2::aes(x = .data[[coords[1]]], y = .data[[coords[2]]], fill = outside_domain),
        na.rm = TRUE
      ) +
      ggplot2::scale_fill_manual(
        values = c("FALSE" = "green", "TRUE" = "red"),
        name = "Outside\nDomain",
        na.value = "transparent"
      )
    } else {
      p <- p + ggplot2::geom_point(
        data = df,
        ggplot2::aes(x = .data[[coords[1]]], y = .data[[coords[2]]], color = outside_domain),
        size = point_size,
        alpha = 0.7,
        na.rm = TRUE
      ) +
      ggplot2::scale_color_manual(
        values = c("FALSE" = "green", "TRUE" = "red"),
        name = "Outside\nDomain",
        na.value = "transparent"
      )
    }
    p <- p + ggplot2::labs(
      x = "Longitude",
      y = "Latitude",
      title = "Applicability domain (red = outside training range)"
    )
  } else {
    stop("Must provide either similarity_scores or flag_outside, and data must have valid values")
  }

  # Add world map first (so it's behind the data)
  if (!is.null(world)) {
    p <- p + ggplot2::geom_polygon(
      data = world,
      ggplot2::aes(x = long, y = lat, group = group),
      fill = "grey90", color = "grey70", alpha = 0.3, inherit.aes = FALSE
    )
  }


  if (plot_data_points && !is.null(dat)) {
    p <- p + ggplot2::geom_point(
      data = dat,
      ggplot2::aes(x = .data[[coords[1]]], y = .data[[coords[2]]]),
      size = point_size,
      alpha = 0.7,
      na.rm = TRUE
    )
  }

  # Add coordinate limits if provided
  if (!is.null(xlim) && !is.null(ylim)) {
    p <- p + ggplot2::coord_cartesian(xlim = xlim, ylim = ylim)
  } else if (!is.null(xlim)) {
    p <- p + ggplot2::coord_cartesian(xlim = xlim)
  } else if (!is.null(ylim)) {
    p <- p + ggplot2::coord_cartesian(ylim = ylim)
  }
  
  p <- p + ggplot2::theme_minimal()
  
  p
}

# ---- Regression kriging (trend + kriged residuals) ----

#' Aggregate to one row per location with median residuals and median value (for regression kriging).
#'
#' @param dat Data frame with coords, "residuals", and value_var columns.
#' @param coords Character vector of coordinate names.
#' @param value_var Name of response variable (kept for observed vs predicted).
#' @return Data frame with one row per unique location: coords, value_var, residuals.
aggregate_residuals_to_locations <- function(dat,
                                             coords = c("longitude", "latitude"),
                                             value_var = "carbon_density_g_c_cm3") {
  dat %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(coords))) %>%
    dplyr::summarise(
      !!value_var := median(.data[[value_var]], na.rm = TRUE),
      residuals = median(.data$residuals, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::filter(!is.na(.data$residuals))
}


#' Run kriging on residuals (residuals ~ 1).
#'
#' @param dat_resid Data frame with coords and "residuals" column.
#' @param vario_fit Fitted vgm from fit_variogram_best.
#' @param prediction_grid Data frame with coordinate columns.
#' @param coords Character vector of coordinate names.
#' @return List with kriged (result), pred (vector), var (vector), success (logical).
run_residual_kriging <- function(dat_resid, vario_fit, prediction_grid,
                                 coords = c("longitude", "latitude")) {
  if (identical(vario_fit, NA) || !is.data.frame(vario_fit)) {
    return(list(kriged = NULL, pred = NULL, var = NULL, success = FALSE))
  }
  locations <- as.formula(paste("~", paste(coords, collapse = " + ")))
  kriged <- tryCatch(
    gstat::krige(residuals ~ 1, locations = locations, data = dat_resid,
                 newdata = prediction_grid, model = vario_fit),
    error = function(e) {
      warning("Residual kriging failed: ", conditionMessage(e))
      return(NULL)
    }
  )
  if (is.null(kriged)) return(list(kriged = NULL, pred = NULL, var = NULL, success = FALSE))
  if (is.data.frame(kriged)) {
    pred <- kriged$var1.pred
    var <- kriged$var1.var
  } else if (!is.null(methods::slotNames(kriged)) && "data" %in% methods::slotNames(kriged)) {
    pred <- kriged@data$var1.pred
    var <- kriged@data$var1.var
  } else {
    return(list(kriged = kriged, pred = NULL, var = NULL, success = FALSE))
  }
  list(kriged = kriged, pred = pred, var = var, success = TRUE)
}


#' Predict trend on grid; if prediction_grid lacks formula covariates, use constant fallback.
trend_at_grid <- function(lm_fit, prediction_grid, fallback) {
  tryCatch({
    p <- stats::predict(lm_fit, newdata = prediction_grid)
    if (length(p) != nrow(prediction_grid) || any(is.na(p)))
      rep(fallback, nrow(prediction_grid))
    else
      as.numeric(p)
  }, error = function(e) rep(fallback, nrow(prediction_grid)))
}


#' Fill prediction_grid with covariate columns from nearest data location (for regression kriging trend on grid).
#' Requires FNN. Modifies prediction_grid in place by adding missing pred_vars; returns prediction_grid.
fill_prediction_grid_covariates <- function(prediction_grid, dat, coords, pred_vars) {
  pred_vars <- intersect(pred_vars, names(dat))
  if (length(pred_vars) == 0) return(prediction_grid)
  missing <- setdiff(pred_vars, names(prediction_grid))
  if (length(missing) == 0) return(prediction_grid)
  if (!requireNamespace("FNN", quietly = TRUE)) return(prediction_grid)
  dat_loc <- dat %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(coords))) %>%
    dplyr::summarise(
      dplyr::across(dplyr::all_of(missing), ~ if (is.numeric(.)) median(., na.rm = TRUE) else .[1]),
      .groups = "drop"
    )
  nn <- FNN::get.knnx(
    as.matrix(dat_loc[, coords]),
    as.matrix(prediction_grid[, coords]),
    k = 1
  )$nn.index[, 1]
  for (v in missing) prediction_grid[[v]] <- dat_loc[[v]][nn]
  prediction_grid
}


#' Fill prediction_grid with smoothly interpolated covariate columns (spatial GAM per covariate).
#' Avoids blocky Voronoi patterns from nearest-neighbour; use when continuous/smooth predictions are required.
#' @param prediction_grid Data frame with coords.
#' @param dat_locations Data frame with one row per location: coords and pred_vars (numeric).
#' @param coords Character vector of coordinate names.
#' @param pred_vars Character vector of covariate names to interpolate.
#' @param k GAM basis dimension for smooth (default 30).
#' @return prediction_grid with pred_vars added/overwritten (smooth spatial interpolation).
fill_prediction_grid_covariates_smooth <- function(prediction_grid, dat_locations, coords, pred_vars, k = 100) {
  pred_vars <- intersect(pred_vars, names(dat_locations))
  if (length(pred_vars) == 0) return(prediction_grid)
  if (!requireNamespace("mgcv", quietly = TRUE)) return(fill_prediction_grid_covariates(prediction_grid, dat_locations, coords, pred_vars))
  k <- min(max(10, k), 50)
  for (v in pred_vars) {
    form <- as.formula(paste(v, "~ s(", coords[1], ",", coords[2], ", k = ", k, ")"))
    fit <- try(mgcv::gam(form, data = dat_locations, method = "REML"), silent = TRUE)
    if (inherits(fit, "try-error")) next
    p <- try(as.numeric(mgcv::predict.gam(fit, newdata = prediction_grid, type = "response")), silent = TRUE)
    if (!inherits(p, "try-error") && length(p) == nrow(prediction_grid)) prediction_grid[[v]] <- p
  }
  prediction_grid
}


#' Full regression kriging pipeline: LM trend + residual variogram + residual kriging + combine.
#'
#' Run with different subsets by passing e.g. dat or dat %>% filter(Region == "Mediterranean Sea").
#' prediction_grid should have coordinate columns; if it also has all formula RHS variables,
#' trend at grid points is from predict(lm_fit); otherwise constant trend (mean of response) is used.
#'
#' @param dat Data frame with coordinate, value, and covariate columns (can have duplicate locations).
#' @param formula Formula for the trend (e.g. carbon_density_g_c_cm3 ~ KD_closest + ...).
#' @param value_var Name of response variable (default "carbon_density_g_c_cm3").
#' @param coords Character vector of coordinate names (default c("longitude", "latitude")).
#' @param prediction_grid Grid for predictions (must have coords; optional covariates for trend).
#' @param smooth_grid_covariates If TRUE (default), interpolate covariates onto the grid with a spatial GAM for a smooth trend surface; if FALSE, use nearest-neighbour (can be blocky).
#' @param make_plots If TRUE, create variogram, spatial, and scatter plots.
#' @param save_plots If TRUE, save plots to files.
#' @param plot_prefix Path prefix for saved plots (e.g. "figures/interpolation/").
#' @param world Optional map data for land overlay (polygon); if NULL, uses map_data("world").
#' @param xlim,ylim Optional plot extent (numeric length 2); if NULL, from prediction_grid.
#' @param verbose If TRUE, print progress and variogram comparison.
#' @return List with lm_fit, dat_resid, vario_resid, fit_resid, resid_kriged, grid_pred, grid_var,
#'   prediction_grid (grid with smoothly filled covariates when used), dat_locations_with_pred (observed + predicted at locations), plots.
regression_kriging_pipeline <- function(dat,
                                       formula,
                                       value_var = "carbon_density_g_c_cm3",
                                       coords = c("longitude", "latitude"),
                                       prediction_grid = NULL,
                                       smooth_grid_covariates = TRUE,
                                       make_plots = TRUE,
                                       save_plots = TRUE,
                                       plot_prefix = "figures/interpolation/",
                                       world = NULL,
                                       xlim = NULL,
                                       ylim = NULL,
                                       verbose = TRUE) {
  if (verbose) cat("Regression kriging: fitting trend model...\n")
  lm_fit <- stats::lm(formula, data = dat)
  # residuals(lm_fit) has length = nrow(model frame); dat may have more rows (NA in formula vars)
  res <- rep(NA_real_, nrow(dat))
  omit <- stats::na.action(lm_fit)
  if (is.null(omit)) {
    res[] <- stats::residuals(lm_fit)
  } else {
    res[-as.integer(omit)] <- stats::residuals(lm_fit)
  }
  dat$residuals <- res

  dat_resid <- aggregate_residuals_to_locations(dat, coords = coords, value_var = value_var)
  if (verbose) cat("Regression kriging: using", nrow(dat_resid), "unique locations for residual variogram and kriging.\n")

  pred_vars_grid <- setdiff(all.vars(formula), value_var)
  pred_vars_grid <- intersect(pred_vars_grid, names(dat))
  if (length(pred_vars_grid) > 0 && any(!pred_vars_grid %in% names(prediction_grid))) {
    dat_loc_cov <- dat %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(coords))) %>%
      dplyr::summarise(
        dplyr::across(dplyr::all_of(pred_vars_grid), ~ if (is.numeric(.)) median(., na.rm = TRUE) else .[1]),
        .groups = "drop"
      )
    if (smooth_grid_covariates) {
      prediction_grid <- fill_prediction_grid_covariates_smooth(prediction_grid, dat_loc_cov, coords, pred_vars_grid)
      if (verbose) cat("Regression kriging: filled prediction_grid with smoothly interpolated covariates (spatial GAM).\n")
    } else {
      prediction_grid <- fill_prediction_grid_covariates(prediction_grid, dat, coords, pred_vars_grid)
      if (verbose) cat("Regression kriging: filled prediction_grid with covariates from nearest data locations.\n")
    }
  }

  locations <- as.formula(paste("~", paste(coords, collapse = " + ")))
  vario_resid <- empirical_variogram(dat_resid, residuals ~ 1, locations)
  fit_resid <- fit_variogram_best(vario_resid, verbose = verbose)
  vario_fit_resid <- fit_resid$vario_fit

  fallback_trend <- mean(dat[[value_var]], na.rm = TRUE)
  trend_grid <- trend_at_grid(lm_fit, prediction_grid, fallback_trend)

  resid_result <- run_residual_kriging(dat_resid, vario_fit_resid, prediction_grid, coords)
  resid_pred <- resid_result$pred
  resid_var <- resid_result$var

  if (resid_result$success && !is.null(resid_pred) && !is.null(resid_var)) {
    grid_pred <- trend_grid + resid_pred
    grid_var <- resid_var
  } else {
    grid_pred <- trend_grid
    grid_var <- rep(NA_real_, nrow(prediction_grid))
  }

  # Predictions at data locations: trend at locations + residual from nearest grid point.
  # dat_resid only has coords, value_var, residuals; predict() needs the formula's predictor columns.
  # Build one row per location with covariates from dat (median per location).
  pred_vars <- setdiff(all.vars(formula), value_var)
  pred_vars <- intersect(pred_vars, names(dat))
  if (length(pred_vars) > 0) {
    dat_loc_covariates <- dat %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(coords))) %>%
      dplyr::summarise(
        dplyr::across(dplyr::all_of(pred_vars), ~ if (is.numeric(.)) median(., na.rm = TRUE) else .[1]),
        .groups = "drop"
      )
    dat_resid_for_trend <- dplyr::left_join(dat_resid, dat_loc_covariates, by = coords)
  } else {
    dat_resid_for_trend <- dat_resid
  }
  trend_at_locs <- tryCatch(
    as.numeric(stats::predict(lm_fit, newdata = dat_resid_for_trend)),
    error = function(e) rep(fallback_trend, nrow(dat_resid))
  )
  if (requireNamespace("FNN", quietly = TRUE)) {
    at_pts <- extract_predictions_at_locations(resid_pred, resid_var, prediction_grid, dat_resid, coords)
    dat_resid$predicted <- trend_at_locs + at_pts$pred
    dat_resid$predicted_error <- sqrt(at_pts$var)
  } else {
    dat_resid$predicted <- trend_at_locs
    dat_resid$predicted_error <- NA_real_
  }

  plots <- list()
  if (make_plots) {
    if (!is.null(fit_resid$best_meta) && !identical(vario_fit_resid, NA) && is.data.frame(vario_fit_resid)) {
      fit_line_resid <- gstat::variogramLine(vario_fit_resid, dist_vector = seq(0, max(vario_resid$dist, na.rm = TRUE), length.out = 200))
      p_fit_resid <- ggplot2::ggplot(vario_resid, ggplot2::aes(x = dist, y = gamma)) +
        ggplot2::geom_point(ggplot2::aes(size = np)) +
        ggplot2::scale_size_continuous(name = "Number of point pairs", range = c(1, 10)) +
        ggplot2::geom_line(data = fit_line_resid, ggplot2::aes(x = dist, y = gamma), color = "red", linewidth = 1) +
        ggplot2::labs(x = "Average bin distance (degrees)", y = "Semi-variance",
                      title = paste0("Fitted residual variogram (best: ", fit_resid$best_meta$model, ")"))
      plots$variogram_residuals_fit <- p_fit_resid
      if (save_plots) ggplot2::ggsave(paste0(plot_prefix, "fitted_variogram_residuals.png"), p_fit_resid, width = 10, height = 10)
    }

    world_plt <- if (is.null(world)) ggplot2::map_data("world") else world
    xlim_plt <- if (is.null(xlim)) range(prediction_grid[[coords[1]]], na.rm = TRUE) else xlim
    ylim_plt <- if (is.null(ylim)) range(prediction_grid[[coords[2]]], na.rm = TRUE) else ylim
    df_map <- data.frame(
      prediction_grid[[coords[1]]],
      prediction_grid[[coords[2]]],
      grid_pred,
      stringsAsFactors = FALSE
    )
    names(df_map) <- c(coords[1], coords[2], value_var)
    p_spatial <- ggplot2::ggplot(df_map, ggplot2::aes(x = .data[[coords[1]]], y = .data[[coords[2]]], fill = .data[[value_var]])) +
      ggplot2::geom_raster() +
      ggplot2::geom_polygon(data = world_plt, ggplot2::aes(x = long, y = lat, group = group), fill = "grey40", color = NA, alpha = 0.1) +
      ggplot2::geom_point(data = dat_resid, ggplot2::aes(x = .data[[coords[1]]], y = .data[[coords[2]]], fill = .data[[value_var]]), shape = 21, color = "black", stroke = 0.5, size = 2) +
      ggplot2::scale_fill_viridis_c(option = "viridis", name = "Carbon density") +
      ggplot2::coord_cartesian(xlim = xlim_plt, ylim = ylim_plt) +
      ggplot2::labs(x = "Longitude", y = "Latitude", title = "Kriged predictions (regression kriging)")
    plots$spatial <- p_spatial
    if (save_plots) ggplot2::ggsave(paste0(plot_prefix, "kriged_predictions_reg_spatial.png"), p_spatial, width = 10, height = 10)

    # Uncertainty map (prediction SE): same structure as predictions (raster + land + points)
    df_unc_reg <- data.frame(
      prediction_grid[[coords[1]]],
      prediction_grid[[coords[2]]],
      prediction_se = sqrt(grid_var),
      stringsAsFactors = FALSE
    )
    names(df_unc_reg)[1:2] <- coords
    p_uncertainty_reg <- ggplot2::ggplot(df_unc_reg, ggplot2::aes(x = .data[[coords[1]]], y = .data[[coords[2]]], fill = prediction_se)) +
      ggplot2::geom_raster() +
      ggplot2::geom_polygon(data = world_plt, ggplot2::aes(x = long, y = lat, group = group), fill = "grey40", color = NA, alpha = 0.1) +
      ggplot2::geom_point(data = dat_resid, ggplot2::aes(x = .data[[coords[1]]], y = .data[[coords[2]]], fill = predicted_error), shape = 21, color = "black", stroke = 0.5, size = 2) +
      ggplot2::scale_fill_viridis_c(option = "plasma", name = "Prediction SE") +
      ggplot2::coord_cartesian(xlim = xlim_plt, ylim = ylim_plt) +
      ggplot2::labs(x = "Longitude", y = "Latitude", title = "Kriged prediction uncertainty (SE, regression kriging)")
    plots$spatial_uncertainty <- p_uncertainty_reg
    if (save_plots) ggplot2::ggsave(paste0(plot_prefix, "kriged_uncertainty_reg_spatial.png"), p_uncertainty_reg, width = 10, height = 10)

    obs <- dat_resid[[value_var]]
    pred <- dat_resid$predicted
    use <- stats::complete.cases(obs) & stats::complete.cases(pred)
    n_use <- sum(use)
    if (n_use > 0) {
      r2 <- stats::cor(obs[use], pred[use])^2
      adj_r2 <- 1 - (1 - r2) * (n_use - 1) / (n_use - 2)
      rmse <- sqrt(mean((obs[use] - pred[use])^2))
      mae <- mean(abs(obs[use] - pred[use]))
      perf_label <- sprintf("R² = %.3f  adj R² = %.3f  RMSE = %.4f  MAE = %.4f  n = %d", r2, adj_r2, rmse, mae, n_use)
    } else {
      perf_label <- "No valid pairs for metrics"
    }
    p_scatter <- ggplot2::ggplot(dat_resid, ggplot2::aes(x = .data[[value_var]], y = predicted)) +
      ggplot2::geom_point(alpha = 0.7, color = "darkblue") +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = predicted - predicted_error, ymax = predicted + predicted_error), width = 0.001, alpha = 0.1) +
      ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
      ggplot2::annotate("label", x = -Inf, y = Inf, label = perf_label, hjust = -0.05, vjust = 1.1, size = 3.5, fill = "white", alpha = 0.9) +
      ggplot2::labs(x = "Observed", y = "Predicted (regression kriging)", title = "Observed vs Predicted (at data locations)") +
      ggplot2::theme_minimal()
    plots$scatter <- p_scatter
    if (save_plots) ggplot2::ggsave(paste0(plot_prefix, "kriged_predictions_reg_scatter.png"), p_scatter, width = 8, height = 8)
  }

  out <- list(
    lm_fit = lm_fit,
    dat_resid = dat_resid,
    vario_resid = vario_resid,
    fit_resid = fit_resid,
    resid_kriged = resid_result$kriged,
    grid_pred = grid_pred,
    grid_var = grid_var,
    prediction_grid = prediction_grid,
    dat_locations_with_pred = dat_resid,
    plots = plots
  )
  return(out)
}
