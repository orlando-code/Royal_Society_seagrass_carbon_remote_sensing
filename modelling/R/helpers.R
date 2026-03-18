# ================================ UTILITY FUNCTIONS ================================

#' Install packages if not already installed
#'
#' @param packages character vector of package names to install
#' @return NULL but prints message to console when packages are installed
install_packages <- function(packages) {
  # if string, convert to list
  if (is.character(packages)) {
    packages <- list(packages)
  }
  for (package in packages) {
    if (!requireNamespace(package, quietly = TRUE)) {
      install.packages(package)
    }
  }
  print(paste("Installed packages:", paste(packages, collapse = ", ")))
}

#' Load packages, installing missing ones; skips already-loaded packages.
#' @param packages character vector of required package names
#' @param optional character vector of package names to warn (not stop) if unavailable
load_packages <- function(packages, optional = NULL) {
  loaded <- character()
  for (pkg in packages) {
    if (isNamespaceLoaded(pkg)) next
    is_opt <- pkg %in% optional
    if (!requireNamespace(pkg, quietly = TRUE)) {
      tryCatch(
        install.packages(pkg, quiet = TRUE),
        error = function(e) if (!is_opt) stop("Failed to install '", pkg, "': ", conditionMessage(e))
      )
    }
    if (suppressWarnings(require(pkg, character.only = TRUE, quietly = TRUE))) {
      loaded <- c(loaded, pkg)
    } else if (is_opt) {
      warning("Optional package '", pkg, "' could not be loaded.")
    } else {
      stop("Required package '", pkg, "' could not be loaded.")
    }
  }
  if (length(loaded) > 0) cat("Loaded package(s):", paste(loaded, collapse = ", "), "\n")
}


safe_write_csv <- function(x, path) {
  tryCatch(
    { write.csv(x, path, row.names = FALSE); cat("Saved", path, "\n") },
    error = function(e) cat("Warning: failed to write", path, ":", e$message, "\n")
  )
}
### CACHING
#' Get or create spatial CV folds with robust caching based on coordinates and parameters.
#'
#' @param dat Data frame with \code{longitude} and \code{latitude} columns
#' @param block_size Block size (metres) for \code{blockCV::cv_spatial}
#' @param n_folds Number of folds
#' @param cache_tag Character tag to distinguish different uses (e.g. "permutation_importance")
#' @param exclude_regions Character vector of Excluded region(s) (for key transparency only)
#' @param progress Logical; passed to \code{blockCV::cv_spatial(progress = ...)}
#' @return List with \code{fold_indices}, \code{cache_path}, and \code{used_cache}
get_cached_spatial_folds <- function(dat,
                                     block_size,
                                     n_folds,
                                     cache_tag,
                                     exclude_regions = character(0),
                                     progress = TRUE) {
  if (!all(c("longitude", "latitude") %in% names(dat))) {
    stop("get_cached_spatial_folds: dat must contain 'longitude' and 'latitude' columns.")
  }

  coord_df <- dat[, c("longitude", "latitude")]
  dat_hash <- digest::digest(coord_df, algo = "sha256")
  cache_key <- list(
    dat_hash = dat_hash,
    exclude_regions = sort(exclude_regions),
    block_size = as.integer(block_size),
    n_folds = as.integer(n_folds)
  )
  key_hash <- digest::digest(cache_key, algo = "sha256")
  cache_dir <- if (requireNamespace("here", quietly = TRUE)) {
    file.path(here::here(), "output/cache")
  } else {
    "output/cache"
  }
  cache_path <- file.path(cache_dir, paste0(cache_tag, "_", key_hash, "_folds.rds"))

  fold_indices <- NULL
  used_cache <- FALSE

  if (file.exists(cache_path)) {
    cached <- tryCatch(readRDS(cache_path), error = function(e) NULL)
    if (is.list(cached) &&
      !is.null(cached$cache_key) &&
      !is.null(cached$fold_indices) &&
      identical(cached$cache_key, cache_key) &&
      length(cached$fold_indices) == nrow(dat)) {
      fold_indices <- cached$fold_indices
      used_cache <- TRUE
      cat("  Using cached spatial folds: ", cache_path, "\n", sep = "")
    } else {
      cat("  Cache found but does not match current data or parameters. Regenerating folds...\n")
    }
  }

  if (is.null(fold_indices)) {
    core_sf <- sf::st_as_sf(dat, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)
    cv_spatial <- blockCV::cv_spatial(
      x = core_sf, k = n_folds, size = block_size,
      selection = "random",
      iteration = 50,
      progress = progress,
      biomod2 = FALSE,
      hexagon = FALSE,
      plot = FALSE
    )
    fold_indices <- cv_spatial$folds_ids
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    saveRDS(list(cache_key = cache_key, fold_indices = fold_indices), cache_path)
    cat("  Cached spatial folds to: ", cache_path, "\n", sep = "")
  }

  list(fold_indices = fold_indices, cache_path = cache_path, used_cache = used_cache)
}

#' Create a generic progress bar e.g. for hyperparameter tuning
#'
#' @param total total number of iterations for the progress bar
#' @param format format string for the progress bar display
#' @param clear whether to clear the progress bar when done
#' @param width width of the progress bar in characters
#' @return a progress bar object from the progress package
progress_bar <- function(total, format = "[:bar] :current/:total (:percent)", clear = FALSE, width = 60) {
  if (!requireNamespace("progress", quietly = TRUE)) {
    install.packages("progress")
  }
  progress::progress_bar$new(total = total, format = format, clear = clear, width = width)
}


# ================================ DATA PREPARATION ================================

#' Process remote sensing covariates by clipping negative values to zero.
#' Call this after loading extracted raster data and before any training or prediction.
#'
#' @param dat Data frame with remote sensing (and possibly other) covariates.
#' @param rs_covariates Character vector of column names to clip (default: common OLCI product names).
#' @return The same data frame with the specified columns clipped to \code{pmax(x, 0)}; columns not present are ignored.
process_rs_covariates <- function(dat, rs_covariates = c(
                                    "S3A_OLCI_ERRNT.20160425_20241231.L3m.CU.CHL.chlor_a.4km",
                                    "S3A_OLCI_ERRNT.20160425_20241231.L3m.CU.IOP.adg_443.4km",
                                    "S3A_OLCI_ERRNT.20160425_20241231.L3m.CU.IOP.bbp_443.4km",
                                    "S3A_OLCI_ERRNT.20160425_20241231.L3m.CU.KD.Kd_490.4km",
                                    "S3A_OLCI_ERRNT.20160425_20241231.L3m.CU.RRS.Rrs_443.4km",
                                    "S3A_OLCI_ERRNT.20160425_20241231.L3m.CU.RRS.Rrs_620.4km"
                                  )) {
  existing_vars <- intersect(rs_covariates, colnames(dat))
  if (length(existing_vars) == 0) {
    return(dat)
  }
  dat <- as.data.frame(dat)
  for (v in existing_vars) {
    if (is.numeric(dat[[v]])) {
      dat[[v]] <- pmax(dat[[v]], 0, na.rm = FALSE)
    }
  }
  dat
}

# ================================ FEATURE SELECTION ================================

#' Check if a variable is a seagrass species variable
#'
#' @param v variable name
#' @return TRUE if the variable is a seagrass species variable, FALSE otherwise
is_species_var <- function(v) v == "seagrass_species" | startsWith(as.character(v), "seagrass_species")

#' Prune predictors by pairwise correlation threshold
#'
#' Removes one variable from each pair with |correlation| > threshold, keeping the one with higher absolute correlation to the target.
#' @param data data.frame with predictors and target
#' @param predictor_vars character vector of predictor names
#' @param target_var name of target variable
#' @param cor_threshold absolute correlation threshold (default 0.8)
#' @return character vector of pruned predictor names
prune_by_correlation <- function(data, predictor_vars, target_var, cor_threshold = 0.8) {
  # Defensive: ensure target_var is a single character string
  if (is.list(target_var) || length(target_var) != 1) {
    target_var <- as.character(target_var)[1]
  }
  dat <- data[, c(target_var, predictor_vars), drop = FALSE]
  dat <- dat[complete.cases(dat), ]
  if (nrow(dat) < 3 || length(predictor_vars) < 2) {
    return(predictor_vars)
  }
  cor_with_target <- sapply(predictor_vars, function(v) {
    if (v %in% names(dat)) {
      stats::cor(dat[[v]], dat[[target_var]], use = "complete.obs")
    } else {
      NA_real_
    }
  })
  cor_df <- data.frame(
    variable = predictor_vars,
    correlation = cor_with_target,
    abs_correlation = abs(cor_with_target),
    stringsAsFactors = FALSE
  )
  cor_df <- cor_df[!is.na(cor_df$correlation), ]
  vars_for_cor <- cor_df$variable
  if (length(vars_for_cor) < 2) {
    return(vars_for_cor)
  }
  covar_data <- dat[, vars_for_cor, drop = FALSE]
  cor_matrix <- stats::cor(covar_data, use = "complete.obs")
  high_cor_pairs <- which(abs(cor_matrix) > cor_threshold & abs(cor_matrix) < 1, arr.ind = TRUE)
  vars_removed <- character()
  if (nrow(high_cor_pairs) > 0) {
    pair_df <- data.frame(
      var1 = rownames(cor_matrix)[high_cor_pairs[, 1]],
      var2 = colnames(cor_matrix)[high_cor_pairs[, 2]],
      correlation = cor_matrix[high_cor_pairs],
      stringsAsFactors = FALSE
    )
    pair_df <- pair_df[pair_df$var1 < pair_df$var2, ]
    pair_df <- pair_df[order(abs(pair_df$correlation), decreasing = TRUE), ]
    for (i in seq_len(nrow(pair_df))) {
      v1 <- pair_df$var1[i]
      v2 <- pair_df$var2[i]
      if (v1 %in% vars_removed || v2 %in% vars_removed) next
      c1 <- cor_df$abs_correlation[cor_df$variable == v1]
      c2 <- cor_df$abs_correlation[cor_df$variable == v2]
      if (length(c1) == 0 || length(c2) == 0) next
      if (round(c1, 3) == round(c2, 3)) next
      if (c1 < c2) {
        vars_removed <- c(vars_removed, v1)
      } else {
        vars_removed <- c(vars_removed, v2)
      }
    }
  }
  setdiff(vars_for_cor, unique(vars_removed))
}

permutation_importance_cv <- function(core_data, target_var, predictor_vars, model_name,
                                      fold_indices, n_permutations = 1L,
                                      verbose = TRUE, log_response = TRUE, hyperparams = NULL) {
  if (is.list(fold_indices)) {
    fvec <- integer(nrow(core_data))
    for (k in seq_along(fold_indices)) fvec[fold_indices[[k]]] <- k
    fold_indices <- fvec
  }
  n_folds <- max(fold_indices)
  if (n_folds < 2L) stop("Need at least 2 folds for CV")
  if (!target_var %in% names(core_data))
    stop("target_var \"", target_var, "\" not found in core_data")
  predictor_vars <- intersect(predictor_vars, names(core_data))
  if (length(predictor_vars) < 2L)
    stop("permutation_importance_cv: need at least 2 predictor columns present in core_data")

  keep_cols <- c(predictor_vars, target_var,
                 intersect(c("longitude", "latitude"), names(core_data)))
  last_error <- NULL
  hp <- hyperparams

  rmse_for_split <- function(train, test) {
    train <- train[complete.cases(train), , drop = FALSE]
    test  <- test[complete.cases(test), , drop = FALSE]
    if (nrow(train) < 2L || nrow(test) < 1L) return(NA_real_)
    test_orig_y <- test[[target_var]]
    train$median_carbon_density <- train[[target_var]]
    test$median_carbon_density  <- test[[target_var]]
    if (log_response) {
      train <- transform_response(train, "median_carbon_density", log = TRUE)
      test  <- transform_response(test, "median_carbon_density", log = TRUE)
    }
    err <- tryCatch({
      prep <- prepare_data_for_model(model_name, train, test, predictor_vars)
      fit <- switch(model_name,
        GPR = fit_gpr(prep$train, prep$predictor_vars, test_data = prep$test, hyperparams = hp),
        RF  = fit_rf(prep$train, prep$test, prep$predictor_vars),
        XGB = fit_xgboost(prep$train, prep$test, prep$predictor_vars, hyperparams = hp),
        GAM = fit_gam(prep$train, prep$test, prep$predictor_vars, include_spatial = FALSE,
                      k_covariate = if (is.list(hp) && !is.null(hp$k_covariate)) hp$k_covariate else 6L),
        stop("Unsupported model: ", model_name)
      )
      pred <- fit$predictions
      if (log_response) pred <- inverse_response_transform(pred, log = TRUE)
      sqrt(mean((test_orig_y - pred)^2, na.rm = TRUE))
    }, error = function(e) {
      last_error <<- conditionMessage(e)
      NA_real_
    })
    err
  }

  if (verbose) cat("  Baseline ", model_name, " RMSE (", n_folds, " folds)... ", sep = "")
  baseline <- vapply(seq_len(n_folds), function(k)
    rmse_for_split(core_data[fold_indices != k, keep_cols, drop = FALSE],
                   core_data[fold_indices == k, keep_cols, drop = FALSE]),
    numeric(1))
  baseline_rmse <- mean(baseline, na.rm = TRUE)
  if (!is.finite(baseline_rmse)) {
    msg <- paste0("Baseline RMSE could not be computed for ", model_name, ".")
    if (!is.null(last_error)) msg <- paste0(msg, " Last error: ", last_error)
    stop(msg)
  }
  if (verbose) cat("done (", round(baseline_rmse, 4), ").\n", sep = "")

  n_vars <- length(predictor_vars)
  if (verbose) { pb <- utils::txtProgressBar(0, n_vars, style = 3, width = 50); on.exit(close(pb), add = TRUE) }

  imp_list <- vector("list", n_vars)
  for (i in seq_len(n_vars)) {
    v <- predictor_vars[i]
    core_perm <- as.data.frame(core_data)
    rmse_inc <- vapply(seq_len(n_permutations), function(perm) {
      core_perm[[v]] <- sample(core_perm[[v]])
      mean(vapply(seq_len(n_folds), function(k)
        rmse_for_split(core_perm[fold_indices != k, keep_cols, drop = FALSE],
                       core_perm[fold_indices == k, keep_cols, drop = FALSE]),
        numeric(1)), na.rm = TRUE) - baseline_rmse
    }, numeric(1))
    imp_list[[i]] <- data.frame(variable = v, rmse_increase = mean(rmse_inc), stringsAsFactors = FALSE)
    if (verbose) utils::setTxtProgressBar(pb, i)
  }

  imp_df <- do.call(rbind, imp_list)
  imp_df[order(-imp_df$rmse_increase), , drop = FALSE]
}

#' Select the top variables from an importance data frame.
#'
#' Keeps variables that together account for \code{coverage} fraction of total
#' importance, capped at \code{max_vars}. Always keeps at least 2.
#'
#' @param df        data.frame with columns \code{variable} and \code{value_col}.
#' @param value_col Name of the numeric importance column.
#' @param max_vars  Maximum number of variables to retain.
#' @param coverage  Cumulative importance coverage threshold (0–1).
#' @return Character vector of selected variable names.
select_top_vars <- function(df, value_col, max_vars = 15L, coverage = 0.99) {
  df  <- df[order(-df[[value_col]]), , drop = FALSE]
  cum <- cumsum(pmax(df[[value_col]], 0))
  tot <- max(cum, na.rm = TRUE)
  cum <- if (tot > 0) cum / tot else seq_along(cum) / length(cum)
  n_keep <- max(2L, min(as.integer(max_vars),
                        sum(cum <= coverage, na.rm = TRUE) + 1L,
                        nrow(df)))
  df$variable[seq_len(n_keep)]
}

#' Print a tidy ranked importance table to the console.
#'
#' @param df         data.frame with columns \code{variable} and \code{value_col}.
#' @param value_col  Name of the numeric importance column.
#' @param method_label Short label for the method (e.g. "Permutation", "SHAP").
#' @param model_name   Model name string, used in the header line.
#' @param max_show     Maximum rows to print (default 15).
print_importance <- function(df, value_col, method_label, model_name,
                             max_show = 15L) {
  df  <- df[order(-df[[value_col]]), , drop = FALSE]
  top <- head(df, as.integer(max_show))
  cat(sprintf("    %s – %s (top %d):\n", method_label, model_name, nrow(top)))
  for (i in seq_len(nrow(top)))
    cat(sprintf("      %2d. %-45s  %+.4f\n",
                i, top$variable[i], top[[value_col]][i]))
}

#' Compute model-agnostic SHAP feature importance via iml::Shapley.
#'
#' Fits a quick model of type \code{model_name} on \code{core_data} and
#' estimates global importance by averaging |phi| over \code{n_points} sampled
#' observations. For GPR, training rows are capped at \code{max_gpr_train} for
#' speed (GPR is O(n^3) to fit).
#'
#' @param core_data  Data frame with target_var, predictor_vars (+ lon/lat for GPR).
#' @param target_var Response variable name.
#' @param predictor_vars Character vector of predictor names.
#' @param model_name One of \code{"XGB"}, \code{"GAM"}, \code{"RF"}, \code{"GPR"}.
#' @param n_points   Observations to sample for SHAP (default 30; higher = slower).
#' @param max_gpr_train Max rows used to fit GPR (subsampled; default 300).
#' @param log_response If TRUE (default), fit and SHAP are on log(response); use same as rest of pipeline.
#' @param hyperparams Optional list of best hyperparameters (e.g. from tuning). XGB: nrounds, max_depth, learning_rate, subsample, colsample_bytree; GAM: k_spatial; GPR: kernel, nug.min, nug.max.
#' @return data.frame(variable, shap_importance) sorted descending, or NULL on failure.
compute_shap_importance <- function(core_data, target_var, predictor_vars,
                                    model_name, n_points = 30L,
                                    max_gpr_train = 300L,
                                    log_response = TRUE,
                                    hyperparams = NULL) {
  load_packages("iml")
  hp <- hyperparams

  keep_cols <- intersect(
    c(target_var, predictor_vars, "longitude", "latitude"), names(core_data)
  )
  core <- core_data[complete.cases(core_data[, keep_cols, drop = FALSE]),
                    keep_cols, drop = FALSE]
  if (nrow(core) < 20L) {
    warning("Too few complete cases for SHAP (", nrow(core), "); returning NULL.")
    return(NULL)
  }

  core$median_carbon_density <- core[[target_var]]
  if (log_response) core <- transform_response(core, "median_carbon_density", log = TRUE)
  y <- core$median_carbon_density

  if (model_name == "XGB") {
    core <- as.data.frame(core)
    prep <- prepare_predictors_train(core, predictor_vars)
    core_sc <- prep$data
    core_sc$median_carbon_density <- y
    pvars <- predictor_vars
    X <- core_sc[, pvars, drop = FALSE]
    fit <- fit_xgboost(core_sc, core_sc, pvars, hyperparams = hp)
    mdl <- fit$model
    pred_fun <- function(m, newdata) as.numeric(predict(m, newdata = as.matrix(newdata)))

  } else if (model_name == "GAM" || model_name == "RF") {
    prep <- prepare_predictors_train_numeric_only(core, predictor_vars)
    core_sc <- prep$data
    core_sc$median_carbon_density <- y
    pvars <- predictor_vars
    if (model_name == "GAM") {
      k_cov <- if (is.list(hp) && !is.null(hp$k_covariate)) hp$k_covariate else 6L
      fit <- fit_gam(core_sc, core_sc, pvars, include_spatial = FALSE, k_covariate = k_cov)
      mdl <- fit$model
      if (is.null(mdl)) { warning("GAM SHAP failed; returning NULL."); return(NULL) }
      pred_fun <- function(m, newdata) {
        as.numeric(mgcv::predict.gam(m, newdata = as.data.frame(newdata), type = "response"))
      }
      X <- core_sc[, pvars, drop = FALSE]
    } else {
      X <- core_sc[, pvars, drop = FALSE]
      fit <- fit_rf(core_sc, core_sc, predictor_vars)
      mdl <- fit$model
      pred_fun <- function(m, newdata) {
        nd <- apply_scaling(as.data.frame(newdata), prep$scale_params, predictor_vars)
        as.numeric(predict(m, newdata = nd))
      }
    }

  } else if (model_name == "GPR") {
    pvars <- predictor_vars
    train_core <- if (nrow(core) > max_gpr_train) {
      core[sample(nrow(core), max_gpr_train), , drop = FALSE]
    } else core
    train_X <- train_core[, predictor_vars, drop = FALSE]
    prep    <- prepare_predictors_train(train_X, predictor_vars)
    train_sc <- prep$data
    train_df <- cbind(median_carbon_density = train_core$median_carbon_density,
                      train_sc)
    X <- train_sc[, pvars, drop = FALSE]
    y_gpr <- train_core$median_carbon_density  # y must match nrow(X) for iml::Predictor
    form_str <- paste("median_carbon_density ~",
                      paste(predictor_vars, collapse = " + "))
    kernel  <- if (is.list(hp) && !is.null(hp$kernel)) hp$kernel else "matern52"
    nug_min <- if (is.list(hp) && !is.null(hp$nug.min)) hp$nug.min else 1e-8
    nug_max <- if (is.list(hp) && !is.null(hp$nug.max)) hp$nug.max else 100
    mdl <- tryCatch(
      GauPro::gpkm(as.formula(form_str), data = train_df,
                   kernel = kernel, nug.min = nug_min, nug.max = nug_max),
      error = function(e) NULL
    )
    if (is.null(mdl)) {
      warning("GPR fitting failed in compute_shap_importance; returning NULL.")
      return(NULL)
    }
    pred_fun  <- function(m, newdata) {
      Xn <- as.matrix(newdata[, predictor_vars, drop = FALSE])
      storage.mode(Xn) <- "double"
      as.numeric(m$pred(Xn, se.fit = FALSE))
    }
    y <- y_gpr
  } else {
    stop("compute_shap_importance: unsupported model_name '", model_name, "'")
  }

  predictor_iml <- iml::Predictor$new(
    model            = mdl,
    data             = X,
    y                = y,
    predict.function = pred_fun
  )

  n_pts <- min(as.integer(n_points), nrow(X))
  idx   <- sample(seq_len(nrow(X)), n_pts)
  shap_mat <- matrix(0, nrow = n_pts, ncol = length(pvars), dimnames = list(NULL, pvars))
  n_ok <- 0L
  for (i in seq_len(n_pts)) {
    sh <- tryCatch(iml::Shapley$new(predictor_iml, x.interest = X[idx[i], , drop = FALSE]),
                   error = function(e) NULL)
    if (is.null(sh)) next
    phi_vals <- as.numeric(sh$results$phi)
    # iml may return feature names that don't match pvars (e.g. make.names: _ to .); match by position
    if (length(phi_vals) == length(pvars)) {
      shap_mat[i, ] <- abs(phi_vals)
      n_ok <- n_ok + 1L
    } else {
      feat_names <- as.character(sh$results$feature)
      common <- intersect(feat_names, pvars)
      if (length(common) == 0L) {
        # try matching after normalizing (dots/underscores)
        norm_pvars <- gsub("[^a-zA-Z0-9]", ".", pvars)
        norm_feat <- gsub("[^a-zA-Z0-9]", ".", feat_names)
        for (j in seq_along(pvars)) {
          match_j <- match(norm_pvars[j], norm_feat)
          if (!is.na(match_j)) shap_mat[i, pvars[j]] <- abs(phi_vals[match_j])
        }
        n_ok <- n_ok + 1L
      } else {
        shap_mat[i, common] <- abs(phi_vals[match(common, feat_names)])
        n_ok <- n_ok + 1L
      }
    }
  }
  if (n_ok == 0L) {
    warning("SHAP: no Shapley computations succeeded for ", model_name, "; returning NULL.")
    return(NULL)
  }
  imp <- colMeans(shap_mat, na.rm = TRUE)
  data.frame(variable = names(imp), shap_importance = as.numeric(imp), stringsAsFactors = FALSE)[order(-imp), , drop = FALSE]
}

#' Combine per-model variable/importance CSVs into a single file.
#'
#' \code{type = "pruned"}: \code{pruned_variables_to_include_<model>.csv} ->
#' \code{pruned_model_variables.csv} (model, variable).
#' \code{type = "permutation"} or \code{"shap"}: \code{importance_<type>_<model>.csv} ->
#' \code{importance_<type>_combined.csv} (model, variable, value column).
#'
#' @param covariate_dir Directory containing the per-model CSVs.
#' @param type One of \code{"pruned"}, \code{"permutation"}, \code{"shap"}.
#' @return Path to the combined CSV, or \code{NULL} if no per-model files found.
combine_pruned_model_variables <- function(covariate_dir,
                                           type = c("pruned", "permutation", "shap")) {
  type <- match.arg(type)
  if (type == "pruned") {
    pat <- "^pruned_variables_to_include_.+\\.csv$"
    out <- file.path(covariate_dir, "pruned_model_variables.csv")
    value_col <- NULL
  } else {
    pat <- paste0("^importance_", type, "_.+\\.csv$")
    out <- file.path(covariate_dir, paste0("importance_", type, "_combined.csv"))
    value_col <- if (type == "permutation") "rmse_increase" else "shap_importance"
  }
  files <- setdiff(list.files(covariate_dir, pattern = pat, full.names = TRUE), out)
  if (length(files) == 0L) return(invisible(NULL))
  combined_list <- lapply(files, function(f) {
    m <- sub(paste0("^(pruned_variables_to_include_|importance_", type, "_)(.+)\\.csv$"),
             "\\2", basename(f))
    d <- read.csv(f, stringsAsFactors = FALSE)
    if (!"variable" %in% names(d)) return(NULL)
    if (is.null(value_col)) return(data.frame(model = m, variable = d$variable, stringsAsFactors = FALSE))
    if (!value_col %in% names(d)) return(NULL)
    out_df <- data.frame(model = m, variable = d$variable, d[[value_col]], stringsAsFactors = FALSE)
    names(out_df)[3L] <- value_col
    out_df
  })
  combined_list <- combined_list[!vapply(combined_list, is.null, logical(1L))]
  if (length(combined_list) == 0L) return(invisible(NULL))
  combined <- do.call(rbind, combined_list)
  write.csv(combined, out, row.names = FALSE)
  invisible(out)
}

#' Compare covariates between models
#'
#' Given a data frame with a 'model' and 'variable' column (for n models), prints and returns
#' a list of variables grouped by how many models share them (e.g. shared by 2, 3 models, etc.),
#' and, for each variable, indicates which models share it.
#'
#' @param model_variables_df data frame with 'model' and 'variable' columns
#' @return named list: each name is number of models that share variables; each entry is a data.frame
#'         with columns 'variable' and 'models' (character vector of model names sharing that variable)
compare_covariates_between_models <- function(model_variables_df) {
  if (!requireNamespace("dplyr", quietly = TRUE))
    stop("Package 'dplyr' is required.")
  out <- model_variables_df %>%
    dplyr::group_by(.data$variable) %>%
    dplyr::summarise(
      n_models = dplyr::n_distinct(.data$model),
      models   = list(sort(unique(as.character(.data$model)))),
      .groups  = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(.data$n_models), .data$variable)

  shared_list <- split(out[, c("variable", "models")], out$n_models)
  for (k in sort(unique(out$n_models), decreasing = TRUE)) {
    vars_tab <- shared_list[[as.character(k)]]
    cat(sprintf("Variables shared by %d model%s:\n", k, if (k == 1) "" else "s"))
    for (i in seq_len(nrow(vars_tab)))
      cat(sprintf("  %s (%s)\n", vars_tab$variable[i],
                  paste(vars_tab$models[[i]], collapse = ", ")))
    cat("\n")
  }
  invisible(shared_list)
}



# ================================ INITIALISATION ================================

#' Index of test rows whose factor predictor levels were present in training
#'
#' Used so CV metrics are computed only on test rows we can fairly evaluate
#' (e.g. species seen in the training fold). Does not change train/test split.
#'
#' @param train training data frame (with factor columns)
#' @param test test data frame (same factor column names)
#' @param factor_vars character vector of factor predictor names (e.g. \code{c("seagrass_species", "region")})
#' @return logical vector of length \code{nrow(test)}: TRUE where all factor_vars (present in both frames) have values in train
test_rows_with_factors_in_train <- function(train, test, factor_vars = c("seagrass_species", "region")) {
  keep <- rep(TRUE, nrow(test))
  for (col in factor_vars) {
    if (col %in% names(train) && col %in% names(test)) {
      train_vals <- unique(train[[col]])
      keep <- keep & (test[[col]] %in% train_vals)
    }
  }
  keep
}

# Original R² metric (squared correlation; kept for reference)
# calculate_metrics <- function(observed, predicted) {
#   r2 <- cor(observed, predicted, use = "complete.obs")^2
#   rmse <- sqrt(mean((observed - predicted)^2, na.rm = TRUE))
#   mae <- mean(abs(observed - predicted), na.rm = TRUE)
#   bias <- mean(predicted - observed, na.rm = TRUE)
#   data.frame(r2 = r2, rmse = rmse, mae = mae, bias = bias)
# }

#' Calculate performance metrics for model predictions
#'
#' R² is computed as the standard coefficient of determination
#'  \(1 - \sum (y - \hat{y})^2 / \sum (y - \bar{y})^2\), which can be negative
#'  when the model performs worse than using the mean.
#'
#' @param observed vector of observed (true) values
#' @param predicted vector of predicted values
#' @return data frame with r2, rmse, mae, and bias metrics
calculate_metrics <- function(observed, predicted) {
  ok <- is.finite(observed) & is.finite(predicted)
  observed <- observed[ok]
  predicted <- predicted[ok]

  if (length(observed) < 2L) {
    return(data.frame(
      r2 = NA_real_, rmse = NA_real_, mae = NA_real_, bias = NA_real_
    ))
  }

  ss_res <- sum((observed - predicted)^2)
  ss_tot <- sum((observed - mean(observed))^2)
  r2 <- if (ss_tot > 0) 1 - ss_res / ss_tot else NA_real_

  rmse <- sqrt(mean((observed - predicted)^2))
  mae <- mean(abs(observed - predicted))
  bias <- mean(predicted - observed)

  data.frame(r2 = r2, rmse = rmse, mae = mae, bias = bias)
}

#' Read per-model predictor variables from a CSV (e.g. pruned_model_variables_shap.csv).
#'
#' @param fn Path to CSV with columns \code{model} and \code{variable}.
#' @param valid_cols Optional character vector (e.g. \code{colnames(core_data)}); if provided, variables are intersected with this.
#' @return Named list: model name -> character vector of variables (or \code{NULL} if file missing/invalid).
read_model_vars <- function(fn, valid_cols = NULL) {
  if (!file.exists(fn)) return(NULL)
  df <- read.csv(fn, stringsAsFactors = FALSE)
  if (!all(c("model", "variable") %in% names(df))) return(NULL)
  out <- lapply(split(df$variable, df$model), function(v) {
    if (length(valid_cols)) intersect(v, valid_cols) else v
  })
  out
}

#' Build per-model variable lists from SHAP and permutation CSV files.
#'
#' @param cov_dir Directory containing \code{pruned_model_variables_shap.csv} and \code{pruned_model_variables_perm.csv}.
#' @param valid_cols Optional character vector to restrict variables (e.g. \code{colnames(core_data)}).
#' @param use_shap_first If \code{TRUE}, prefer SHAP file when both exist; else prefer permutation.
#' @return List with elements \code{shap} and \code{perm}, each a named list (model -> variables).
get_per_model_vars <- function(cov_dir, valid_cols = NULL, use_shap_first = TRUE) {
  shap_file <- file.path(cov_dir, "pruned_model_variables_shap.csv")
  perm_file <- file.path(cov_dir, "pruned_model_variables_perm.csv")
  list(
    shap = read_model_vars(shap_file, valid_cols),
    perm = read_model_vars(perm_file, valid_cols)
  )
}

#' Return the predictor variable vector for a model from pre-built SHAP/perm lists.
#'
#' @param model_name Model name (e.g. \code{"GPR"}, \code{"GAM"}, \code{"XGB"}).
#' @param per_model_vars List from \code{get_per_model_vars()} with elements \code{shap} and \code{perm}.
#' @param use_shap_first If \code{TRUE}, use SHAP vars when available and sufficient; else try permutation first.
#' @return Character vector of variable names (at least 2).
load_model_vars <- function(model_name, per_model_vars, use_shap_first = TRUE) {
  src_order <- if (use_shap_first) c("shap", "perm") else c("perm", "shap")
  for (s in src_order) {
    src_list <- per_model_vars[[s]]
    if (!is.list(src_list)) next
    v <- src_list[[model_name]]
    if (!is.null(v) && length(v) >= 2L) return(v)
  }
  stop("No pruned variables for ", model_name, " in SHAP or permutation CSV. Run covariate selection first.")
}

# ================================ HYPERPARAMETER TUNING ================================

#' Unified function to tune hyperparameters for different model types
#'
#' @param model_type character string specifying model type: "rf", "xgb", or "nn"
#' @param train_data training data frame containing predictor variables and response
#' @param predictor_vars character vector of predictor variable names
#' @param n_folds number of folds for cross-validation during tuning (default: 3)
#' @param verbose whether to print detailed tuning progress (default: false)
#' @param fold_indices optional integer vector of fold IDs; if NULL, random folds
#' @param block_size optional block size (m) for spatial CV
#' @param cache_tag tag for spatial fold cache
#' @param exclude_regions passed to get_cached_spatial_folds when block_size is used
#' @return list of best hyperparameters for the specified model type
tune_hyperparameters <- function(model_type, train_data, predictor_vars, n_folds = 3, verbose = FALSE,
                                 fold_indices = NULL, block_size = NULL, cache_tag = NULL,
                                 exclude_regions = character(0)) {
  cat("Tuning", model_type, "hyperparameters...\n")
  tag <- cache_tag %||% paste0("tune_", model_type)
  if (model_type == "rf" || model_type == "random_forest") {
    return(tune_rf(train_data, predictor_vars, n_folds, verbose, fold_indices, block_size, tag, exclude_regions))
  } else if (model_type == "xgb" || model_type == "xgboost") {
    return(tune_xgboost(train_data, predictor_vars, n_folds, verbose, fold_indices, block_size, tag, exclude_regions))
  } else if (model_type == "nn" || model_type == "neural_network") {
    return(tune_nn(train_data, predictor_vars, n_folds, verbose, fold_indices, block_size, tag, exclude_regions))
  } else {
    stop("Unknown model type. Use 'rf', 'xgb', or 'nn'")
  }
}

#' Resolve fold indices: use provided, or spatial (block_size), or random.
#' @return integer vector of fold IDs, length = nrow(dat)
resolve_fold_indices <- function(dat, n_folds, fold_indices = NULL, block_size = NULL,
                                  cache_tag = "tune", exclude_regions = character(0)) {
  if (!is.null(fold_indices)) {
    if (is.list(fold_indices)) {
      fv <- integer(nrow(dat))
      for (k in seq_along(fold_indices)) fv[fold_indices[[k]]] <- k
      return(fv)
    }
    return(as.integer(fold_indices))
  }
  if (!is.null(block_size) && all(c("longitude", "latitude") %in% names(dat))) {
    fi <- get_cached_spatial_folds(dat, block_size, n_folds, cache_tag, exclude_regions, progress = FALSE)
    return(fi$fold_indices)
  }
  sample(rep(seq_len(n_folds), length.out = nrow(dat)))
}

#' Tune random forest hyperparameters using cross-validation
#'
#' @param train_data training data frame containing predictor variables and response
#' @param predictor_vars character vector of predictor variable names
#' @param n_folds number of folds for cross-validation during tuning (default: 3)
#' @param verbose whether to print detailed tuning progress (default: false)
#' @param fold_indices optional integer vector of fold IDs (length = nrow(train_data)); if NULL, random folds
#' @param block_size optional block size (m) for spatial CV; if set with lon/lat in data, uses get_cached_spatial_folds
#' @param cache_tag tag for spatial fold cache (default: "tune")
#' @param exclude_regions passed to get_cached_spatial_folds when block_size is used
#' @return list of best hyperparameters (mtry, ntree, nodesize)
tune_rf <- function(train_data, predictor_vars, n_folds = 3, verbose = FALSE,
                    fold_indices = NULL, block_size = NULL, cache_tag = "tune_rf",
                    exclude_regions = character(0)) {
  formula_str <- paste("median_carbon_density ~", paste(predictor_vars, collapse = " + "))
  formula_obj <- as.formula(formula_str)

  mtry_vals <- c(2, 4, 6, floor(sqrt(length(predictor_vars))), length(predictor_vars))
  ntree_vals <- c(300, 500, 700)
  nodesize_vals <- c(5, 10, 20)

  best_rmse <- Inf
  best_params <- NULL
  folds <- resolve_fold_indices(train_data, n_folds, fold_indices, block_size, cache_tag, exclude_regions)
  n_folds <- max(folds)

  pb <- progress_bar(
    total = length(mtry_vals) * length(ntree_vals) * length(nodesize_vals)
  )
  for (mtry in mtry_vals) {
    for (ntree in ntree_vals) {
      for (nodesize in nodesize_vals) {
        cv_rmse <- c()
        for (fold in seq_len(n_folds)) {
          train_fold <- train_data[folds != fold, ]
          val_fold <- train_data[folds == fold, ]

          model <- randomForest(formula_obj,
            data = train_fold,
            mtry = mtry, ntree = ntree, nodesize = nodesize,
            importance = FALSE
          )
          pred <- predict(model, val_fold)
          cv_rmse <- c(cv_rmse, sqrt(mean((val_fold$median_carbon_density - pred)^2)))
        }

        mean_rmse <- mean(cv_rmse)
        if (verbose) cat(sprintf("  mtry=%d, ntree=%d, nodesize=%d: RMSE=%.4f\n", mtry, ntree, nodesize, mean_rmse))

        if (mean_rmse < best_rmse) {
          best_rmse <- mean_rmse
          best_params <- list(mtry = mtry, ntree = ntree, nodesize = nodesize)
        }
        pb$tick()
      }
    }
  }

  return(best_params)
}

#' Tune xgboost hyperparameters using cross-validation
#'
#' @param train_data training data frame containing predictor variables and response
#' @param predictor_vars character vector of predictor variable names
#' @param n_folds number of folds for cross-validation during tuning (default: 3)
#' @param verbose whether to print detailed tuning progress (default: false)
#' @param fold_indices optional integer vector of fold IDs; if NULL, random folds
#' @param block_size optional block size (m) for spatial CV
#' @param cache_tag tag for spatial fold cache (default: "tune_xgboost")
#' @param exclude_regions passed to get_cached_spatial_folds when block_size is used
#' @return list of best hyperparameters (max_depth, learning_rate, subsample, colsample_bytree, nrounds)
tune_xgboost <- function(train_data, predictor_vars, n_folds = 3, verbose = FALSE,
                         fold_indices = NULL, block_size = NULL, cache_tag = "tune_xgboost",
                         exclude_regions = character(0)) {
  X_train <- as.matrix(train_data[, predictor_vars])
  y_train <- train_data$median_carbon_density

  max_depth_vals <- c(3, 6, 9)
  learning_rate_vals <- c(0.1, 0.3)
  subsample_vals <- c(0.8, 1.0)
  colsample_bytree_vals <- c(0.8, 1.0)
  nrounds_vals <- c(100, 200)

  best_rmse <- Inf
  best_params <- NULL
  folds <- resolve_fold_indices(train_data, n_folds, fold_indices, block_size, cache_tag, exclude_regions)
  n_folds <- max(folds)

  total_combos <- length(max_depth_vals) * length(learning_rate_vals) * length(subsample_vals) *
    length(colsample_bytree_vals) * length(nrounds_vals)
  pb <- progress_bar(
    total = total_combos
  )
  for (max_depth in max_depth_vals) {
    for (learning_rate in learning_rate_vals) {
      for (subsample in subsample_vals) {
        for (colsample_bytree in colsample_bytree_vals) {
          for (nrounds in nrounds_vals) {
            cv_rmse <- c()

            for (fold in seq_len(n_folds)) {
              X_fold_train <- X_train[folds != fold, ]
              y_fold_train <- y_train[folds != fold]
              X_fold_val <- X_train[folds == fold, ]
              y_fold_val <- y_train[folds == fold]

              model <- xgboost(
                x = X_fold_train, y = y_fold_train,
                nrounds = nrounds, max_depth = max_depth,
                learning_rate = learning_rate, subsample = subsample,
                colsample_bytree = colsample_bytree,
                nthread = 1L,
                objective = "reg:squarederror"
              )
              pred <- predict(model, X_fold_val)
              cv_rmse <- c(cv_rmse, sqrt(mean((y_fold_val - pred)^2)))
            }

            mean_rmse <- mean(cv_rmse)
            if (mean_rmse < best_rmse) {
              best_rmse <- mean_rmse
              best_params <- list(
                max_depth = max_depth, learning_rate = learning_rate,
                subsample = subsample, colsample_bytree = colsample_bytree,
                nrounds = nrounds
              )
            }
            pb$tick()
          }
        }
      }
    }
  }

  return(best_params)
}

#' Tune neural network hyperparameters using cross-validation
#'
#' @param train_data training data frame containing predictor variables and response
#' @param predictor_vars character vector of predictor variable names
#' @param n_folds number of folds for cross-validation during tuning (default: 3)
#' @param verbose whether to print detailed tuning progress (default: false)
#' @param fold_indices optional integer vector of fold IDs; if NULL, random folds
#' @param block_size optional block size (m) for spatial CV
#' @param cache_tag tag for spatial fold cache (default: "tune_nn")
#' @param exclude_regions passed to get_cached_spatial_folds when block_size is used
#' @return list of best hyperparameters (hidden, learningrate)
tune_nn <- function(train_data, predictor_vars, n_folds = 3, verbose = FALSE,
                    fold_indices = NULL, block_size = NULL, cache_tag = "tune_nn",
                    exclude_regions = character(0)) {
  train_data <- as.data.frame(train_data)
  train_numeric <- train_data %>%
    dplyr::select(dplyr::all_of(predictor_vars), median_carbon_density) %>%
    dplyr::mutate_all(as.numeric) %>%
    as.data.frame()
  idx <- complete.cases(train_numeric)
  train_numeric <- train_numeric[idx, ]
  folds <- resolve_fold_indices(train_data, n_folds, fold_indices, block_size, cache_tag, exclude_regions)
  folds <- folds[idx]
  n_folds <- max(folds)

  # scale data
  scale_params <- list(
    means = colMeans(train_numeric[, predictor_vars, drop = FALSE]),
    sds = apply(train_numeric[, predictor_vars, drop = FALSE], 2, sd)
  )
  train_scaled <- train_numeric
  for (var in predictor_vars) {
    if (scale_params$sds[var] > 0) {
      train_scaled[[var]] <- (train_numeric[[var]] - scale_params$means[var]) / scale_params$sds[var]
    } else {
      train_scaled[[var]] <- 0
    }
  }
  response_mean <- mean(train_numeric$median_carbon_density)
  response_sd <- sd(train_numeric$median_carbon_density)
  train_scaled$median_carbon_density <- (train_numeric$median_carbon_density - response_mean) / response_sd

  formula_str <- paste("median_carbon_density ~", paste(predictor_vars, collapse = " + "))
  formula_obj <- as.formula(formula_str)

  # parameter grid
  hidden_layers <- list(c(5), c(10), c(5, 5), c(10, 5))
  learningrate_vals <- c(0.01, 0.1)

  best_rmse <- Inf
  best_params <- NULL
  pb <- progress_bar(
    total = length(hidden_layers) * length(learningrate_vals)
  )
  for (hidden in hidden_layers) {
    for (learningrate in learningrate_vals) {
      cv_rmse <- c()

      for (fold in seq_len(n_folds)) {
        train_fold <- train_scaled[folds != fold, ]
        val_fold <- train_scaled[folds == fold, ]

        tryCatch(
          {
            model <- neuralnet(formula_obj,
              data = train_fold,
              hidden = hidden, learningrate = learningrate,
              linear.output = TRUE, threshold = 0.01, stepmax = 1e+05
            )
            pred_scaled <- predict(model, val_fold[, predictor_vars, drop = FALSE])
            if (is.matrix(pred_scaled)) pred_scaled <- pred_scaled[, 1]
            pred <- pred_scaled * response_sd + response_mean
            obs <- val_fold$median_carbon_density * response_sd + response_mean
            cv_rmse <- c(cv_rmse, sqrt(mean((obs - pred)^2)))
          },
          error = function(e) {
            cv_rmse <<- c(cv_rmse, Inf) # skip if model fails
          }
        )
      }

      mean_rmse <- mean(cv_rmse)
      if (verbose) {
        cat(sprintf(
          "  hidden=%s, lr=%.2f: RMSE=%.4f\n",
          paste(hidden, collapse = ","), learningrate, mean_rmse
        ))
      }

      if (mean_rmse < best_rmse && !is.infinite(mean_rmse)) {
        best_rmse <- mean_rmse
        best_params <- list(hidden = hidden, learningrate = learningrate)
      }
      pb$tick()
    }
  }

  return(best_params)
}

#' Tune GPR hyperparameters using cross-validation
#'
#' @param train_data training data frame containing predictor variables and response (median_carbon_density)
#' @param predictor_vars character vector of predictor variable names
#' @param n_folds number of folds for cross-validation during tuning (default: 3)
#' @param verbose whether to print detailed tuning progress (default: false)
#' @param value_var response column name (default: median_carbon_density)
#' @param fold_indices optional integer vector of fold IDs; if NULL, random folds
#' @param block_size optional block size (m) for spatial CV
#' @param cache_tag tag for spatial fold cache (default: "tune_gpr")
#' @param exclude_regions passed to get_cached_spatial_folds when block_size is used
#' @return list of best hyperparameters (kernel, nug.min, nug.max, nug.est) for fit_gpr
tune_gpr <- function(train_data, predictor_vars, n_folds = 3, verbose = FALSE,
                     value_var = "median_carbon_density",
                     fold_indices = NULL, block_size = NULL, cache_tag = "tune_gpr",
                     exclude_regions = character(0)) {
  if (!"median_carbon_density" %in% names(train_data) && value_var %in% names(train_data)) {
    train_data$median_carbon_density <- train_data[[value_var]]
  }
  need_cols <- c(value_var, predictor_vars)
  if (!all(need_cols %in% names(train_data))) {
    stop("tune_gpr: train_data must contain ", paste(need_cols, collapse = ", "))
  }
  idx <- complete.cases(train_data[, need_cols, drop = FALSE])
  dat <- train_data[idx, need_cols, drop = FALSE]
  if (nrow(dat) < 2L * n_folds) {
    return(list(kernel = "matern52", nug.min = 1e-8, nug.max = 100, nug.est = TRUE))
  }
  folds <- resolve_fold_indices(train_data, n_folds, fold_indices, block_size, cache_tag, exclude_regions)
  folds <- folds[idx]
  n_folds <- max(folds)
  kernels <- c("matern52", "matern32", "gaussian", "exponential")
  nug_min_vals <- c(1e-8, 1e-6, 1e-4)
  nug_max_vals <- c(10, 100)
  best_rmse <- Inf
  best_params <- list(kernel = "matern52", nug.min = 1e-8, nug.max = 100, nug.est = TRUE)
  total <- length(kernels) * length(nug_min_vals) * length(nug_max_vals)
  pb <- progress_bar(total = total)
  for (kernel in kernels) {
    for (nug_min in nug_min_vals) {
      for (nug_max in nug_max_vals) {
        if (nug_min >= nug_max) next
        cv_rmse <- numeric(n_folds)
        for (fold in seq_len(n_folds)) {
          train_fold <- dat[folds != fold, , drop = FALSE]
          test_fold <- dat[folds == fold, , drop = FALSE]
          res <- fit_gpr(
            train_data = train_fold,
            test_data = test_fold,
            predictor_vars = predictor_vars,
            value_var = value_var,
            hyperparams = list(kernel = kernel, nug.min = nug_min, nug.max = nug_max)
          )
          cv_rmse[fold] <- if (!is.null(res$rmse) && is.finite(res$rmse)) res$rmse else NA_real_
        }
        mean_rmse <- mean(cv_rmse, na.rm = TRUE)
        if (verbose) cat(sprintf("  kernel=%s, nug.min=%s, nug.max=%g: RMSE=%.4f\n",
                                kernel, format(nug_min, scientific = TRUE), nug_max, mean_rmse))
        if (is.finite(mean_rmse) && mean_rmse < best_rmse) {
          best_rmse <- mean_rmse
          best_params <- list(kernel = kernel, nug.min = nug_min, nug.max = nug_max, nug.est = TRUE)
        }
        pb$tick()
      }
    }
  }
  best_params
}


# ================================ FOLD ASSIGNMENTS FOR CV STRATEGIES ================================

folds_to_vector <- function(folds, n) {
  if (is.list(folds)) {
    fv <- integer(n)
    for (k in seq_along(folds)) fv[folds[[k]]] <- k
    fv
  } else as.integer(folds)
}
#' Build pixel-grouped fold assignments: points with identical covariate
#' vectors are placed in the same fold so the model never sees an exact
#' duplicate input in the test set during training.
#'
#' @param dat       data.frame with at least the columns in \code{covariate_cols}
#' @param covariate_cols character vector of column names to hash
#' @param n_folds   number of folds
#' @param seed      random seed for reproducibility (NULL = no set.seed)
#' @return list(fold_indices, n_groups) where fold_indices is an integer vector
#'   of length nrow(dat) and n_groups is the number of unique covariate vectors
make_pixel_grouped_folds <- function(dat, covariate_cols, n_folds, seed = 42L) {
  cov_mat <- dat[, covariate_cols, drop = FALSE]
  pixel_id <- as.integer(factor(do.call(paste, c(cov_mat, sep = "|"))))
  unique_ids <- unique(pixel_id)
  if (!is.null(seed)) set.seed(seed)
  grp_fold <- sample(rep(seq_len(n_folds), length.out = length(unique_ids)))
  list(
    fold_indices = grp_fold[match(pixel_id, unique_ids)],
    n_groups     = length(unique_ids)
  )
}

#' Create CV fold assignments based on the global cv_type setting.
#'
#' Dispatches to the appropriate fold construction method:
#'   "random"           — naive i.i.d. row-level split
#'   "location_grouped" — group by identical (lon, lat)
#'   "pixel_grouped"    — group by identical covariate vector
#'   "spatial"          — spatial blocks via blockCV (requires cv_blocksize)
#'
#' @param dat            data.frame (must contain longitude, latitude, and covariates)
#' @param covariate_cols character vector of raster covariate column names
#' @param n_folds        integer number of folds
#' @param cv_type        one of "random", "location_grouped", "pixel_grouped", "spatial"
#' @param cv_blocksize   block size in metres (only used when cv_type = "spatial")
#' @param exclude_regions character vector of regions to exclude (spatial only)
#' @param cache_tag      cache tag string for spatial folds (spatial only)
#' @param seed           random seed (NULL = no set.seed)
#' @return list with fold_indices (integer vector, length nrow(dat)),
#'         method_name (character), and n_groups (integer, number of groups)
make_cv_folds <- function(dat, covariate_cols, n_folds, cv_type,
                          cv_blocksize = NULL, exclude_regions = character(0),
                          cache_tag = "cv_folds", seed = 42L) {
  if (!is.null(seed)) set.seed(seed)

  if (identical(cv_type, "random")) {
    fi <- sample(rep(seq_len(n_folds), length.out = nrow(dat)))
    cat("  CV folds: random split (", nrow(dat), " rows -> ", n_folds, " folds).\n", sep = "")
    return(list(fold_indices = fi, method_name = "random_split", n_groups = nrow(dat)))
  }

  if (identical(cv_type, "location_grouped")) {
    loc_id <- as.integer(factor(paste(dat$longitude, dat$latitude)))
    unique_ids <- unique(loc_id)
    grp_fold <- sample(rep(seq_len(n_folds), length.out = length(unique_ids)))
    fi <- grp_fold[match(loc_id, unique_ids)]
    cat("  CV folds: location-grouped random (", length(unique_ids),
        " unique locations -> ", n_folds, " folds).\n", sep = "")
    return(list(fold_indices = fi, method_name = "location_grouped_random", n_groups = length(unique_ids)))
  }

  if (identical(cv_type, "pixel_grouped")) {
    pf <- make_pixel_grouped_folds(dat, covariate_cols, n_folds, seed = NULL)
    cat("  CV folds: pixel-grouped random (", pf$n_groups,
        " unique covariate vectors -> ", n_folds, " folds).\n", sep = "")
    return(list(fold_indices = pf$fold_indices, method_name = "pixel_grouped_random", n_groups = pf$n_groups))
  }

  if (identical(cv_type, "spatial")) {
    if (is.null(cv_blocksize)) stop("cv_type = 'spatial' requires cv_blocksize to be set.")
    fi_info <- get_cached_spatial_folds(
      dat = dat, block_size = cv_blocksize, n_folds = n_folds,
      cache_tag = cache_tag, exclude_regions = exclude_regions, progress = TRUE
    )
    method <- paste0("spatial_block_", cv_blocksize, "m")
    n_actual <- max(fi_info$fold_indices, na.rm = TRUE)
    cat("  CV folds: spatial blocks (", cv_blocksize, " m, ",
        n_actual, " folds).\n", sep = "")
    return(list(fold_indices = fi_info$fold_indices, method_name = method, n_groups = n_actual))
  }

  stop("Unknown cv_type '", cv_type, "'. Use 'random', 'location_grouped', 'pixel_grouped', or 'spatial'.")
}

method_to_display <- function(m) {
  if (grepl("^spatial_block_.+m$", m)) {
    size_str <- sub("^spatial_block_(.+)m$", "\\1", m)
    size_m <- as.numeric(size_str)
    size_km <- if (!is.na(size_m) && size_m >= 1000) round(size_m / 1000) else size_m
    if (!is.na(size_km)) {
      return(paste("Spatial block", size_km, "km"))
    }
  }
  if (grepl("^region_stratified_.+m$", m)) {
    size_str <- sub("^region_stratified_(.+)m$", "\\1", m)
    size_m <- as.numeric(size_str)
    size_km <- if (!is.na(size_m) && size_m >= 1000) round(size_m / 1000) else size_m
    if (!is.na(size_km)) {
      return(paste("Region-strat.", size_km, "km"))
    }
  }
  if (grepl("^distance_buffer_", m)) {
    km <- sub("distance_buffer_|km", "", m)
    return(paste("Distance buffer", km, "km"))
  }
  switch(m,
    random_split = "Random split",
    location_grouped_random = "Location-grouped random",
    pixel_grouped_random = "Pixel-grouped random",
    env_cluster = "Env. cluster",
    paste0(toupper(substring(m, 1, 1)), substring(m, 2))
  )
}
# Order methods: random first (0), then spatial blocks by ascending size (1km, 5km, ...)
method_block_size_km <- function(m) {
  if (m == "random_split") return(-1)
  if (m == "location_grouped_random") return(-0.5)
  if (m == "pixel_grouped_random") return(0)
  if (grepl("^spatial_block_.+m$", m)) {
    size_str <- sub("^spatial_block_(.+)m$", "\\1", m)
    size_m <- as.numeric(size_str)
    if (!is.na(size_m)) return(round(size_m / 1000))
  }
  if (grepl("^region_stratified_.+m$", m)) {
    size_str <- sub("^region_stratified_(.+)m$", "\\1", m)
    size_m <- as.numeric(size_str)
    if (!is.na(size_m)) return(round(size_m / 1000) + 0.5)
  }
  NA_real_
}
order_methods_by_size <- function(methods) {
  method_sizes <- vapply(methods, method_block_size_km, numeric(1))
  order_idx <- order(ifelse(is.na(method_sizes), 999, method_sizes), na.last = TRUE)
  methods[order_idx]
}


# ================================ CROSS-VALIDATION ================================

#' Run cross-validation for multiple models with optional hyperparameter tuning
#'
#' @param cv_method_name character string describing the cv method (e.g., "random_split", "spatial_split")
#' @param fold_indices vector of fold assignments for each observation (or list of test indices per fold)
#' @param core_data data frame containing predictor variables and response (median_carbon_density)
#' @param predictor_vars character vector of predictor variable names
#' @param tune_hyperparams logical, whether to tune hyperparameters (default: false)
#' @param nested_tuning logical, if true uses nested cv (tune per fold), if false tunes once on first fold (default: true)
#' @param verbose logical, whether to print detailed progress during tuning (default: false)
#' @param buffer_m optional; if set with core_sf, training set = points not in test and > buffer_m from all test points (distance-buffered CV)
#' @param core_sf optional sf object of core locations (required when buffer_m is set; used for kriging if provided)
#' @param return_predictions if TRUE, return list(metrics = ..., predictions = ...) with observed/predicted/longitude/latitude per test point
#' @param skip_inla if TRUE, do not run INLA (avoids slow/hanging fits when INLA is slow or problematic)
#' @param models optional character vector of model names to run; NULL = all models. Names: "Random Forest", "XGBoost", "SVM", "Neural Network", "GPR", "Ordinary Kriging", "Universal Kriging (env drift)", "GAM", "INLA".
#' @param predictor_vars_by_model optional named list: model -> character vector of vars; when set, each model uses its own vars instead of predictor_vars
#' @param log_response if TRUE (default from \code{log_transform_target} in .GlobalEnv), fit on log(response) and back-transform predictions for metrics
#' @return data frame with performance metrics (r2, rmse, mae, bias) for each model and fold; or list(metrics, predictions) if return_predictions
run_cv <- function(cv_method_name, fold_indices, core_data, predictor_vars,
                   verbose = FALSE, buffer_m = NA_real_, core_sf = NULL,
                   return_predictions = FALSE, models = NULL,
                   tune_hyperparams = FALSE, nested_tuning = FALSE,
                   predictor_vars_by_model = NULL,
                   hyperparams_by_model = NULL,
                   log_response = get0("log_transform_target", envir = .GlobalEnv, ifnotfound = TRUE)) {
  models <- intersect(models %||% c("GPR", "GAM", "XGB", "RF"),
                      c("GPR", "GAM", "XGB", "RF"))
  if (length(models) == 0L) stop("run_cv: no valid models selected")
  if (!"median_carbon_density" %in% names(core_data))
    stop("run_cv: core_data must contain 'median_carbon_density' (or set it from target_var before calling)")

  # Normalise fold_indices to integer vector
  if (is.list(fold_indices)) {
    fv <- integer(nrow(core_data))
    for (k in seq_along(fold_indices)) fv[fold_indices[[k]]] <- k
    fold_indices <- fv
  }
  fold_indices <- as.integer(fold_indices)
  stopifnot(length(fold_indices) == nrow(core_data))

  use_buffer <- !is.na(buffer_m) && buffer_m > 0 && inherits(core_sf, "sf")
  cat("=== CV:", cv_method_name, "(", paste(models, collapse = ", "), ") ===\n")

  results_list <- list()
  preds_list   <- list()
  n_folds      <- max(fold_indices)

  for (fold in seq_len(n_folds)) {
    test_idx  <- which(fold_indices == fold)
    train_idx <- which(fold_indices != fold)
    if (use_buffer && length(test_idx) > 0L && length(train_idx) > 0L) {
      dists <- sf::st_distance(core_sf[train_idx, ], core_sf[test_idx, ])
      train_idx <- train_idx[apply(dists, 1, min) > buffer_m]
    }
    if (length(test_idx) == 0L || length(train_idx) == 0L) next

    train_raw <- as.data.frame(core_data[train_idx, ])
    test_raw  <- as.data.frame(core_data[test_idx, ])
    observed_orig <- test_raw$median_carbon_density
    if (log_response) {
      train_raw <- transform_response(train_raw, "median_carbon_density", log = TRUE)
      test_raw  <- transform_response(test_raw, "median_carbon_density", log = TRUE)
    }

    for (m in models) {
      pvars <- if (is.list(predictor_vars_by_model) && m %in% names(predictor_vars_by_model))
        predictor_vars_by_model[[m]] else predictor_vars
      res <- tryCatch({
        prep <- prepare_data_for_model(m, train_raw, test_raw, pvars)
        hp <- if (is.list(hyperparams_by_model) && m %in% names(hyperparams_by_model))
          hyperparams_by_model[[m]] else NULL
        fit <- switch(m,
          GPR = fit_gpr(
            prep$train, prep$predictor_vars,
            test_data = prep$test,
            hyperparams = hp
          ),
          RF  = fit_rf(
            prep$train, prep$test, prep$predictor_vars,
            hyperparams = hp
          ),
          XGB = fit_xgboost(
            prep$train, prep$test, prep$predictor_vars,
            hyperparams = hp
          ),
          GAM = fit_gam(
            prep$train, prep$test, prep$predictor_vars,
            include_spatial = FALSE,
            k_covariate = if (is.list(hp) && !is.null(hp$k_covariate)) hp$k_covariate else 6L
          ),
          stop("Unsupported model: ", m)
        )
        if (all(is.na(fit$predictions))) return(NULL)
        pred <- fit$predictions
        if (log_response) pred <- inverse_response_transform(pred, log = TRUE)
        keep <- test_rows_with_factors_in_train(train_raw, test_raw, intersect(pvars, c("seagrass_species", "region")))
        obs_sub <- observed_orig[keep]
        pred_sub <- pred[keep]
        metrics <- if (sum(keep) < 2L) data.frame(r2 = NA_real_, rmse = NA_real_, mae = NA_real_, bias = NA_real_) else calculate_metrics(obs_sub, pred_sub)
        list(
          metrics = data.frame(method = cv_method_name, fold = fold, model = m, metrics),
          preds   = if (return_predictions) data.frame(
            method = cv_method_name, fold = fold, model = m,
            observed = observed_orig, predicted = pred,
            longitude = test_raw$longitude, latitude = test_raw$latitude,
            stringsAsFactors = FALSE
          )
        )
      }, error = function(e) {
        cat("    ", m, " fold ", fold, ": ", e$message, "\n", sep = "")
        NULL
      })
      if (is.null(res)) next
      results_list[[length(results_list) + 1]] <- res$metrics
      if (!is.null(res$preds)) preds_list[[length(preds_list) + 1]] <- res$preds
    }
  }

  results_df <- dplyr::bind_rows(results_list)
  if (return_predictions && length(preds_list) > 0)
    return(list(metrics = results_df, predictions = dplyr::bind_rows(preds_list)))
  results_df
}

# ================================ APPLICABILITY DOMAIN ================================

#' Compute environmental similarity scores between prediction locations and training data.
#'
#' @param training_data Data frame with training observations and predictor variables.
#' @param prediction_data Data frame with prediction locations and predictor variables.
#' @param predictor_vars Character vector of predictor variable names.
#' @param method Similarity metric: "euclidean" (scaled), "mahalanobis", or "range_overlap".
#' @param verbose If TRUE, print diagnostic information.
#' @return List with similarity_scores, flag_outside_range, summary, method, and data_similarity_scores.
compute_environmental_similarity <- function(training_data, prediction_data, predictor_vars,
                                            method = "euclidean", verbose = FALSE) {
  pv <- intersect(intersect(predictor_vars, names(training_data)), names(prediction_data))
  if (length(pv) == 0) {
    warning("No common predictor variables found")
    return(list(similarity_scores = NULL, flag_outside_range = NULL, summary = NULL,
                data_similarity_scores = NULL))
  }

  X_train <- as.matrix(training_data[, pv, drop = FALSE])
  X_pred <- as.matrix(prediction_data[, pv, drop = FALSE])
  X_train[!is.finite(X_train)] <- NA_real_
  X_pred[!is.finite(X_pred)] <- NA_real_

  train_complete <- stats::complete.cases(X_train)
  n_train_complete <- sum(train_complete)
  if (verbose) {
    cat("Environmental similarity: ", nrow(training_data), " train (", n_train_complete, " complete), ",
        nrow(prediction_data), " pred, ", length(pv), " predictors\n", sep = "")
  }
  if (n_train_complete == 0) {
    warning("No complete cases in training data for similarity computation")
    n_pred <- nrow(prediction_data)
    return(list(similarity_scores = rep(NA_real_, n_pred), flag_outside_range = rep(FALSE, n_pred), summary = NULL,
                data_similarity_scores = rep(NA_real_, n_train_complete)))
  }

  X_train_complete <- X_train[train_complete, , drop = FALSE]
  train_cols_valid <- colSums(!is.na(X_train_complete)) > 0
  if (sum(train_cols_valid) == 0) {
    warning("No valid predictor columns in training data")
    n_pred <- nrow(prediction_data)
    return(list(similarity_scores = rep(NA_real_, n_pred), flag_outside_range = rep(FALSE, n_pred), summary = NULL,
                data_similarity_scores = rep(NA_real_, nrow(X_train_complete))))
  }
  X_train_complete <- X_train_complete[, train_cols_valid, drop = FALSE]
  X_pred <- X_pred[, train_cols_valid, drop = FALSE]
  n_train <- nrow(X_train_complete)
  n_pred <- nrow(X_pred)

  # Helper: convert distance vector to similarity in [0,1]
  dist_to_sim <- function(d) {
    finite <- is.finite(d)
    if (!any(finite)) return(list(scores = rep(NA_real_, length(d)), finite = finite))
    max_d <- max(d[finite], na.rm = TRUE)
    if (!is.finite(max_d) || max_d <= 0) return(list(scores = rep(NA_real_, length(d)), finite = finite))
    s <- rep(NA_real_, length(d))
    s[finite] <- pmax(0, pmin(1, 1 - (d[finite] / (max_d + 1e-10))))
    list(scores = s, finite = finite)
  }

  # Chunk size: avoid large matrices. D_sq is chunk_size * n_train * 8 bytes; keep ~50MB max per chunk.
  chunk_size <- min(5000L, max(1000L, as.integer(5e7 / (n_train * 8L))))

  similarity_scores <- NULL
  data_similarity_scores <- NULL

  if (method == "euclidean") {
    train_mean <- colMeans(X_train_complete, na.rm = TRUE)
    train_sd <- apply(X_train_complete, 2, sd, na.rm = TRUE)
    train_sd[train_sd == 0 | !is.finite(train_sd)] <- 1
    train_mean[!is.finite(train_mean)] <- 0
    X_train_scaled <- scale(X_train_complete, center = train_mean, scale = train_sd)
    X_train_scaled[!is.finite(X_train_scaled)] <- 0
    train_sq <- rowSums(X_train_scaled^2)
    pred_all_na <- rowSums(is.finite(X_pred)) == 0L

    # For predictions
    min_dists <- numeric(n_pred)
    for (start in seq(1L, n_pred, by = chunk_size)) {
      end <- min(start + chunk_size - 1L, n_pred)
      P <- scale(X_pred[start:end, , drop = FALSE], center = train_mean, scale = train_sd)
      P[!is.finite(P)] <- 0
      D_sq <- rowSums(P^2) + matrix(train_sq, nrow(P), n_train, byrow = TRUE) - 2 * (P %*% t(X_train_scaled))
      min_dists[start:end] <- sqrt(pmax(apply(D_sq, 1L, min, na.rm = TRUE), 0))
    }
    min_dists[pred_all_na] <- NA_real_
    out <- dist_to_sim(min_dists)
    similarity_scores <- out$scores

    # For training data points (self-similarity / environmental similarity to other data points)
    # Compute for each train point the minimum distance to any other (leave-one-out)
    train_min_dists <- numeric(n_train)
    if (n_train > 1) {
      # Precompute full distance matrix (can chunk if n_train too big)
      Dtt_sq <- matrix(NA_real_, n_train, n_train)
      for (i in seq(1L, n_train, by = chunk_size)) {
        i_end <- min(i + chunk_size - 1L, n_train)
        x1 <- X_train_scaled[i:i_end, , drop = FALSE]
        D_chunk <- rowSums(x1^2) + matrix(train_sq, nrow(x1), n_train, byrow = TRUE) - 2 * (x1 %*% t(X_train_scaled))
        diag(D_chunk) <- Inf # exclude self
        train_min_dists[i:i_end] <- sqrt(pmax(apply(D_chunk, 1L, min, na.rm = TRUE), 0))
      }
    } else {
      train_min_dists[] <- NA_real_
    }
    train_out <- dist_to_sim(train_min_dists)
    data_similarity_scores <- train_out$scores
  } else if (method == "mahalanobis") {
    train_mean <- colMeans(X_train_complete, na.rm = TRUE)
    train_cov <- cov(X_train_complete, use = "complete.obs")
    train_mean[!is.finite(train_mean)] <- 0
    train_cov[!is.finite(train_cov)] <- 0
    if (det(train_cov) < 1e-10 || !is.finite(det(train_cov))) {
      if (verbose) cat("  Covariance singular, falling back to euclidean\n")
      return(compute_environmental_similarity(training_data, prediction_data, predictor_vars, "euclidean", verbose))
    }
    Sigma_inv <- solve(train_cov)
    d_sq <- numeric(n_pred)
    for (start in seq(1L, n_pred, by = chunk_size)) {
      end <- min(start + chunk_size - 1L, n_pred)
      X_chunk <- X_pred[start:end, , drop = FALSE]
      X_cent <- sweep(X_chunk, 2L, train_mean, "-")
      X_cent[!is.finite(X_cent)] <- 0
      d_sq[start:end] <- rowSums((X_cent %*% Sigma_inv) * X_cent)
    }
    similarity_scores <- dist_to_sim(sqrt(pmax(d_sq, 0)))$scores

    # training data similarity (leave-one-out Mahalanobis)
    train_d_sq <- numeric(n_train)
    if (n_train > 1) {
      for (start in seq(1L, n_train, by = chunk_size)) {
        end <- min(start + chunk_size - 1L, n_train)
        X_chunk <- X_train_complete[start:end, , drop = FALSE]
        X_cent <- sweep(X_chunk, 2L, train_mean, "-")
        X_cent[!is.finite(X_cent)] <- 0
        # For each row, calculate Mahalanobis to all others and take min (but exclude self)
        D_chunk <- matrix(NA_real_, nrow(X_chunk), n_train)
        for (j in seq_len(nrow(X_chunk))) {
          diffs <- t(X_train_complete) - X_chunk[j, ]
          D_chunk[j, ] <- sqrt(colSums((Sigma_inv %*% diffs) * diffs))
          D_chunk[j, start + j - 1L] <- Inf # exclude self
        }
        train_d_sq[start:end] <- apply(D_chunk, 1L, min, na.rm = TRUE)
      }
    } else {
      train_d_sq[] <- NA_real_
    }
    data_similarity_scores <- dist_to_sim(train_d_sq)$scores
  } else if (method == "range_overlap") {
    train_ranges <- apply(X_train_complete, 2L, range, na.rm = TRUE)
    train_ranges[!is.finite(train_ranges)] <- 0
    similarity_scores <- numeric(n_pred)
    for (start in seq(1L, n_pred, by = chunk_size)) {
      end <- min(start + chunk_size - 1L, n_pred)
      X_chunk <- X_pred[start:end, , drop = FALSE]
      in_range <- (X_chunk >= train_ranges[1, ]) & (X_chunk <= train_ranges[2, ])
      n_valid <- rowSums(is.finite(X_chunk))
      vec <- rowSums(in_range & is.finite(X_chunk), na.rm = TRUE) / pmax(n_valid, 1L)
      vec[n_valid == 0L] <- NA_real_
      similarity_scores[start:end] <- vec
    }
    # For training data, all should be in range, so set to 1 if finite
    data_similarity_scores <- rep(1, nrow(X_train_complete))
    data_similarity_scores[!is.finite(rowSums(X_train_complete))] <- NA_real_
  } else {
    stop("Unknown method: ", method)
  }

  finite_scores <- is.finite(similarity_scores)
  flag_outside_range <- finite_scores & (similarity_scores < 0.5)
  n_finite <- sum(finite_scores)
  summary_stats <- list(
    mean_similarity = if (n_finite > 0) mean(similarity_scores[finite_scores], na.rm = TRUE) else NA_real_,
    min_similarity = if (n_finite > 0) min(similarity_scores[finite_scores], na.rm = TRUE) else NA_real_,
    max_similarity = if (n_finite > 0) max(similarity_scores[finite_scores], na.rm = TRUE) else NA_real_,
    n_outside_range = sum(flag_outside_range, na.rm = TRUE),
    pct_outside_range = 100 * mean(flag_outside_range, na.rm = TRUE),
    n_finite_scores = n_finite,
    n_total_pred = n_pred,
    pct_finite = 100 * n_finite / n_pred
  )
  if (verbose && n_finite > 0) {
    cat("  Similarity range: ", round(summary_stats$min_similarity, 3), " - ", round(summary_stats$max_similarity, 3),
        ", ", summary_stats$n_outside_range, " outside range\n", sep = "")
  }
  list(
    similarity_scores = similarity_scores,
    flag_outside_range = flag_outside_range,
    summary = summary_stats,
    method = method,
    data_similarity_scores = data_similarity_scores
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
#' @param plot_data_points If TRUE, plot original data points (optional).
#' @param dat Optional data frame with training data points to overlay.
#' @param flag_outside Logical vector indicating outside-domain locations (optional).
#' @param coords Character vector for coordinate columns (default c("longitude", "latitude")).
#' @param world Optional world map data. If NULL, map is retrieved.
#' @param use_raster If TRUE, use geom_raster; else use geom_point.
#' @param point_size Size of points when use_raster = FALSE (default 0.3).
#' @param xlim, ylim Optional axis limits.
#' @return ggplot2 object showing environmental similarity or applicability domain.
plot_applicability_domain <- function(
  prediction_grid,
  similarity_scores = NULL,
  flag_outside = NULL,
  plot_data_points = FALSE,
  dat = NULL,
  data_similarity_scores = NULL,
  coords = c("longitude", "latitude"),
  world = NULL,
  use_raster = FALSE,
  point_size = 0.3,
  xlim = NULL,
  ylim = NULL
) {
  stopifnot(requireNamespace("ggplot2", quietly = TRUE))

  # Ensure world map is a plain data frame; if we can't coerce safely, drop it
  if (is.null(world) && requireNamespace("maps", quietly = TRUE)) {
    world <- ggplot2::map_data("world")
  }
  if (!is.null(world)) {
    world_df <- tryCatch(as.data.frame(world), error = function(e) NULL)
    if (is.null(world_df) || length(dim(world_df)) != 2L) {
      warning("world object is not data.frame-like; skipping background polygons.")
      world <- NULL
    } else {
      world <- world_df
    }
  }

  df <- data.frame(prediction_grid[, coords, drop = FALSE], stringsAsFactors = FALSE)

  if (!is.null(similarity_scores)) df$similarity <- similarity_scores
  if (!is.null(flag_outside)) df$outside_domain <- flag_outside

  # Filter for finite/defined values
  if (!is.null(similarity_scores)) df <- df[is.finite(df$similarity), , drop = FALSE]
  if (!is.null(flag_outside)) df <- df[!is.na(df$outside_domain), , drop = FALSE]

  p <- ggplot2::ggplot()

  # 1) Prediction-grid similarity or applicability domain
  if (!is.null(similarity_scores) && nrow(df) > 0) {
    if (use_raster) {
      p <- p + ggplot2::geom_raster(
        data = df,
        ggplot2::aes(x = .data[[coords[1]]], y = .data[[coords[2]]], fill = similarity),
        na.rm = TRUE
      )
    } else {
      p <- p + ggplot2::geom_point(
        data = df,
        ggplot2::aes(x = .data[[coords[1]]], y = .data[[coords[2]]], fill = similarity),
        shape = 21,
        colour = "black",
        stroke = 0.2,
        size = point_size,
        alpha = 0.7,
        na.rm = TRUE
      )
    }
    # Single shared colour scale for similarity (raster + data points)
    p <- p + ggplot2::scale_fill_viridis_c(option = "plasma", name = "Similarity", na.value = "transparent")
    p <- p + ggplot2::labs(
      x = coords[1], y = coords[2],
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
          name = "Outside\nDomain", na.value = "transparent"
        )
    } else {
      p <- p + ggplot2::geom_point(
        data = df,
        ggplot2::aes(x = .data[[coords[1]]], y = .data[[coords[2]]], color = outside_domain),
        size = point_size, alpha = 0.7, na.rm = TRUE
      ) +
        ggplot2::scale_color_manual(
          values = c("FALSE" = "green", "TRUE" = "red"),
          name = "Outside\nDomain", na.value = "transparent"
        )
    }
    p <- p + ggplot2::labs(
      x = coords[1], y = coords[2],
      title = "Applicability domain (red = outside training range)"
    )
  } else {
    stop("Provide either similarity_scores or flag_outside with valid values.")
  }

  # plot worldmap
  required_cols <- c("long", "lat", "group")
  if (!is.null(world) && all(required_cols %in% names(world))) {
    p <- p + ggplot2::geom_polygon(
      data = world,
      ggplot2::aes(x = long, y = lat, group = group),
      fill = "grey90", color = "grey70", alpha = 1, inherit.aes = FALSE
    )
  }


  # overlay original data points (training data), coloured by their similarity
  if (plot_data_points && !is.null(dat)) {
    df_dat <- data.frame(dat[, coords, drop = FALSE], stringsAsFactors = FALSE)
    if (!is.null(data_similarity_scores)) {
      # Align length defensively
      if (length(data_similarity_scores) != nrow(df_dat)) {
        len <- min(length(data_similarity_scores), nrow(df_dat))
        df_dat <- df_dat[seq_len(len), , drop = FALSE]
        data_similarity_scores <- data_similarity_scores[seq_len(len)]
      }
      df_dat$similarity <- data_similarity_scores
      df_dat <- df_dat[is.finite(df_dat$similarity), , drop = FALSE]
      if (nrow(df_dat) > 0) {
        p <- p + ggplot2::geom_point(
          data = df_dat,
          ggplot2::aes(
            x = .data[[coords[1]]],
            y = .data[[coords[2]]],
            fill = similarity
          ),
          shape = 21,
          colour = "black",
          stroke = 0.2,
          size = point_size * 1.8,
          alpha = 0.9,
          na.rm = TRUE
        )
      }
    } else {
      # Fallback: uncoloured points if no similarity for data was supplied
      p <- p + ggplot2::geom_point(
        data = df_dat,
        ggplot2::aes(x = .data[[coords[1]]], y = .data[[coords[2]]]),
        color = "black", size = point_size * 1.5, alpha = 0.7, na.rm = TRUE
      )
    }
  }

  # Axis limits
  if (!is.null(xlim) || !is.null(ylim)) {
    p <- p + ggplot2::coord_cartesian(xlim = xlim, ylim = ylim)
  }

  p + ggplot2::theme_minimal()
}
# Load ML helpers after all utility functions are defined (avoids circular source with ml.R)
source("modelling/R/ml.R")
