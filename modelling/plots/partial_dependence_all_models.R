# Partial Dependence Plots for XGB, GAM, and GPR final models
#
# Uses iml::FeatureEffect (method = "pdp") to compute and plot the marginal
# effect of each predictor for each model's final fitted object.
#
# Requires: output/final_models/<model>_final.rds (from fit_final_models.R)
# Outputs:  output/covariate_selection/pdp_<model>.png
#
# Usage: sourced from run_paper_figures.R (step 5b), or run standalone.

setwd(here::here())
source("modelling/R/helpers.R")
# source("modelling/R/gpr_funs.R")
source("modelling/R/extract_covariates_from_rasters.R")
source("modelling/R/assign_region_from_latlon.R")
load_packages(c("here", "dplyr", "ggplot2", "patchwork"))
if (!requireNamespace("iml", quietly = TRUE))
  stop("Package 'iml' is required; install with install.packages('iml').")

# ---------------------------------------------------------------------------
# Rebuild core_data (same logic as fit_final_models.R)
# ---------------------------------------------------------------------------
dat <- readr::read_rds("data/all_extracted_new.rds")
target_var     <- "median_carbon_density_100cm"
predictor_vars_all <- raster_covariates[raster_covariates %in% colnames(dat)]

excl <- get0("exclude_regions", envir = .GlobalEnv, ifnotfound = character(0))
if (length(excl) > 0L) {
  if (!"region" %in% names(dat)) dat <- assign_region_from_latlon(dat)
  dat <- dat[is.na(dat$region) | !dat$region %in% excl, ]
}

core_data <- if ("random_core_variable" %in% colnames(dat)) {
  dat %>%
    dplyr::group_by(.data$random_core_variable) %>%
    dplyr::summarise(
      median_carbon_density = median(.data[[target_var]], na.rm = TRUE),
      dplyr::across(dplyr::all_of(c(predictor_vars_all, "longitude", "latitude")),
                    dplyr::first),
      .groups = "drop"
    ) %>%
    dplyr::filter(!is.na(.data$median_carbon_density))
} else {
  dat %>%
    dplyr::mutate(median_carbon_density = .data[[target_var]]) %>%
    dplyr::select(longitude, latitude, median_carbon_density,
                  dplyr::all_of(predictor_vars_all))
}
core_data <- as.data.frame(core_data[complete.cases(core_data), , drop = FALSE])

out_dir <- "output/covariate_selection"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Build iml::Predictor for a final model
# ---------------------------------------------------------------------------
make_predictor <- function(model_name) {
  rds_path <- file.path("output/final_models",
                        paste0(model_name, "_final.rds"))
  if (!file.exists(rds_path)) {
    cat("Skipping", model_name, ": no final model at", rds_path, "\n")
    return(NULL)
  }
  obj   <- readRDS(rds_path)
  pvars <- obj$predictor_vars

  if (model_name == "XGB") {
    sp <- obj$scale_params
    enc <- obj$encoding
    if (!is.null(enc)) {
      data_enc <- apply_categorical_encoding(core_data[, pvars, drop = FALSE], enc, pvars)
      X <- apply_scaling(data_enc, sp, obj$encoded_names %||% pvars)
      pvars_use <- obj$encoded_names %||% pvars
    } else {
      X <- apply_scaling(core_data[, pvars, drop = FALSE], sp, pvars)
      pvars_use <- pvars
    }
    pred_fun <- function(m, newdata) {
      nd <- apply_scaling(as.data.frame(newdata), sp, pvars_use)
      as.numeric(predict(m, newdata = as.matrix(nd[, pvars_use, drop = FALSE])))
    }
    iml::Predictor$new(
      model            = obj$model,
      data             = X[, pvars_use, drop = FALSE],
      y                = core_data$median_carbon_density,
      predict.function = pred_fun
    )

  } else if (model_name == "GAM") {
    sp <- obj$scale_params
    X_gam <- core_data[, c(pvars, "longitude", "latitude"), drop = FALSE]
    pred_fun <- function(m, newdata) {
      as.numeric(mgcv::predict.gam(m,
        newdata = apply_scaling(as.data.frame(newdata), sp, pvars),
        type = "response"))
    }
    iml::Predictor$new(
      model            = obj$model,
      data             = X_gam,
      y                = core_data$median_carbon_density,
      predict.function = pred_fun
    )

  } else if (model_name == "GPR") {
    # Reconstruct the gpr list that predict_gpr expects
    gpr_list <- list(
      model          = obj$model,
      predictor_vars = obj$predictor_vars,
      scale_params   = obj$scale_params,
      encoding       = obj$encoding,
      encoded_names  = obj$encoded_names
    )
    pred_fun <- function(m, newdata) {
      res <- predict_gpr(m, newdata[, pvars, drop = FALSE], se = FALSE)
      as.numeric(res$mean)
    }
    iml::Predictor$new(
      model            = gpr_list,
      data             = core_data[, pvars, drop = FALSE],
      y                = core_data$median_carbon_density,
      predict.function = pred_fun
    )

  } else {
    cat("Unsupported model for PDPs:", model_name, "\n")
    NULL
  }
}

# ---------------------------------------------------------------------------
# Plot PDPs for one model's top variables
# ---------------------------------------------------------------------------
plot_pdp_for_model <- function(model_name, max_vars = 6L) {
  pred <- make_predictor(model_name)
  if (is.null(pred)) return(invisible(NULL))

  per_path <- file.path("output/covariate_selection",
    paste0("pruned_variables_to_include_", tolower(model_name), ".csv"))
  obj <- readRDS(file.path("output/final_models", paste0(model_name, "_final.rds")))
  pvars <- if (model_name == "XGB" && !is.null(obj$encoded_names)) obj$encoded_names else obj$predictor_vars
  if (file.exists(per_path)) {
    pv_df <- read.csv(per_path, stringsAsFactors = FALSE)
    pvars <- intersect(pv_df$variable, pvars)
  }
  pvars <- pvars[seq_len(min(length(pvars), as.integer(max_vars)))]

  cat("  Computing PDPs for", model_name, "(", length(pvars), "vars)...\n")
  plots <- lapply(pvars, function(v) {
    tryCatch({
      fe <- iml::FeatureEffect$new(pred, feature = v, method = "pdp",
                                   grid.size = 20)
      fe$plot() +
        ggplot2::labs(title = v, x = v, y = "Predicted") +
        ggplot2::theme_minimal(base_size = 10) +
        ggplot2::theme(plot.title = ggplot2::element_text(size = 9))
    }, error = function(e) {
      cat("    ! PDP error for", v, ":", conditionMessage(e), "\n")
      NULL
    })
  })
  plots <- plots[!vapply(plots, is.null, logical(1))]
  if (length(plots) == 0) return(invisible(NULL))

  combined <- patchwork::wrap_plots(plots, ncol = 2) +
    patchwork::plot_annotation(
      title    = paste("Partial Dependence Plots –", model_name),
      subtitle = paste(length(pvars), "top predictors"),
      theme    = ggplot2::theme(
        plot.title    = ggplot2::element_text(size = 13, face = "bold"),
        plot.subtitle = ggplot2::element_text(size = 10)
      )
    )

  dpi_val  <- get0("dpi", envir = .GlobalEnv, ifnotfound = 150)
  out_path <- file.path(out_dir, paste0("pdp_", model_name, ".png"))
  ggplot2::ggsave(out_path, combined,
                  width = 10, height = ceiling(length(pvars) / 2) * 3,
                  dpi = dpi_val, limitsize = FALSE)
  cat("  Saved PDPs for", model_name, "->", out_path, "\n")
}

for (m in c("XGB", "GAM", "GPR")) {
  plot_pdp_for_model(m, max_vars = 6L)
}
cat("Partial dependence plots complete.\n")
