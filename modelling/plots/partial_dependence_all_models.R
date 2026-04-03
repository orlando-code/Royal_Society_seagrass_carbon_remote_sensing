# Partial Dependence Plots for XGB, GAM, GPR, and LR final models
#
# Uses iml::FeatureEffect (method = "pdp") to compute and plot the marginal
# effect of each predictor for each model's final fitted object.
#
# Requires: output/<cv_regime>/final_models/<model>_final.rds (from fit_final_models.R)
# Outputs:  output/<cv_regime>/covariate_selection/pdp_<model>.png
#
# Usage: sourced from run_paper.R (step 5b), or run standalone.

setwd(here::here())
source("modelling/R/helpers.R")
source("modelling/R/plot_config.R")
source("modelling/R/extract_covariates_from_rasters.R")
source("modelling/R/assign_region_from_latlon.R")
source("modelling/config/pipeline_config.R")
cfg <- get_pipeline_config()
apply_pipeline_defaults(
  cfg,
  c(
    "cv_output_dir", "target_var", "exclude_regions", "model_list", "n_lons", "n_lats", "dpi", "show_titles", "prediction_map_use_log_fill", "prediction_map_log_fill_epsilon"
  ),
  envir = .GlobalEnv
)

load_packages(c("here", "dplyr", "ggplot2", "patchwork"))
if (!requireNamespace("iml", quietly = TRUE))
  stop("Package 'iml' is required; install with install.packages('iml').")

# ---------------------------------------------------------------------------
# Rebuild core_data (align with fit_final_models.R)
# ---------------------------------------------------------------------------
dat <- readr::read_rds("data/all_extracted_new.rds")
target_var        <- "median_carbon_density_100cm"
predictor_vars_all <- raster_covariates[raster_covariates %in% colnames(dat)]

excl <- get0("exclude_regions", envir = .GlobalEnv, ifnotfound = character(0))
if (length(excl) > 0L) {
  if (!"region" %in% names(dat)) dat <- assign_region_from_latlon(dat)
  dat <- dat[is.na(dat$region) | !dat$region %in% excl, ]
}

core_data <- dat %>%
  dplyr::mutate(median_carbon_density = .data[[target_var]]) %>%
  dplyr::select(
    longitude, latitude, median_carbon_density,
    dplyr::all_of(predictor_vars_all),
    dplyr::all_of(intersect(c("seagrass_species", "region"), names(.)))
  ) %>%
  dplyr::filter(complete.cases(.)) %>%
  as.data.frame()

cv_out  <- get0("cv_output_dir", envir = .GlobalEnv, ifnotfound = "output")
out_dir <- file.path(cv_out, "covariate_selection")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Build iml::Predictor for a final model
# ---------------------------------------------------------------------------
make_predictor <- function(model_name) {
  rds_path <- file.path(cv_out, "final_models",
                        paste0(model_name, "_final.rds"))
  if (!file.exists(rds_path)) {
    cat("Skipping", model_name, ": no final model at", rds_path, "\n")
    return(NULL)
  }
  obj   <- readRDS(rds_path)
  pvars <- obj$predictor_vars
  pvars <- intersect(pvars, names(core_data))
  if (length(pvars) == 0L) {
    cat("Skipping", model_name, ": no overlap between predictor_vars and core_data.\n")
    return(NULL)
  }

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
  } else if (model_name == "LR") {
    pred_fun <- function(m, newdata) {
      as.numeric(stats::predict(m, newdata = as.data.frame(newdata)))
    }
    iml::Predictor$new(
      model            = obj$model,
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
get_pdp_vars_for_model <- function(model_name, max_vars = 6L) {
  per_path <- file.path(cv_out, "covariate_selection",
    paste0("pruned_variables_to_include_", tolower(model_name), ".csv"))
  rds_path <- file.path(cv_out, "final_models", paste0(model_name, "_final.rds"))
  if (!file.exists(rds_path)) return(character(0))
  obj <- readRDS(rds_path)
  pvars <- if (model_name == "XGB" && !is.null(obj$encoded_names)) obj$encoded_names else obj$predictor_vars
  if (is.null(pvars) || length(pvars) == 0L) return(character(0))
  pvars <- intersect(pvars, names(core_data))
  if (file.exists(per_path)) {
    pv_df <- read.csv(per_path, stringsAsFactors = FALSE)
    pvars <- intersect(pv_df$variable, pvars)
  }
  pvars[seq_len(min(length(pvars), as.integer(max_vars)))]
}

compute_pdp_plots <- function(pred, pvars, model_name) {
  if (is.null(pred) || length(pvars) == 0L) return(list())
  lapply(pvars, function(v) {
    tryCatch({
      fe <- iml::FeatureEffect$new(pred, feature = v, method = "pdp", grid.size = 20)
      fe$plot() +
        ggplot2::labs(title = label_vars(v), y = "Predicted y") +
        ggplot2::theme_minimal(base_size = 10) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(size = 9),
          axis.title.x = ggplot2::element_blank()
        )
    }, error = function(e) {
      cat("    ! PDP error for", model_name, "::", v, ":", conditionMessage(e), "\n")
      NULL
    })
  })
}

plot_pdp_for_model <- function(model_name, max_vars = 6L) {
  pred <- make_predictor(model_name)
  if (is.null(pred)) return(invisible(NULL))
  pvars <- get_pdp_vars_for_model(model_name, max_vars = max_vars)
  cat("  Computing PDPs for", model_name, "(", length(pvars), "vars)...\n")
  plots <- compute_pdp_plots(pred, pvars, model_name)
  plots <- plots[!vapply(plots, is.null, logical(1))]
  if (length(plots) == 0) return(invisible(NULL))

  combined <- patchwork::wrap_plots(plots, ncol = 2) +
  if (show_titles) {
    patchwork::plot_annotation(
      title    = paste("Partial Dependence Plots –", model_name),
      subtitle = paste(length(pvars), "top predictors"),
      theme    = ggplot2::theme(
        plot.title    = ggplot2::element_text(size = 13, face = "bold"),
        plot.subtitle = ggplot2::element_text(size = 10)
      )
    )
  }

  out_path <- file.path(out_dir, paste0("pdp_", model_name, ".png"))
  ggplot2::ggsave(out_path, combined,
                  width = 10, height = ceiling(length(pvars) / 2) * 3,
                  dpi = dpi, limitsize = FALSE)
  cat("  Saved PDPs for", model_name, "->", out_path, "\n")
}

top_vars_by_model <- list()
all_vars_by_model <- list()
for (m in model_list) {
  # top 6 model covariates (for per-model plots)
  top_vars_by_model[[m]] <- get_pdp_vars_for_model(m, max_vars = 6L)
  # all model covariates (for shared-variable plots)
  all_vars_by_model[[m]] <- get_pdp_vars_for_model(m, max_vars = 42L)
  plot_pdp_for_model(m, max_vars = 6L)
}

predictors_by_model <- stats::setNames(
  lapply(model_list, make_predictor),
  model_list
)

# Shared-variable PDPs across models
shared_vars <- Reduce(intersect, all_vars_by_model)
if (length(shared_vars) > 0L) {
  cat("  Computing shared-variable PDPs (", length(shared_vars), " vars shared by configured models)...\n", sep = "")
  shared_plots <- list()
  n_models <- length(model_list)
  n_shared <- length(shared_vars)
  for (i in seq_along(shared_vars)) {
    v <- shared_vars[[i]]
    for (j in seq_along(model_list)) {
      m <- model_list[[j]]
      pred <- predictors_by_model[[m]]
      if (is.null(pred)) next
      p <- compute_pdp_plots(pred, v, m)[[1]]
      if (!is.null(p)) {
        # Top row, first column: title = model, subtitle = predictor. Other rows, first column: title = predictor.
        # Title/subtitle only on j == 1 (first plot per row); other columns have no panel title.
        vlab_display <- if (identical(v, "seagrass_species")) {
          "Seagrass species"
        } else {
          label_vars(v)
        }
        panel_title <- if (j == 1L) {
          if (i == 1L) m else vlab_display
        } else {
          NULL
        }
        panel_subtitle <- if (j == 1L && i == 1L) vlab_display else NULL
        # Bottom row: x-axis title = predictor; no x title for seagrass_species (discrete axis).
        panel_x <- if (i == n_shared && !identical(v, "seagrass_species")) vlab_display else NULL
        y_axis_theme <- if (j == 1L) {
          ggplot2::theme(
            axis.title.y = ggplot2::element_text(size = 7),
            axis.text.y = ggplot2::element_text(size = 6),
            axis.ticks.y = ggplot2::element_line()
          )
        } else {
          ggplot2::theme(
            axis.title.y = ggplot2::element_blank(),
            axis.text.y = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank()
          )
        }
        p_panel <- p +
          ggplot2::labs(
            title = panel_title,
            subtitle = panel_subtitle,
            x = panel_x,
            y = if (j == 1L) "Predicted y" else NULL
          ) +
          ggplot2::theme(
            plot.title = ggplot2::element_text(size = 8, face = "bold", hjust = 0),
            plot.subtitle = ggplot2::element_text(size = 7, hjust = 0),
            axis.title.x = ggplot2::element_text(size = 7),
            axis.text.x = ggplot2::element_text(size = 6),
            axis.text.y = ggplot2::element_text(size = 6),
            axis.ticks.x = ggplot2::element_line()
          ) +
          y_axis_theme
        if (identical(v, "seagrass_species")) {
          p_panel <- p_panel +
            ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(n.dodge = 2)) +
            ggplot2::theme(
              axis.text.x = ggplot2::element_text(size = 6, angle = 30, hjust = 1)
            )
        }
        shared_plots[[length(shared_plots) + 1L]] <- p_panel
      }
    }
  }
  if (length(shared_plots) > 0L) {
    combined_shared <- patchwork::wrap_plots(shared_plots, ncol = length(model_list))
    if (show_titles) {
      combined_shared <- combined_shared +
        patchwork::plot_annotation(
          title = "Partial Dependence - Shared Top Variables Across Models",
          subtitle = paste0("Shared predictors (n=", length(shared_vars), ")"),
          theme = ggplot2::theme(
            plot.title = ggplot2::element_text(size = 13, face = "bold"),
            plot.subtitle = ggplot2::element_text(size = 9)
          )
        )
    }
    dpi <- get0("dpi", envir = .GlobalEnv, ifnotfound = 150)
    out_shared <- file.path(out_dir, "pdp_shared_top_vars.png")
    ggplot2::ggsave(
      out_shared, combined_shared,
      width = 12, height = max(4, length(shared_vars) * 2.5),
      dpi = dpi, limitsize = FALSE
    )
    cat("  Saved shared-variable PDPs ->", out_shared, "\n")
  }
} else {
  cat("  No shared top variables across configured models; skipping shared-variable PDP figure.\n")
}

cat("Partial dependence plots complete.\n")
