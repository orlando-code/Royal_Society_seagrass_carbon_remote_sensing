# Compare prediction-location environmental vectors to the model training set.
# Used by scripts in review/; relies on helpers from modelling/R/helpers.R.

#' @param training_data Training observations (e.g. all_extracted_new.rds contents).
#' @param prediction_data Locations with extracted predictor values.
#' @param predictor_vars Model predictor names.
#' @param output_dir Directory for CSV (and optional plot) outputs.
#' @param similarity_method Passed to compute_environmental_similarity().
#' @param similarity_threshold Scores below this are flagged as low similarity.
#' @param strict_range_flag If TRUE, flag when any predictor is outside training range.
#' @param make_plots Write diagnostic ggplot2 figures when TRUE.
#' @param point_id_cols Columns used to label points in outputs.
compare_prediction_env_to_training <- function(
  training_data,
  prediction_data,
  predictor_vars,
  output_dir,
  similarity_method = "euclidean",
  similarity_threshold = 0.5,
  strict_range_flag = TRUE,
  make_plots = FALSE,
  point_id_cols = c("longitude", "latitude", "seagrass_species"),
  verbose = TRUE
) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  env_vars <- setdiff(
    predictor_vars,
    c("longitude", "latitude", "seagrass_species")
  )
  env_vars <- intersect(env_vars, names(training_data))
  env_vars <- intersect(env_vars, names(prediction_data))
  if (length(env_vars) == 0L) {
    stop("No numeric environmental predictors available for training comparison.")
  }

  training_env <- process_rs_covariates(training_data)

  if (verbose) {
    cat(
      "\nComparing ", nrow(prediction_data), " prediction location(s) to ",
      nrow(training_env), " training observation(s) across ",
      length(env_vars), " environmental predictor(s)...\n",
      sep = ""
    )
  }

  similarity <- compute_environmental_similarity(
    training_data = training_env,
    prediction_data = prediction_data,
    predictor_vars = env_vars,
    method = similarity_method,
    verbose = verbose
  )

  outside_strict <- flag_outside_applicability_domain(
    training_env,
    prediction_data,
    env_vars,
    strict = TRUE
  )
  outside_relaxed <- flag_outside_applicability_domain(
    training_env,
    prediction_data,
    env_vars,
    strict = strict_range_flag
  )

  train_ranges <- stats::setNames(
    lapply(env_vars, function(v) range(training_env[[v]], na.rm = TRUE)),
    env_vars
  )
  train_mat <- as.matrix(training_env[, env_vars, drop = FALSE])
  train_mean <- colMeans(train_mat, na.rm = TRUE)
  train_sd <- apply(train_mat, 2, stats::sd, na.rm = TRUE)
  train_sd[!is.finite(train_sd) | train_sd == 0] <- 1

  predictor_long <- vector("list", length = nrow(prediction_data) * length(env_vars))
  idx <- 1L
  n_outside <- integer(nrow(prediction_data))
  for (i in seq_len(nrow(prediction_data))) {
    outside_n <- 0L
    for (v in env_vars) {
      pred_val <- prediction_data[[v]][i]
      tr <- train_ranges[[v]]
      outside <- is.finite(pred_val) && (pred_val < tr[1] || pred_val > tr[2])
      if (isTRUE(outside)) {
        outside_n <- outside_n + 1L
      }
      predictor_long[[idx]] <- data.frame(
        prediction_row = i,
        variable = v,
        train_min = tr[1],
        train_max = tr[2],
        train_mean = train_mean[[v]],
        train_sd = train_sd[[v]],
        pred_value = pred_val,
        pred_zscore = (pred_val - train_mean[[v]]) / train_sd[[v]],
        outside_training_range = outside,
        flag_extreme_zscore = is.finite(pred_val) &&
          is.finite(train_mean[[v]]) &&
          abs((pred_val - train_mean[[v]]) / train_sd[[v]]) > 2,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
    n_outside[i] <- outside_n
  }
  predictor_long_df <- dplyr::bind_rows(predictor_long)

  id_cols <- intersect(point_id_cols, names(prediction_data))
  point_summary <- if (length(id_cols) > 0L) {
    prediction_data[, id_cols, drop = FALSE]
  } else {
    data.frame(prediction_row = seq_len(nrow(prediction_data)))
  }
  if (!"prediction_row" %in% names(point_summary)) {
    point_summary$prediction_row <- seq_len(nrow(prediction_data))
  }

  low_similarity <- is.finite(similarity$similarity_scores) &
    similarity$similarity_scores < similarity_threshold

  point_summary$env_similarity_score <- similarity$similarity_scores
  point_summary$env_flag_low_similarity <- low_similarity
  point_summary$env_flag_any_predictor_outside_range <- outside_strict
  point_summary$env_flag_majority_predictors_outside_range <- outside_relaxed
  point_summary$env_n_predictors_outside_range <- n_outside
  point_summary$env_n_predictors_compared <- length(env_vars)
  point_summary$env_pct_predictors_outside_range <-
    100 * n_outside / length(env_vars)

  if ("seagrass_species" %in% names(prediction_data) &&
      "seagrass_species" %in% names(training_env)) {
    train_species <- unique(as.character(training_env$seagrass_species))
    point_summary$env_flag_novel_species <-
      !as.character(prediction_data$seagrass_species) %in% train_species
  }

  point_summary$env_flag_extrapolation <-
    point_summary$env_flag_low_similarity |
    point_summary$env_flag_any_predictor_outside_range |
    if ("env_flag_novel_species" %in% names(point_summary)) {
      point_summary$env_flag_novel_species
    } else {
      FALSE
    }

  predictor_summary <- predictor_long_df %>%
    dplyr::group_by(.data$variable) %>%
    dplyr::summarise(
      train_min = dplyr::first(.data$train_min),
      train_max = dplyr::first(.data$train_max),
      train_mean = dplyr::first(.data$train_mean),
      train_sd = dplyr::first(.data$train_sd),
      pred_min = min(.data$pred_value, na.rm = TRUE),
      pred_max = max(.data$pred_value, na.rm = TRUE),
      n_points_outside_range = sum(.data$outside_training_range, na.rm = TRUE),
      pct_points_outside_range = 100 * mean(.data$outside_training_range, na.rm = TRUE),
      max_abs_zscore = max(abs(.data$pred_zscore), na.rm = TRUE),
      .groups = "drop"
    )

  if (length(id_cols) > 0L) {
    ids <- prediction_data[, id_cols, drop = FALSE]
    ids$prediction_row <- seq_len(nrow(prediction_data))
    predictor_long_df <- dplyr::left_join(predictor_long_df, ids, by = "prediction_row")
  }

  readr::write_csv(
    point_summary,
    file.path(output_dir, "env_training_similarity_by_point.csv")
  )
  readr::write_csv(
    predictor_long_df,
    file.path(output_dir, "env_training_predictor_ranges_by_point.csv")
  )
  readr::write_csv(
    predictor_summary,
    file.path(output_dir, "env_training_predictor_summary.csv")
  )

  if (verbose) {
    cat("\nEnvironmental comparison summary:\n")
    cat(
      "  Similarity (", similarity_method, "): ",
      "mean = ", sprintf("%.3f", similarity$summary$mean_similarity),
      ", min = ", sprintf("%.3f", similarity$summary$min_similarity),
      ", max = ", sprintf("%.3f", similarity$summary$max_similarity), "\n",
      sep = ""
    )
    cat(
      "  Points flagged low similarity (< ", similarity_threshold, "): ",
      sum(low_similarity, na.rm = TRUE), " / ", nrow(point_summary), "\n",
      sep = ""
    )
    cat(
      "  Points with any predictor outside training range: ",
      sum(outside_strict, na.rm = TRUE), " / ", nrow(point_summary), "\n",
      sep = ""
    )
    if ("env_flag_novel_species" %in% names(point_summary)) {
      cat(
        "  Points with species not seen in training: ",
        sum(point_summary$env_flag_novel_species, na.rm = TRUE),
        " / ", nrow(point_summary), "\n",
        sep = ""
      )
    }
    cat(
      "  Points flagged for extrapolation (any criterion): ",
      sum(point_summary$env_flag_extrapolation, na.rm = TRUE),
      " / ", nrow(point_summary), "\n",
      sep = ""
    )

    flagged <- point_summary[point_summary$env_flag_extrapolation, , drop = FALSE]
    if (nrow(flagged) > 0L) {
      cat("\n  Flagged prediction location(s):\n")
      show_cols <- intersect(
        c(id_cols, "prediction_row", "env_similarity_score",
          "env_n_predictors_outside_range", "env_flag_novel_species"),
        names(flagged)
      )
      print(flagged[, show_cols, drop = FALSE], row.names = FALSE)
      outside_vars <- predictor_long_df %>%
        dplyr::filter(.data$outside_training_range) %>%
        dplyr::select(dplyr::any_of(c(id_cols, "prediction_row", "variable",
                                      "train_min", "train_max", "pred_value")))
      if (nrow(outside_vars) > 0L) {
        cat("\n  Predictor(s) outside training range:\n")
        print(outside_vars, row.names = FALSE)
      }
    }

    cat("\nWrote environmental comparison tables to:\n  ", output_dir, "\n", sep = "")
  }

  if (isTRUE(make_plots) && requireNamespace("ggplot2", quietly = TRUE)) {
    clamp_plot_inches <- function(x, min_val = 4, max_val = 48) {
      pmax(min_val, pmin(max_val, x))
    }
    plot_height_for_n_points <- function(n) {
      row_in <- if (n > 100L) 0.12 else if (n > 50L) 0.2 else 0.35
      clamp_plot_inches(row_in * n + 2)
    }
    plot_width_for_n_vars <- function(n) {
      clamp_plot_inches(0.45 * n + 4, min_val = 8)
    }

    label_col <- if ("Site" %in% names(point_summary)) {
      "Site"
    } else if (all(c("longitude", "latitude") %in% names(point_summary))) {
      NULL
    } else {
      "prediction_row"
    }

    plot_df <- point_summary
    if (!is.null(label_col)) {
      plot_df$point_label <- as.character(plot_df[[label_col]])
    } else if (all(c("longitude", "latitude") %in% names(plot_df))) {
      plot_df$point_label <- sprintf(
        "(%.3f, %.3f)",
        plot_df$longitude,
        plot_df$latitude
      )
    } else {
      plot_df$point_label <- as.character(plot_df$prediction_row)
    }
    plot_df$point_label <- factor(
      plot_df$point_label,
      levels = rev(unique(plot_df$point_label))
    )

    p_sim <- ggplot2::ggplot(
      plot_df,
      ggplot2::aes(
        x = .data$env_similarity_score,
        y = .data$point_label,
        fill = .data$env_flag_extrapolation
      )
    ) +
      ggplot2::geom_col(width = 0.7) +
      ggplot2::geom_vline(
        xintercept = similarity_threshold,
        linetype = "dashed",
        colour = "firebrick"
      ) +
      ggplot2::scale_fill_manual(
        values = c("TRUE" = "#D6604D", "FALSE" = "#2166AC"),
        labels = c("TRUE" = "Extrapolation flagged", "FALSE" = "Within training envelope"),
        name = NULL
      ) +
      ggplot2::labs(
        title = "Environmental similarity to training data",
        subtitle = paste0(
          "Dashed line = similarity threshold (", similarity_threshold, ")"
        ),
        x = paste0("Similarity score (", similarity_method, ")"),
        y = NULL
      ) +
      ggplot2::theme_minimal()

    ggplot2::ggsave(
      file.path(output_dir, "env_training_similarity_by_point.png"),
      p_sim,
      width = 8,
      height = plot_height_for_n_points(nrow(plot_df)),
      dpi = 150,
      limitsize = FALSE
    )

    z_plot_df <- predictor_long_df
    if (!"point_label" %in% names(z_plot_df)) {
      z_plot_df <- dplyr::left_join(
        z_plot_df,
        plot_df[, c("prediction_row", "point_label"), drop = FALSE],
        by = "prediction_row"
      )
    }
    z_plot_df$point_label <- factor(
      z_plot_df$point_label,
      levels = levels(plot_df$point_label)
    )

    p_z <- ggplot2::ggplot(
      z_plot_df,
      ggplot2::aes(
        x = .data$variable,
        y = .data$point_label,
        fill = .data$pred_zscore
      )
    ) +
      ggplot2::geom_tile(colour = "white", linewidth = 0.2) +
      ggplot2::scale_fill_gradient2(
        low = "#2166AC",
        mid = "grey95",
        high = "#D6604D",
        midpoint = 0,
        name = "Z-score vs training"
      ) +
      ggplot2::labs(
        title = "Environmental predictors relative to training distribution",
        subtitle = "Values are standardised using training mean and SD",
        x = "Predictor",
        y = NULL
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

    ggplot2::ggsave(
      file.path(output_dir, "env_training_predictor_zscores.png"),
      p_z,
      width = plot_width_for_n_vars(length(env_vars)),
      height = plot_height_for_n_points(nrow(plot_df)),
      dpi = 150,
      limitsize = FALSE
    )

    if (verbose) {
      cat(
        "  Saved plots: env_training_similarity_by_point.png, ",
        "env_training_predictor_zscores.png\n",
        sep = ""
      )
    }
  }

  list(
    point_summary = point_summary,
    predictor_summary = predictor_summary,
    predictor_long = predictor_long_df,
    similarity = similarity
  )
}
