# Model permutation pruning: plot importance and kept variables from cached CSV
#
# Reads cv_model_permutation_pruning_*.csv (from model_permutation_pruning.R),
# plots permutation importance by model and variable, and which variables are kept.
# Usage: source("modelling/model_permutation_pruning_plots.R") or Rscript modelling/model_permutation_pruning_plots.R

if (requireNamespace("here", quietly = TRUE)) setwd(here::here())
if (!requireNamespace("here", quietly = TRUE)) install.packages("here", quiet = TRUE)
library(here)
library(ggplot2)
library(dplyr)
library(tidyr)

out_dir <- "figures/cv_pipeline_output"
search_dirs <- c("figures/covariate_selection")
if (exists("dpi", envir = .GlobalEnv)) dpi <- get("dpi", envir = .GlobalEnv) else dpi <- 150
if (exists("show_titles", envir = .GlobalEnv)) show_titles <- get("show_titles", envir = .GlobalEnv) else show_titles <- TRUE

# Find the permutation pruning CSV (any block size)
csv_path <- NULL
for (d in search_dirs) {
  if (!dir.exists(d)) next
  f <- list.files(d, pattern = "^pruned_model_variables*\\.csv$", full.names = TRUE)
  if (length(f) > 0L) {
    csv_path <- f[1L]
    break
  }
}
if (is.null(csv_path)) {
  stop("No pruned_model_variables*.csv found in ", paste(search_dirs, collapse = " or "),
       ". Run model_permutation_pruning.R first.")
}

cat("Reading:", csv_path, "\n")
imp <- read.csv(csv_path, stringsAsFactors = FALSE)
if (!all(c("model", "variable", "rmse_increase", "keep") %in% names(imp))) {
  stop("CSV must contain columns: model, variable, rmse_increase, keep. Got: ", paste(names(imp), collapse = ", "))
}

# Ensure keep is logical for plotting
imp$keep <- as.logical(imp$keep)
imp$keep[is.na(imp$keep)] <- FALSE

# Short, readable labels for variables (avoid overlapping unreadable names)
short_label <- function(x, max_len = 32) {
  out <- x
  for (i in seq_along(x)) {
    if (nchar(x[i]) > max_len) {
      # Prefer last meaningful part (e.g. after last dot) or truncate
      if (grepl("\\.", x[i])) {
        parts <- strsplit(x[i], "\\.")[[1]]
        out[i] <- parts[length(parts)]
      }
      if (nchar(out[i]) > max_len)
        out[i] <- paste0(substr(out[i], 1L, max_len - 3L), "...")
    }
  }
  out
}
imp$variable_label <- short_label(imp$variable)

# Order variables by mean importance (across models); use variable_label for display
var_order <- imp %>%
  group_by(variable, variable_label) %>%
  summarise(mean_inc = mean(rmse_increase, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_inc))
imp$variable_label <- factor(imp$variable_label, levels = rev(unique(var_order$variable_label)))

# Optional: show only top N variables per model to reduce clutter (set to Inf to show all)
top_n_per_model <- 14
imp_plot <- imp %>%
  group_by(model) %>%
  arrange(desc(rmse_increase), .by_group = TRUE) %>%
  slice_head(n = top_n_per_model) %>%
  ungroup()
# Order variable_label by mean importance (so most important at top of bar chart)
lvl <- imp_plot %>%
  group_by(variable_label) %>%
  summarise(m = mean(rmse_increase, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(m)) %>%
  pull(variable_label) %>%
  rev()
imp_plot$variable_label <- factor(imp_plot$variable_label, levels = lvl)

# 1) Permutation importance by model and variable; colour by keep.
#    sqrt scale so the full ranking is visible (not dominated by top 2-3 bars).
p_importance <- ggplot(imp_plot, aes(x = variable_label, y = rmse_increase, fill = keep)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~ model, scales = "free_x", ncol = 3) +
  scale_y_sqrt() +
  scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "grey70"), name = "Kept") +
  labs(
    x = "Predictor",
    y = "RMSE increase (permutation importance; sqrt scale)",
    title = if (show_titles) "Model-wise permutation importance (1 km spatial block CV)" else NULL,
    subtitle = if (show_titles) paste0("Top ", top_n_per_model, " predictors per model; blue = retained. Sqrt scale shows full ranking.") else NULL
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 11),
    plot.subtitle = element_text(hjust = 0.5, size = 8)
  )

# 2) Number of variables kept per model (summary bar)
n_kept <- imp %>%
  group_by(model) %>%
  summarise(n_kept = sum(keep, na.rm = TRUE), .groups = "drop")
p_n_kept <- ggplot(n_kept, aes(x = reorder(model, n_kept), y = n_kept, fill = model)) +
  geom_col(show.legend = FALSE) +
  coord_flip() +
  labs(x = NULL, y = "Number of predictors kept", title = if (show_titles) "Predictors kept per model (cumulative importance)" else NULL) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Save
if (requireNamespace("patchwork", quietly = TRUE)) {
  p_combined <- p_n_kept + p_importance + patchwork::plot_layout(ncol = 1, heights = c(0.3, 0.7))
  ggsave(file.path(out_dir, "model_permutation_pruning_importance.png"), p_combined, width = 12, height = 2 + 0.22 * length(unique(imp_plot$variable_label)), dpi = dpi, limitsize = FALSE)
} else {
  ggsave(file.path(out_dir, "model_permutation_pruning_importance.png"), p_importance, width = 12, height = 2 + 0.22 * length(unique(imp_plot$variable_label)), dpi = dpi, limitsize = FALSE)
}
cat("Saved", file.path(out_dir, "model_permutation_pruning_importance.png"), "\n")

ggsave(file.path(out_dir, "model_permutation_pruning_n_kept.png"), p_n_kept, width = 6, height = 3, dpi = dpi)
cat("Saved", file.path(out_dir, "model_permutation_pruning_n_kept.png"), "\n")

# ---------------------------------------------------------------------------
# GPR-only importance (variable SELECTION step): importance from spatial CV
# used to choose the 14 predictors. Not the same as "best fitted model" importance
# (that is in gpr_feature_importance_best_model.png, also sqrt scale).
# ---------------------------------------------------------------------------
source("modelling/R/plot_config.R")
imp_gpr <- imp %>% filter(model == "GPR")
if (nrow(imp_gpr) > 0L) {
  imp_gpr$variable_label <- label_vars(imp_gpr$variable)
  imp_gpr <- imp_gpr %>% arrange(desc(rmse_increase))
  imp_gpr$variable_label <- factor(imp_gpr$variable_label, levels = rev(unique(imp_gpr$variable_label)))
  p_gpr <- ggplot(imp_gpr, aes(x = variable_label, y = rmse_increase, fill = keep)) +
    geom_col() +
    coord_flip() +
    scale_y_sqrt() +
    scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "grey70"), name = "Retained") +
    labs(x = NULL, y = "RMSE increase (permutation importance; sqrt scale)") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom", axis.text.y = element_text(size = rel(0.9)))
  ggsave(file.path(out_dir, "model_permutation_pruning_importance_gpr_paper.png"), p_gpr,
         width = 8, height = max(5, 0.32 * nrow(imp_gpr)), dpi = dpi)
  cat("Saved", file.path(out_dir, "model_permutation_pruning_importance_gpr_paper.png"), " (GPR-only, for paper)\n")
}

cat("Done.\n")
