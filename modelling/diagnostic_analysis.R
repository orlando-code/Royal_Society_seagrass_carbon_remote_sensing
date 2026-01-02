### Diagnostic Analysis for Seagrass Carbon Stock Model
### This script identifies potential issues with the modeling approach

library(tidyverse)
library(mgcv)

# Load data
setwd("/Users/rt582/Library/CloudStorage/OneDrive-UniversityofCambridge/cambridge/phd/Paper_Conferences/seagrass/Fw_ Following up on modeling conversation at the Royal Society Meeting/SG_Cstock_modeling")
data <- read_rds("data/For_modeling_df_shallow_carbon_density.rds")

cat("=== DIAGNOSTIC ANALYSIS ===\n\n")

# 1. Check pseudo-replication
cat("1. PSEUDO-REPLICATION CHECK\n")
cat("   Total samples:", nrow(data), "\n")
cat("   Unique cores:", length(unique(data$random_core_variable)), "\n")
cat("   Unique locations (lat/lon):", nrow(unique(data[, c("latitude", "longitude")])), "\n")
cat("   Mean samples per core:", mean(table(data$random_core_variable)), "\n")
cat("   Max samples per core:", max(table(data$random_core_variable)), "\n\n")

# 2. Check remote sensing resolution
cat("2. REMOTE SENSING RESOLUTION CHECK\n")
rs_vars <- c("KD_closest", "RRS443_closest", "po4_mean_1.5m_mmol_m3_closest", 
             "wave_height_VHM0_p95_m_closest", "pH_mean_1.5m_closest")
for(var in rs_vars) {
  n_unique <- length(unique(data[[var]]))
  cat(sprintf("   %s: %d unique values for %d samples (%.1f%% unique)\n", 
              var, n_unique, nrow(data), 100*n_unique/nrow(data)))
}

# Check combinations
rs_combos <- data %>%
  select(all_of(rs_vars)) %>%
  distinct()
cat("\n   Unique combinations of key RS variables:", nrow(rs_combos), "\n")

shared_rs <- data %>%
  group_by(across(all_of(rs_vars))) %>%
  summarise(n_samples = n(), 
            n_cores = n_distinct(random_core_variable),
            .groups = "drop")
cat("   Max samples per RS combination:", max(shared_rs$n_samples), "\n")
cat("   Max cores per RS combination:", max(shared_rs$n_cores), "\n\n")

# 2b. Check environmental variable constancy within cores (PSEUDOREPLICATION TEST)
cat("2b. ENVIRONMENTAL VARIABLE CONSTANCY WITHIN CORES (PSEUDOREPLICATION TEST)\n")
cat("   This tests whether environmental variables are constant within cores.\n")
cat("   If they are, samples within cores are NOT independent observations.\n\n")

# Get all environmental predictor variables from the model
env_vars <- c("KD_closest", "RRS443_closest", "wave_height_VHM0_p95_m_closest",
              "po4_mean_1.5m_mmol_m3_closest", "pH_mean_1.5m_closest",
              "bottomT_p95_C_closest", "vo_p90_1.5m_m_s_closest",
              "uo_mean_1.5m_m_s_closest", "Surf_fgco2_p95_molC_m2_yr_closest")

# Check for each variable: how many cores have constant values across samples?
within_core_constancy <- data %>%
  group_by(random_core_variable) %>%
  summarise(
    across(all_of(env_vars), 
           function(x) length(unique(x[!is.na(x)])) == 1,  # TRUE if all values are the same
           .names = "{.col}_constant"),
    n_samples = n(),
    .groups = "drop"
  )

# Calculate proportion of cores with constant values for each variable
constancy_summary <- within_core_constancy %>%
  summarise(
    across(ends_with("_constant"), 
           function(x) sum(x, na.rm = TRUE) / length(x),
           .names = "{.col}_prop")
  )

cat("   Proportion of cores where environmental variables are CONSTANT across all samples:\n")
for(var in env_vars) {
  prop_constant <- constancy_summary[[paste0(var, "_constant_prop")]]
  n_cores_constant <- sum(within_core_constancy[[paste0(var, "_constant")]], na.rm = TRUE)
  cat(sprintf("   %-35s: %.1f%% (%d/%d cores)\n", 
              var, prop_constant * 100, n_cores_constant, nrow(within_core_constancy)))
}

# Calculate within-core vs between-core variance
cat("\n   Variance decomposition (within-core vs between-core):\n")
variance_decomp <- map_dfr(env_vars, function(var) {
  # Between-core variance (variance of core means)
  core_means <- data %>%
    group_by(random_core_variable) %>%
    summarise(mean_val = mean(.data[[var]], na.rm = TRUE), .groups = "drop")
  between_var <- var(core_means$mean_val, na.rm = TRUE)
  
  # Within-core variance (average variance within each core)
  within_vars <- data %>%
    group_by(random_core_variable) %>%
    summarise(within_var = var(.data[[var]], na.rm = TRUE), .groups = "drop")
  within_var <- mean(within_vars$within_var, na.rm = TRUE)
  
  # Total variance
  total_var <- var(data[[var]], na.rm = TRUE)
  
  # Proportion of variance that is within-core
  prop_within <- ifelse(total_var > 0, within_var / total_var, 0)
  
  data.frame(
    variable = var,
    between_core_var = between_var,
    within_core_var = within_var,
    total_var = total_var,
    prop_variance_within = prop_within
  )
})

for(i in seq_len(nrow(variance_decomp))) {
  vd <- variance_decomp[i, ]
  cat(sprintf("   %-35s: %.1f%% of variance is within-core\n", 
              vd$variable, vd$prop_variance_within * 100))
}

# Count cores with ALL environmental variables constant
cores_all_constant <- within_core_constancy %>%
  filter(if_all(ends_with("_constant"), ~ .x == TRUE)) %>%
  nrow()
cat(sprintf("\n   Cores with ALL environmental variables constant: %d/%d (%.1f%%)\n",
            cores_all_constant, nrow(within_core_constancy),
            100 * cores_all_constant / nrow(within_core_constancy)))

cat("\n   ⚠️  WHY THIS PSEUDOREPLICATION IS PROBLEMATIC:\n")
cat("   1. If environmental variables are constant within cores, samples from the same\n")
cat("      core are NOT independent - they share identical predictor values.\n")
cat("   2. The model treats these as separate observations, artificially inflating the\n")
cat("      effective sample size from 382 (cores) to 3,680 (samples).\n")
cat("   3. This violates the independence assumption of regression models.\n")
cat("   4. Standard errors and p-values will be underestimated (false precision).\n")
cat("   5. R² values are inflated because the model can 'memorize' core-specific patterns\n")
cat("      through the random effect, while environmental variables appear predictive\n")
cat("      simply because they're constant within cores.\n")
cat("   6. For spatial prediction, you need to predict NEW cores, not interpolate\n")
cat("      within existing cores. The current approach doesn't test this capability.\n\n")

# 3. Check core-level correlations with RS variables
cat("3. CORE-LEVEL CORRELATIONS WITH RS VARIABLES\n")
core_means <- data %>%
  group_by(random_core_variable) %>%
  summarise(
    mean_carbon = mean(carbon_density_g_c_cm3),
    KD = first(KD_closest),
    RRS443 = first(RRS443_closest),
    po4 = first(po4_mean_1.5m_mmol_m3_closest),
    .groups = "drop"
  )

cat("   Correlation (core mean carbon ~ KD):", 
    round(cor(core_means$mean_carbon, core_means$KD, use = "complete.obs"), 3), "\n")
cat("   Correlation (core mean carbon ~ RRS443):", 
    round(cor(core_means$mean_carbon, core_means$RRS443, use = "complete.obs"), 3), "\n")
cat("   Correlation (core mean carbon ~ po4):", 
    round(cor(core_means$mean_carbon, core_means$po4, use = "complete.obs"), 3), "\n\n")

# 4. Assess random effect contribution
cat("4. RANDOM EFFECT CONTRIBUTION\n")
cat("   Fitting model with ONLY random effect...\n")
model_re_only <- gam(carbon_density_g_c_cm3 ~ 
                     s(random_core_variable, bs = "re"),
                     family = Gamma(link = "log"),
                     data = data, method = "REML")
summary_re <- summary(model_re_only)
cat("   R² (adj) with only random effect:", round(summary_re$r.sq, 3), "\n")
cat("   Deviance explained:", round(summary_re$dev.expl * 100, 1), "%\n\n")

# 5. Compare to full model
cat("5. COMPARISON WITH ENVIRONMENTAL VARIABLES\n")
cat("   Fitting reduced model with environmental variables...\n")
# model_full <- gam(carbon_density_g_c_cm3 ~ 
#                   s(KD_closest) + s(RRS443_closest) + s(wave_height_VHM0_p95_m_closest) +
#                   s(po4_mean_1.5m_mmol_m3_closest) + s(pH_mean_1.5m_closest) + 
#                   s(bottomT_p95_C_closest) + s(vo_p90_1.5m_m_s_closest) + 
#                   s(uo_mean_1.5m_m_s_closest) + s(Surf_fgco2_p95_molC_m2_yr_closest) +
#                   s(sediment_mean_depth_cm) + seagrass_species + 
#                   s(random_core_variable, bs = "re"),
#                   family = Gamma(link = "log"),
#                   data = data, method = "REML")
# summary_full <- summary(model_full)
# cat("   R² (adj) with environmental variables:", round(summary_full$r.sq, 3), "\n")
# cat("   Deviance explained:", round(summary_full$dev.expl * 100, 1), "%\n")
# cat("   Additional R² from environmental variables:", 
#     round(summary_full$r.sq - summary_re$r.sq, 3), "\n\n")

# 6. Check spatial clustering
cat("6. SPATIAL CLUSTERING\n")
spatial_summary <- data %>%
  group_by(latitude, longitude) %>%
  summarise(
    n_cores = n_distinct(random_core_variable),
    n_samples = n(),
    .groups = "drop"
  )
cat("   Locations with multiple cores:", sum(spatial_summary$n_cores > 1), "\n")
cat("   Max cores per location:", max(spatial_summary$n_cores), "\n")
cat("   Locations with >10 cores:", sum(spatial_summary$n_cores > 10), "\n\n")

# 7. Core-level model performance
cat("7. CORE-LEVEL MODEL PERFORMANCE\n")
# Predict using full model
data$predicted <- predict(model_full, newdata = data, type = "response")

# Aggregate to core level
core_performance <- data %>%
  group_by(random_core_variable) %>%
  summarise(
    obs_mean = mean(carbon_density_g_c_cm3),
    pred_mean = mean(predicted),
    n_samples = n(),
    .groups = "drop"
  )

core_r2 <- cor(core_performance$obs_mean, core_performance$pred_mean)^2
core_rmse <- sqrt(mean((core_performance$obs_mean - core_performance$pred_mean)^2))

cat("   Core-level R²:", round(core_r2, 3), "\n")
cat("   Core-level RMSE:", round(core_rmse, 4), "\n")
cat("   (Compare to sample-level R²:", round(summary_full$r.sq, 3), ")\n\n")

# 8. Check for cross-validation
cat("8. CROSS-VALIDATION STATUS\n")
cat("   ⚠️  WARNING: No cross-validation found in current code\n")
cat("   Recommendation: Implement core-level cross-validation\n\n")

cat("=== END OF DIAGNOSTIC ANALYSIS ===\n")

# Save diagnostic results
diagnostic_results <- list(
  n_samples = nrow(data),
  n_cores = length(unique(data$random_core_variable)),
  n_locations = nrow(unique(data[, c("latitude", "longitude")])),
  rs_resolution = sapply(rs_vars, function(v) length(unique(data[[v]]))),
  env_var_constancy = constancy_summary,
  variance_decomposition = variance_decomp,
  cores_all_env_constant = cores_all_constant,
  core_level_correlations = c(
    KD = cor(core_means$mean_carbon, core_means$KD, use = "complete.obs"),
    RRS443 = cor(core_means$mean_carbon, core_means$RRS443, use = "complete.obs"),
    po4 = cor(core_means$mean_carbon, core_means$po4, use = "complete.obs")
  ),
  r2_random_effect_only = summary_re$r.sq,
  r2_full_model = summary_full$r.sq,
  r2_increment = summary_full$r.sq - summary_re$r.sq,
  core_level_r2 = core_r2,
  core_level_rmse = core_rmse
)

saveRDS(diagnostic_results, "diagnostic_results.rds")
cat("Diagnostic results saved to diagnostic_results.rds\n")

