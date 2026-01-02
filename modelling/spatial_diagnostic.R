# ============================================
# DIAGNOSTIC: Why do environmental variables work alone but add little with random effect?
# ============================================

cat("\n=== DIAGNOSTIC: Environmental vs Random Effect Confounding ===\n\n")


# 1. Can environmental variables predict which core a sample comes from?
cat("1. CAN ENVIRONMENTAL VARIABLES PREDICT CORE IDENTITY?\n")
cat("   If yes, they're redundant with the random effect.\n\n")

# Create a dataset with core-level environmental variables
core_env <- data %>%
  group_by(random_core_variable) %>%
  summarise(
    KD = first(KD_closest),
    RRS443 = first(RRS443_closest),
    po4 = first(po4_mean_1.5m_mmol_m3_closest),
    pH = first(pH_mean_1.5m_closest),
    wave = first(wave_height_VHM0_p95_m_closest),
    bottomT = first(bottomT_p95_C_closest),
    vo = first(vo_p90_1.5m_m_s_closest),
    uo = first(uo_mean_1.5m_m_s_closest),
    fgco2 = first(Surf_fgco2_p95_molC_m2_yr_closest),
    species = first(seagrass_species),
    lat = first(latitude),
    lon = first(longitude),
    .groups = "drop"
  )

# How many unique environmental "fingerprints" are there?
env_fingerprint <- core_env %>%
  select(KD, RRS443, po4, pH, wave, bottomT, vo, uo, fgco2, species) %>%
  distinct()

cat("   Unique environmental 'fingerprints' (combinations of env vars + species):", 
    nrow(env_fingerprint), "\n")
cat("   Number of cores:", nrow(core_env), "\n")
cat("   Average cores per unique fingerprint:", round(nrow(core_env) / nrow(env_fingerprint), 2), "\n")
cat("   If this is close to 1, environmental variables uniquely identify cores.\n\n")

# 2. Check spatial structure of environmental variables
cat("2. SPATIAL STRUCTURE OF ENVIRONMENTAL VARIABLES\n")
cat("   Are environmental variables spatially autocorrelated?\n\n")

# Calculate spatial distances between cores
# install.packages("geosphere")
library(geosphere)
core_coords <- core_env %>%
  select(random_core_variable, lat, lon) %>%
  distinct()

# Sample a subset for distance calculation (can be slow with many cores)
if(nrow(core_coords) <= 500) {
  dist_matrix <- distm(core_coords[, c("lon", "lat")], 
                       core_coords[, c("lon", "lat")]) / 1000  # km
  
  # For each environmental variable, check spatial autocorrelation
  env_vars_spatial <- c("KD", "RRS443", "po4", "pH", "wave", "bottomT", "vo", "uo", "fgco2")
  
  for(var in env_vars_spatial) {
    if(var %in% names(core_env)) {
      # Calculate correlation between environmental variable and spatial distance
      # (simplified: check if nearby cores have similar values)
      var_vals <- core_env[[var]]
      # Get mean absolute difference for cores within 10km vs >100km
      nearby_diffs <- c()
      far_diffs <- c()
      
      for(i in 1:min(100, nrow(core_coords))) {  # Sample for speed
        for(j in (i+1):min(100, nrow(core_coords))) {
          dist_km <- dist_matrix[i, j]
          val_diff <- abs(var_vals[i] - var_vals[j])
          
          if(dist_km < 10) {
            nearby_diffs <- c(nearby_diffs, val_diff)
          } else if(dist_km > 100) {
            far_diffs <- c(far_diffs, val_diff)
          }
        }
      }
      
      if(length(nearby_diffs) > 0 && length(far_diffs) > 0) {
        mean_nearby <- mean(nearby_diffs, na.rm = TRUE)
        mean_far <- mean(far_diffs, na.rm = TRUE)
        cat(sprintf("   %-20s: Nearby cores diff = %.4f, Far cores diff = %.4f\n",
                    var, mean_nearby, mean_far))
      }
    }
  }
} else {
  cat("   Too many cores for full distance calculation. Skipping.\n")
}

cat("\n")

# 3. Compare predictions: environmental vs random effect
cat("3. PREDICTION COMPARISON\n")
cat("   How well do environmental variables predict core-level means?\n\n")

# Fit environmental-only model
model_env_only <- gam(carbon_density_g_c_cm3 ~ 
                      s(KD_closest) + s(RRS443_closest) + s(wave_height_VHM0_p95_m_closest) +
                      s(po4_mean_1.5m_mmol_m3_closest) + s(pH_mean_1.5m_closest) + 
                      s(bottomT_p95_C_closest) + s(vo_p90_1.5m_m_s_closest) + 
                      s(uo_mean_1.5m_m_s_closest) + s(Surf_fgco2_p95_molC_m2_yr_closest) +
                      s(sediment_mean_depth_cm) + seagrass_species,
                      family = Gamma(link = "log"),
                      data = data, method = "REML")

# Fit random effect only model
model_re_only <- gam(carbon_density_g_c_cm3 ~ 
                     s(random_core_variable, bs = "re"),
                     family = Gamma(link = "log"),
                     data = data, method = "REML")

# Predict at core level
# Predict at core level -- safer and bug-free: compute predictions for all data, then aggregate
data$pred_env <- predict(model_env_only, newdata = data, type = "response")
data$pred_re  <- predict(model_re_only,  newdata = data, type = "response")

core_predictions <- data %>%
  group_by(random_core_variable) %>%
  summarise(
    obs_mean  = mean(carbon_density_g_c_cm3),
    pred_env  = mean(pred_env),
    pred_re   = mean(pred_re),
    .groups = "drop"
  )

# Core-level R²
r2_env_cores <- cor(core_predictions$obs_mean, core_predictions$pred_env)^2
r2_re_cores <- cor(core_predictions$obs_mean, core_predictions$pred_re)^2

cat("   Core-level R² (environmental variables only):", round(r2_env_cores, 3), "\n")
cat("   Core-level R² (random effect only):", round(r2_re_cores, 3), "\n")
cat("   Difference:", round(r2_re_cores - r2_env_cores, 3), "\n\n")

# 4. Check if environmental variables vary between cores that share same location
cat("4. ENVIRONMENTAL VARIABLES AT SAME LOCATION\n")
cat("   Do cores at the same lat/lon have different environmental values?\n\n")

same_location_cores <- data %>%
  group_by(latitude, longitude) %>%
  filter(n_distinct(random_core_variable) > 1) %>%
  group_by(latitude, longitude, random_core_variable) %>%
  summarise(
    KD = first(KD_closest),
    RRS443 = first(RRS443_closest),
    .groups = "drop"
  ) %>%
  group_by(latitude, longitude) %>%
  summarise(
    n_cores = n(),
    KD_unique = n_distinct(KD),
    RRS443_unique = n_distinct(RRS443),
    .groups = "drop"
  )

if(nrow(same_location_cores) > 0) {
  cat("   Locations with multiple cores:", nrow(same_location_cores), "\n")
  cat("   Of these, locations where cores have different KD values:", 
      sum(same_location_cores$KD_unique > 1), "\n")
  cat("   Of these, locations where cores have different RRS443 values:", 
      sum(same_location_cores$RRS443_unique > 1), "\n")
  cat("   If cores at same location have same env values, they're spatially confounded.\n\n")
} else {
  cat("   No locations with multiple cores found.\n\n")
}

# 5. Correlation between random effect estimates and environmental variables
cat("5. CORRELATION: RANDOM EFFECT ESTIMATES vs ENVIRONMENTAL VARIABLES\n")
cat("   If high correlation, they're capturing the same information.\n\n")

# Extract random effect coefficients
re_coefs <- coef(model_re_only)
re_names <- names(re_coefs)[grepl("random_core_variable", names(re_coefs))]

# Match to core environmental variables
re_df <- data.frame(
  random_core_variable = gsub("s\\(random_core_variable\\)\\.", "", re_names),
  re_estimate = re_coefs[re_names]
)

core_env_re <- core_env %>%
  left_join(re_df, by = "random_core_variable")

# Correlations
cat("   Correlation (RE estimate ~ KD):", 
    round(cor(core_env_re$re_estimate, core_env_re$KD, use = "complete.obs"), 3), "\n")
cat("   Correlation (RE estimate ~ RRS443):", 
    round(cor(core_env_re$re_estimate, core_env_re$RRS443, use = "complete.obs"), 3), "\n")
cat("   Correlation (RE estimate ~ po4):", 
    round(cor(core_env_re$re_estimate, core_env_re$po4, use = "complete.obs"), 3), "\n")
cat("   If correlations are high (>0.5), RE and env vars are confounded.\n\n")

cat("=== END DIAGNOSTIC ===\n\n")