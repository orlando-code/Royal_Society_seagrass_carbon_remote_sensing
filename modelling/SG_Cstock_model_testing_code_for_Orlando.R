### Segarass carbon stock model for testing
### Code prepared by Natalya Gallo to share with Orlando Timmerman
### December 17, 2025

### Script uses outputs of previous scripts: SG_carbonstock_model.R (where multiple model formulations were tested) 
### and dataset produced in script Environmental_covariates_matching.R

# Refresh working directory
rm(list = ls())

# install.packages("visreg")
# install.packages("randomForest")


## Open needed packages
library(here)
library(mgcv)
library(tidyverse)
library(ggplot2)
library(visreg) 
library(randomForest)
library(corrplot)

# Set working directory
# setwd(here::here())
setwd("/Users/rt582/Library/CloudStorage/OneDrive-UniversityofCambridge/cambridge/phd/Paper_Conferences/seagrass/Fw_ Following up on modeling conversation at the Royal Society Meeting/SG_Cstock_modeling")


# Open dataset prepared for modeling (this is the same dataset used by Katalin Blix for machine learning models and
# myself for non-machine learning models)
For_modeling_df_shallow_carbon_density <- read_rds("data/For_modeling_df_shallow_carbon_density.rds")
summary(For_modeling_df_shallow_carbon_density)


# export rds to csv for investigation in Python
write_csv(For_modeling_df_shallow_carbon_density, "data/For_modeling_df_shallow_carbon_density.csv")

# Note: the random core variable was created by pasting the location_name, station_id, core_id, year, month, day,
# latitude, longitude, and water_depth_m (separated by "_"), so should represent each unique core in the dataset 
# (since several samples can come from different depths of each core)
str(For_modeling_df_shallow_carbon_density$random_core_variable) #382 levels
# Can have a closer look at this: summarize the number of samples per core 
summary_cores <- For_modeling_df_shallow_carbon_density %>%
  group_by(random_core_variable) %>%
  summarise(sample_count = n())
summary(summary_cores) #mean: 9.6, median: 4, max: 86
hist(summary_cores$sample_count) #I am worried about overfitting related to having a high sample number (3680) but a 
# low number of actual unique cores (382)



# ================================ Experimenting
carbon_density_re_only <- gam(carbon_density_g_c_cm3 ~ s(random_core_variable, bs = "re"), data = For_modeling_df_shallow_carbon_density, method = "REML")
carbon_density_re_only_log <- gam(carbon_density_g_c_cm3 ~ s(random_core_variable, bs = "re"), family = Gamma(link = "log"), data = For_modeling_df_shallow_carbon_density, method = "REML")

print("Summary of random core only model with identity link function:")
summary(carbon_density_re_only) # R-sq.(adj) = 0.906, Deviance explained = 91.5%
print("Summary of random core only model with log link function:")
summary(carbon_density_re_only_log) # R-sq.(adj) = 0.883, Deviance explained = 76.7%
# ie. incorporating random core term alone does better than when supplemented with remote sensing variables


# what r2 can you get without random core term?
# ================================

# This was the top reduced model we arrived at as a final recommended model, and that I presented the main results from
# during the Royal Society meeting
# The model can be open directly, or rerun:
GAM_top_reduced_SGstock <- readRDS("GAM_top_reduced_model_SGstock.rds")

GAM_top_reduced_SGstock_no_random <- gam(carbon_density_g_c_cm3 ~ 
                                 s(KD_closest) + s(RRS443_closest) + s(wave_height_VHM0_p95_m_closest) +
                                 s(po4_mean_1.5m_mmol_m3_closest) + s(pH_mean_1.5m_closest) + 
                                 s(bottomT_p95_C_closest) + s(vo_p90_1.5m_m_s_closest) + 
                                 s(uo_mean_1.5m_m_s_closest) + s(Surf_fgco2_p95_molC_m2_yr_closest) +
                                 s(sediment_mean_depth_cm) + seagrass_species,
                                  family = Gamma(link = "log"),
                               data = For_modeling_df_shallow_carbon_density, method = "REML")

summary(GAM_top_reduced_SGstock_no_random)

GAM_top_reduced_SGstock <- gam(carbon_density_g_c_cm3 ~ 
                                 s(KD_closest) + s(RRS443_closest) + s(wave_height_VHM0_p95_m_closest) +
                                 s(po4_mean_1.5m_mmol_m3_closest) + s(pH_mean_1.5m_closest) + 
                                 s(bottomT_p95_C_closest) + s(vo_p90_1.5m_m_s_closest) + 
                                 s(uo_mean_1.5m_m_s_closest) + s(Surf_fgco2_p95_molC_m2_yr_closest) +
                                 s(sediment_mean_depth_cm) + seagrass_species + s(random_core_variable, bs = "re"),
                                  family = Gamma(link = "log"),
                               data = For_modeling_df_shallow_carbon_density, method = "REML")

summary(GAM_top_reduced_SGstock) # Model summary
#R-sq.(adj) = 0.905  Deviance explained = 76.7%
visreg(GAM_top_reduced_SGstock) #Note some high outlier points for phosphate and pH. I think I should drop these and rerun. # TODO: how does this show this?

# Full model structure tested for all modeling approaches in first iteration:
# carbon_density_g_c_cm3 ~ all_env_covariates_closest + sediment_mean_depth_cm + random_core_variable + seagrass_species + Region 
# random_core_variable, Region, and seagrass_species are factors
# "all_env_covariates_closest" represents 42 environmental predictor variables tested
# sediment_mean_depth_cm is the decompacted middle depth for each sample (midway between start and end depth)

# Full generalized additive model with all 42 environmental predictor variables
# (Note: lots of correlation between variables, so I tested that the reduced model did not have issues with 
# correlation with the retained variables)

gam_model_full <- bam(carbon_density_g_c_cm3 ~ s(KD_closest) + s(PBS443_closest) + s(Chla_closest) +
                                   s(RRS620_closest) + s(CDOM_closest) + s(RRS443_closest) +
                                   s(bottomT_mean_C_closest) + s(bottomT_p95_C_closest) + s(sal_mean_0.5m_closest) +
                                   s(sal_p05_0.5m_closest) + s(SST_daily_mean_C_closest) + s(SST_daily_p95_C_closest) +
                                   s(Surf_talk_mean_umol_kg_closest) + s(Surf_talk_p95_umol_kg_closest) + s(Surf_spco2_mean_uatm_closest) +
                                   s(Surf_spco2_p95_uatm_closest) + s(Surf_fgco2_mean_molC_m2_yr_closest) + s(Surf_fgco2_p95_molC_m2_yr_closest) +
                                   s(fe_mean_1.5m_mmol_m3_closest) + s(no3_mean_1.5m_mmol_m3_closest) + s(po4_mean_1.5m_mmol_m3_closest) +
                                   s(phyc_mean_0.5m_mmol_m3_closest) + s(nppv_mean_0.5m_mg_m3_day_closest) + s(pH_mean_1.5m_closest) +
                                   s(pH_p05_1.5m_closest) + s(spco2_monthly_mean_Pa_closest) + s(spco2_monthly_p95_Pa_closest) +
                                   s(ox_mean_2.7m_mmol_m3_closest) + s(ox_p05_2.7m_mmol_m3_closest) + s(depth_0.083deg_m_closest) +
                                   s(uo_mean_1.5m_m_s_closest) + s(uo_p90_1.5m_m_s_closest) + s(vo_mean_1.5m_m_s_closest) +
                                   s(vo_p90_1.5m_m_s_closest) + s(wave_height_VHM0_mean_m_closest) + s(wave_height_VHM0_p95_m_closest) +
                                   s(SW_KE_mean_cm2_s2_closest) + s(SW_KE_p95_cm2_s2_closest) + s(Daily_max_wave_height_mean_m_closest) +
                                   s(Daily_mean_wave_height_mean_m_closest) + s(Daily_std_wave_height_mean_m_closest) + s(No_sig_wave_heights_mean_closest) +
                                   s(sediment_mean_depth_cm) +
                                   s(random_core_variable, bs = "re") + seagrass_species + Region,
                                 family = Gamma(link = "log"), data = For_modeling_df_shallow_carbon_density, 
                                 method = "fREML", select = TRUE)
summary(gam_model_full_no_species) #R-sq.(adj) =  0.801   Deviance explained = 81.9%  # TODO: this is not present

# Random forest model using all predictors
# Note, random_core_variable not included in random forest modeling because Random Forests don't support random effects
# directly, so it could only be modeled as a regular categorical variable, which would then make the model not applicable
# for applying to a new dataset with different cores. 
rf_carbon_density_shallow <- randomForest(carbon_density_g_c_cm3 ~ KD_closest + PBS443_closest + Chla_closest +
                                            RRS620_closest + CDOM_closest + RRS443_closest + bottomT_mean_C_closest +
                                            bottomT_p95_C_closest + sal_mean_0.5m_closest + sal_p05_0.5m_closest +
                                            SST_daily_mean_C_closest + SST_daily_p95_C_closest + Surf_talk_mean_umol_kg_closest +
                                            Surf_talk_p95_umol_kg_closest + Surf_spco2_mean_uatm_closest + Surf_spco2_p95_uatm_closest +
                                            Surf_fgco2_mean_molC_m2_yr_closest + Surf_fgco2_p95_molC_m2_yr_closest + fe_mean_1.5m_mmol_m3_closest +
                                            no3_mean_1.5m_mmol_m3_closest + po4_mean_1.5m_mmol_m3_closest + phyc_mean_0.5m_mmol_m3_closest +
                                            nppv_mean_0.5m_mg_m3_day_closest + pH_mean_1.5m_closest + pH_p05_1.5m_closest +
                                            spco2_monthly_mean_Pa_closest + spco2_monthly_p95_Pa_closest + ox_mean_2.7m_mmol_m3_closest +
                                            ox_p05_2.7m_mmol_m3_closest + depth_0.083deg_m_closest + uo_mean_1.5m_m_s_closest +
                                            uo_p90_1.5m_m_s_closest + vo_mean_1.5m_m_s_closest + vo_p90_1.5m_m_s_closest +
                                            wave_height_VHM0_mean_m_closest + wave_height_VHM0_p95_m_closest + SW_KE_mean_cm2_s2_closest +
                                            SW_KE_p95_cm2_s2_closest + Daily_max_wave_height_mean_m_closest + Daily_mean_wave_height_mean_m_closest +
                                            Daily_std_wave_height_mean_m_closest + No_sig_wave_heights_mean_closest +
                                            seagrass_species + Region + sediment_mean_depth_cm, data = For_modeling_df_shallow_carbon_density)
print(rf_carbon_density_shallow) 

# Have a look at correlation across environmental predictors
Env_corr <- For_modeling_df_shallow_carbon_density %>% select(-c(number_id_final_version, latitude,
                                                                 longitude, dry_bulk_density_g_cm3, 
                                                                 organic_carbon_percent, carbon_density_g_c_cm3,
                                                                 seagrass_species, sediment_mean_depth_cm, Region,
                                                                 random_core_variable))
cor_matrix <- cor(Env_corr, use = "complete.obs", method = "pearson")
print(cor_matrix)
corrplot(cor_matrix, 
         method = "color",
         type = "upper",
         order = "hclust",
         col = colorRampPalette(c("blue", "white", "red"))(200), # Custom color palette
         #addCoef.col = "black", # Add correlation coefficients
         tl.col = "black",      # Text label color
         tl.srt = 90,           # Text label rotation
         diag = FALSE           # Remove diagonal elements
)

# Get the upper triangle of the correlation matrix
upper_tri <- cor_matrix
upper_tri[lower.tri(upper_tri, diag = TRUE)] <- NA
# Find variable pairs with correlation > 0.8
high_corr_pairs <- which(abs(upper_tri) > 0.8, arr.ind = TRUE)
# Create a data frame of the results
high_corr_df <- data.frame(
  Var1 = rownames(upper_tri)[high_corr_pairs[, 1]],
  Var2 = colnames(upper_tri)[high_corr_pairs[, 2]],
  Correlation = upper_tri[high_corr_pairs]
)
#Order the data frame by absolute correlation value, descending
high_corr_df <- high_corr_df[order(-abs(high_corr_df$Correlation)), ]
# View the sorted result
print(high_corr_df)

#Remove all but the 9 key predictor variables
#Carbon density (gC cm⁻³) ~ KD490 + RRS443 + VHM0_p95 + Phosphate_mean + pH_mean +
#Bottom_T_p95 + Vo_p90 + Uo_mean + fgCO₂_p95 + Sediment depth + Seagrass species
Env_corr_reduced <- Env_corr %>% select(-c(RRS620_closest, CDOM_closest, SW_KE_mean_cm2_s2_closest,
                                           SW_KE_p95_cm2_s2_closest, Surf_talk_p95_umol_kg_closest, Surf_talk_mean_umol_kg_closest,
                                           ox_p05_2.7m_mmol_m3_closest, ox_mean_2.7m_mmol_m3_closest, Chla_closest,
                                           sal_mean_0.5m_closest, sal_p05_0.5m_closest, SST_daily_mean_C_closest,
                                           SST_daily_p95_C_closest, wave_height_VHM0_mean_m_closest, Daily_mean_wave_height_mean_m_closest,
                                           Daily_std_wave_height_mean_m_closest, vo_mean_1.5m_m_s_closest, 
                                           uo_p90_1.5m_m_s_closest, Surf_fgco2_mean_molC_m2_yr_closest, Daily_max_wave_height_mean_m_closest,
                                           No_sig_wave_heights_mean_closest,  depth_0.083deg_m_closest, spco2_monthly_p95_Pa_closest,
                                           spco2_monthly_mean_Pa_closest, pH_p05_1.5m_closest, nppv_mean_0.5m_mg_m3_day_closest,
                                           phyc_mean_0.5m_mmol_m3_closest, no3_mean_1.5m_mmol_m3_closest, fe_mean_1.5m_mmol_m3_closest,
                                           Surf_spco2_p95_uatm_closest, Surf_spco2_mean_uatm_closest, bottomT_mean_C_closest,
                                           PBS443_closest))
cor_matrix2 <- cor(Env_corr_reduced, use = "complete.obs", method = "pearson")
print(cor_matrix2)

corrplot(cor_matrix2, 
         method = "color",
         type = "upper",
         order = "hclust",
         col = colorRampPalette(c("blue", "white", "red"))(200), # Custom color palette
         #addCoef.col = "black", # Add correlation coefficients
         tl.col = "black",    # Text label color
         tl.srt = 90,           # Text label rotation
         diag = FALSE           # Remove diagonal elements
)


# Get the upper triangle of the correlation matrix
upper_tri <- cor_matrix2
upper_tri[lower.tri(upper_tri, diag = TRUE)] <- NA
# Find variable pairs with correlation > 0.8
high_corr_pairs <- which(abs(upper_tri) > 0.8, arr.ind = TRUE)
# Create a data frame of the results
high_corr_df <- data.frame(
  Var1 = rownames(upper_tri)[high_corr_pairs[, 1]],
  Var2 = colnames(upper_tri)[high_corr_pairs[, 2]],
  Correlation = upper_tri[high_corr_pairs]
)
#Order the data frame by absolute correlation value, descending
high_corr_df <- high_corr_df[order(-abs(high_corr_df$Correlation)), ]
# View the sorted result
print(high_corr_df)

#Note: phosphate and mean pH are highly correlated (0.897) none of the other retained variables are
