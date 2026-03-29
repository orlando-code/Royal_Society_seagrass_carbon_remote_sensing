# Plot configuration: human-friendly labels for covariates used in figures.
# 
# If a variable is not in VAR_LABELS, label_vars() simply returns the
# original name, so it is safe to use everywhere.


REGION_COLOURS <- c(
  "Mediterranean Sea" = "#9EF8EE",
  "North European Atlantic" = "#F24405",
  "South European Atlantic" = "#FA7F08",
  "Baltic Sea" = "#22bb52",
  "Black Sea" = "#348888"
)


SPECIES_COLOURS <- c(
    "Cymodocea nodosa" = "#3366CC", # blue
    "Halophila stipulacea" = "#FF9900", # orange
    "Posidonia oceanica" = "#00CED1", # dark cyan
    "Zostera marina" = "#32CD32", # lime green
    "Zostera noltei" = "#CCCC00", # yellow-green
    "Zostera marina and Zostera noltei" = "#006400", # dark green
    "Zostera marina and Cymodocea nodosa" = "#8B0000", # dark red
    "Unspecified" = "#808080" # gray
  )


MODEL_COLOURS <- c(
  "Ordinary Kriging" = "#08306b", # Dark blue
  "Universal Kriging (env drift)" = "#3182bd", # Light blue
  "GAM" = "#99000d", # Dark red
  "INLA" = "#cb181d", # Light red
  "XGBoost" = "#006d2c", # Dark green
  "XGB" = "#006d2c", # Alias used by run_cv (same as XGBoost)
  "Random Forest" = "#41ab5d", # Light green
  "RF" = "#41ab5d", # Alias used by run_cv (same as Random Forest)
  "Neural Network" = "#54278f", # Dark purple
  "NN" = "#54278f", # Alias used by run_cv (same as Neural Network)
  "SVM" = "#9e9ac8", # Light purple
  "GPR" = "#fd8d3c" # Orange
)

# Standard line styles for pooled vs mean-fold metrics across plots.
METRIC_LINESTYLES <- c(
  "Pooled R2" = "solid",
  "Mean-fold R2" = "22",
  "Pooled RMSE" = "solid",
  "Mean-fold RMSE" = "22"
)


VAR_LABELS <- c(
  # Response variable
  median_carbon_density_100cm = "Median Carbon Density (g C/cm³)",
  # Bathymetry / depth
  cmems_mod_glo_phy_my_0.083deg_static_deptho_9.00w_34.00e_34.00n_61.00n =
    "Bathymetry (depth, m)",
  gebco_2025_n61.0_s34.0_w_10.0_e35.0 = "Bathymetry (depth, m)",

  # Bottom temperature
  bottomt_mean_daily_c = "Bottom temperature (mean daily, °C)",
  bottomt_p95_daily_c  = "Bottom temperature (95th percentile daily, °C)",

  # Waves
  daily_max_wave_height_mean_m = "Max. significant wave height (daily mean, m)",
  daily_mean_wave_height_mean_m = "Significant wave height (daily mean, m)",
  daily_std_wave_height_mean_m = "Significant wave height (daily SD, m)",
  no_sig_wave_heights_mean = "Number of significant wave heights (mean)",
  sw_ke_mean_cm2_s2 = "Surface wave kinetic energy (mean, cm²/s²)",
  sw_ke_p95_cm2_s2  = "Surface wave kinetic energy (95th percentile, cm²/s²)",
  wave_height_mean_m = "Significant wave height (mean, m)",
  wave_height_p95_m  = "Significant wave height (95th percentile, m)",

  # Nutrients / biogeochemistry
  fe_mean_monthly_1.5m_mmol_m3 = "Dissolved iron (mean, mmol/m³)",
  no3_mean_monthly_1.5m_mmol_m3 = "Nitrate (mean, mmol/m³)",
  nppv_mean_monthly_0.5m_mg_m3_day = "Net primary production (mean, mg C m³/d)",
  ox_mean_daily_2.7m_mmol_m3 = "Dissolved oxygen (mean daily, mmol/m³)",
  ox_p05_daily_2.7m_mmol_m3  = "Dissolved oxygen (5th percentile daily, mmol/m³)",
  ph_mean_monthly_1.5m = "pH (mean, surface)",
  ph_p05_monthly_1.5m  = "pH (5th percentile, surface)",
  phyc_mean_monthly_0.5m_mmol_m3 = "Phytoplanktic carbon (mean, mmol/m³)",
  po4_mean_monthly_1.5m_mmol_m3  = "Phosphate (mean, mmol/m³)",

  # Ocean colour / optical properties (Sentinel-3 OLCI)
  "s3a_olci_errnt.20160425_20241231.l3m.cu.chl.chlor_a.4km" =
    "Chlorophyll-a (OLCI, 4 km)",
  "s3a_olci_errnt.20160425_20241231.l3m.cu.iop.adg_443.4km" =
    "Adg(443) (OLCI, 4 km)",
  "s3a_olci_errnt.20160425_20241231.l3m.cu.iop.bbp_443.4km" =
    "bbp(443) (OLCI, 4 km)",
  "s3a_olci_errnt.20160425_20241231.l3m.cu.kd.kd_490.4km" =
    "Kd(490) diffuse attenuation (OLCI, 4 km)",
  "s3a_olci_errnt.20160425_20241231.l3m.cu.rrs.rrs_443.4km" =
    "Rrs(443) (OLCI, 4 km)",
  "s3a_olci_errnt.20160425_20241231.l3m.cu.rrs.rrs_620.4km" =
    "Rrs(620) (OLCI, 4 km)",

  # Salinity
  sal_mean_monthly_0.5m = "Sea surface salinity (mean)",
  sal_p05_monthly_0.5m  = "Sea surface salinity (5th percentile)",

  # Sea surface temperature
  sst_daily_mean_k = "Sea surface temperature (mean daily, K)",
  sst_daily_p95_k  = "Sea surface temperature (95th percentile daily, K)",

  # CO2 fluxes and related
  spco2_monthly_mean_pa = "Surface pCO2 (mean, Pa)",
  spco2_monthly_p95_pa  = "Surface pCO2 (95th percentile, Pa)",
  surf_fgco2_mean_molc_m2_yr = "Air–sea CO2 flux (mean, mol C/m²yr)",
  surf_fgco2_p95_molc_m2_yr  = "Air–sea CO2 flux (95th percentile, mol C/m²yr)",
  surf_spco2_mean_uatm = "Surface pCO2 (mean, µatm)",
  surf_spco2_p95_uatm  = "Surface pCO2 (95th percentile, µatm)",

  # Alkalinity
  surf_talk_mean_umol_kg = "Total alkalinity (mean, µmol/kg)",
  surf_talk_p95_umol_kg  = "Total alkalinity (95th percentile, µmol/kg)",

  # Currents
  uo_mean_1.5m_m_s = "Zonal current (mean, m/s)",
  uo_p90_1.5m_m_s  = "Zonal current (90th percentile, m/s)",
  vo_mean_1.5m_m_s = "Meridional current (mean, m/s)",
  vo_p90_1.5m_m_s  = "Meridional current (90th percentile, m/s)"
)

# Model colours for CV and tuning plots (used by cv_pipeline_plots.R and cv_pipeline.R)
model_colours <- c(
  "Ordinary Kriging" = "#08306b",
  "Universal Kriging (env drift)" = "#3182bd",
  "GAM" = "#99000d",
  "INLA" = "#cb181d",
  "XGBoost" = "#006d2c",
  "XGB" = "#006d2c",
  "Random Forest" = "#41ab5d",
  "RF" = "#41ab5d",
  "Neural Network" = "#54278f",
  "NN" = "#54278f",
  "SVM" = "#9e9ac8",
  "GPR" = "#fd8d3c"
)

# Helper to map variable names to labels, falling back to the original
# name if no mapping exists.
label_vars <- function(x) {
  labs <- VAR_LABELS[x]
  labs[is.na(labs)] <- x[is.na(labs)]
  labs
}