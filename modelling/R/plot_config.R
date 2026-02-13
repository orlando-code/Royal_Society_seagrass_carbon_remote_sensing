# Plot configuration: human-friendly labels for covariates used in figures.
# 
# Usage:
#   source("modelling/plot_config.R")
#   label_vars(c("bottomt_mean_daily_c", "wave_height_p95_m"))
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


VAR_LABELS <- c(
  # Bathymetry / depth
  cmems_mod_glo_phy_my_0.083deg_static_deptho_9.00w_34.00e_34.00n_61.00n =
    "Bathymetry (depth, m)",
  
  # Bottom temperature
  bottomt_mean_daily_c = "Bottom temperature (mean daily, °C)",
  bottomt_p95_daily_c  = "Bottom temperature (95th percentile daily, °C)",

  # Waves
  daily_max_wave_height_mean_m = "Max. significant wave height (daily mean, m)",
  daily_mean_wave_height_mean_m = "Significant wave height (daily mean, m)",
  daily_std_wave_height_mean_m = "Significant wave height (daily SD, m)",
  no_sig_wave_heights_mean = "Number of significant wave heights (mean)",
  sw_ke_mean_cm2_s2 = "Surface wave kinetic energy (mean, cm² s⁻²)",
  sw_ke_p95_cm2_s2  = "Surface wave kinetic energy (95th percentile, cm² s⁻²)",
  wave_height_mean_m = "Significant wave height (mean, m)",
  wave_height_p95_m  = "Significant wave height (95th percentile, m)",

  # Nutrients / biogeochemistry
  fe_mean_monthly_1.5m_mmol_m3 = "Dissolved iron (mean, mmol m⁻³)",
  no3_mean_monthly_1.5m_mmol_m3 = "Nitrate (mean, mmol m⁻³)",
  nppv_mean_monthly_0.5m_mg_m3_day = "Net primary production (mean, mg C m⁻³ d⁻¹)",
  ox_mean_daily_2.7m_mmol_m3 = "Dissolved oxygen (mean daily, mmol m⁻³)",
  ox_p05_daily_2.7m_mmol_m3  = "Dissolved oxygen (5th percentile daily, mmol m⁻³)",
  ph_mean_monthly_1.5m = "pH (mean, surface)",
  ph_p05_monthly_1.5m  = "pH (5th percentile, surface)",
  phyc_mean_monthly_0.5m_mmol_m3 = "Phytoplanktic carbon (mean, mmol m⁻³)",
  po4_mean_monthly_1.5m_mmol_m3  = "Phosphate (mean, mmol m⁻³)",

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
  spco2_monthly_mean_pa = "Surface pCO₂ (mean, Pa)",
  spco2_monthly_p95_pa  = "Surface pCO₂ (95th percentile, Pa)",
  surf_fgco2_mean_molc_m2_yr = "Air–sea CO₂ flux (mean, mol C m⁻² yr⁻¹)",
  surf_fgco2_p95_molc_m2_yr  = "Air–sea CO₂ flux (95th percentile, mol C m⁻² yr⁻¹)",
  surf_spco2_mean_uatm = "Surface pCO₂ (mean, µatm)",
  surf_spco2_p95_uatm  = "Surface pCO₂ (95th percentile, µatm)",

  # Alkalinity
  surf_talk_mean_umol_kg = "Total alkalinity (mean, µmol kg⁻¹)",
  surf_talk_p95_umol_kg  = "Total alkalinity (95th percentile, µmol kg⁻¹)",

  # Currents
  uo_mean_1.5m_m_s = "Zonal current (mean, m s⁻¹)",
  uo_p90_1.5m_m_s  = "Zonal current (90th percentile, m s⁻¹)",
  vo_mean_1.5m_m_s = "Meridional current (mean, m s⁻¹)",
  vo_p90_1.5m_m_s  = "Meridional current (90th percentile, m s⁻¹)"
)

# Model colours for CV and tuning plots (used by cv_pipeline_plots.R and cv_pipeline.R)
model_colours <- c(
  "Ordinary Kriging" = "#08306b",
  "Universal Kriging (env drift)" = "#3182bd",
  "GAM" = "#99000d",
  "INLA" = "#cb181d",
  "XGBoost" = "#006d2c",
  "Random Forest" = "#41ab5d",
  "Neural Network" = "#54278f",
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

