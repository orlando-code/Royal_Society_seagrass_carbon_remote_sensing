# Plot configuration: human-friendly labels for covariates used in figures.
# 
# If a variable is not in VAR_LABELS, label_vars() simply returns the
# original name, so it is safe to use everywhere.

EUROPE_LAT_RANGE <- c(34.65, 60.11)
EUROPE_LON_RANGE <- c(-8.85, 33.30)

# Colours for region shapes in maps.
REGION_COLOURS <- c(
  "Mediterranean Sea" = "#9EF8EE",
  "North European Atlantic" = "#F24405",
  "South European Atlantic" = "#FA7F08",
  "Baltic Sea" = "#22bb52",
  "Black Sea" = "#348888"
)

# Colours for seagrass species in maps.
SPECIES_COLOURS <- c(
    "Cymodocea nodosa" = "#1b9e77", # blue
    "Halophila stipulacea" = "#d95f02", # orange
    "Posidonia oceanica" = "#66a61e", # dark cyan
    "Zostera marina" = "#a6761d", # lime green
    "Zostera noltei" = "#e7298a", # yellow-green
    "Zostera marina and Zostera noltei" = "#666666", # dark green
    "Zostera marina and Cymodocea nodosa" = "#e6ab02", # dark red
    "Unspecified" = "black" # gray
  )

# Colours for machine learning models in plots.
MODEL_COLOURS <- c(
  "LR" = "#3182bd", # Light blue
  "GAM" = "#99000d", # Dark red
  "XGB" = "#006d2c", # Alias used by run_cv (same as XGBoost)
  "RF" = "#50f01b", # Alias used by run_cv (same as Random Forest)
  "GPR" = "#fd8d3c" # Orange
)

# Standard line styles for pooled vs mean-fold metrics across plots.
METRIC_LINESTYLES <- c(
  "Pooled R2" = "solid",
  "Mean-fold R2" = "22",
  "Pooled RMSE" = "solid",
  "Mean-fold RMSE" = "22"
)

# Human-friendly labels for covariates.
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


# Helper to map variable names to labels, falling back to the original
# name if no mapping exists.
label_vars <- function(x) {
  labs <- VAR_LABELS[x]
  labs[is.na(labs)] <- x[is.na(labs)]
  labs
}

# -----------------------------------------------------------------------------
# Default ggplot theme (manuscript / supplement figures)
# -----------------------------------------------------------------------------

#' Minimal theme with title emphasis, legend at bottom, no minor grid.
#' @param base_size Base font size passed to [ggplot2::theme_minimal()].
theme_paper <- function(base_size = 11) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = base_size + 1),
      axis.title = ggplot2::element_text(size = base_size),
      legend.position = "bottom",
      legend.key.width = grid::unit(2, "cm"),
      panel.grid.minor = ggplot2::element_blank()
    )
}


# Default for all ggplot2 figures after this file is sourced.
if (requireNamespace("ggplot2", quietly = TRUE) && requireNamespace("grid", quietly = TRUE)) {
  ggplot2::theme_set(theme_paper())
}

