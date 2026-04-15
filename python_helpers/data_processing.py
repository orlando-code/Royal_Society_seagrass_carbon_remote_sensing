from __future__ import annotations

from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import xarray as xa
from rasterio.features import rasterize
from rasterio.transform import from_origin
from shapely import wkt
from tqdm.auto import tqdm

from . import spatial_join, stocks


def load_or_prepare_eov_geometries(
    *,
    all_poly_fp: Path,
    processed_polygon_fp: Path,
    processed_points_fp: Path,
    poly_cache_fp: Path,
    points_cache_fp: Path,
    inv_species_index_map: dict,
    lons=(-20, 30),
    lats=(30, 60),
):
    """Load or prepare seagrass EOV geometries from CSV files and cache them as GeoJSON files.

    Args:
        all_poly_fp (Path): Path to the CSV file containing the seagrass EOV geometries.
        processed_polygon_fp (Path): Path to the CSV file containing the processed seagrass EOV geometries.
        processed_points_fp (Path): Path to the CSV file containing the processed seagrass EOV points.
        poly_cache_fp (Path): Path to the GeoJSON file to cache the processed seagrass EOV geometries.
        points_cache_fp (Path): Path to the GeoJSON file to cache the processed seagrass EOV points.
        inv_species_index_map (dict): Dictionary mapping seagrass species indices to seagrass species names.
        lons (tuple, optional): Longitude range to crop the seagrass EOV geometries to. Defaults to (-20, 30).
        lats (tuple, optional): Latitude range to crop the seagrass EOV geometries to. Defaults to (30, 60).

    Returns:
        gpd.GeoDataFrame: GeoDataFrame containing the processed seagrass EOV geometries.
        gpd.GeoDataFrame: GeoDataFrame containing the processed seagrass EOV points.
    """
    # Ensure cache parent directories exist before any reads/writes.
    poly_cache_fp.parent.mkdir(parents=True, exist_ok=True)
    points_cache_fp.parent.mkdir(parents=True, exist_ok=True)

    if poly_cache_fp.exists() and points_cache_fp.exists():
        try:
            return gpd.read_parquet(poly_cache_fp), gpd.read_parquet(points_cache_fp)
        except (FileNotFoundError, OSError):
            # If a stale/partial cache exists, rebuild from source CSVs.
            pass

    poly_df = pd.read_csv(all_poly_fp)
    poly_df["geometry"] = poly_df["geom"].apply(wkt.loads)
    poly_gdf = gpd.GeoDataFrame(poly_df, geometry="geometry").drop(columns=["geom"])
    poly_gdf.set_crs(crs="epsg:3857", inplace=True)
    poly_gdf.to_crs(epsg=4326, inplace=True)
    poly_gdf = poly_gdf.cx[min(lons) : max(lons), min(lats) : max(lats)]

    processed_polygon_df = pd.read_csv(processed_polygon_fp, low_memory=False)
    id_and_species = processed_polygon_df[["id", "seagrass_species"]]
    poly_gdf = poly_gdf.merge(id_and_species, on="id", how="left")
    poly_gdf = poly_gdf[poly_gdf["seagrass_species"].notna()].copy()
    poly_gdf["prediction_species_index"] = poly_gdf["seagrass_species"].map(
        inv_species_index_map
    )

    points_df = pd.read_csv(processed_points_fp)
    points_df["geometry"] = points_df["geom"].apply(wkt.loads)
    points_gdf = gpd.GeoDataFrame(points_df, geometry="geometry").drop(columns=["geom"])
    points_gdf.set_crs(crs="epsg:3857", inplace=True)
    points_gdf.to_crs(epsg=4326, inplace=True)
    points_gdf["lat"] = points_gdf["geometry"].y
    points_gdf["lon"] = points_gdf["geometry"].x
    points_gdf = points_gdf[points_gdf["seagrass_species"].notna()].copy()
    points_gdf["prediction_species_index"] = points_gdf["seagrass_species"].map(
        inv_species_index_map
    )

    poly_gdf.to_parquet(poly_cache_fp, index=False)
    points_gdf.to_parquet(points_cache_fp, index=False)
    return poly_gdf, points_gdf


def attribute_geometries_from_prediction_grid(
    poly_gdf: gpd.GeoDataFrame,
    points_gdf: gpd.GeoDataFrame,
    pred_ds: xa.Dataset,
    *,
    req_col: str = "prediction_species_index",
    polygon_stat: str = "mean",
    all_touched: bool = True,
):
    """Attribute polygon and point geometries from a species-specific prediction grid.

    Args:
        poly_gdf (gpd.GeoDataFrame): Polygon geometries with species index column.
        points_gdf (gpd.GeoDataFrame): Point geometries with species index column.
        pred_ds (xa.Dataset): Prediction dataset with ``prediction`` and ``prediction_se``.
        req_col (str, optional): Column containing species indices. Defaults to ``"prediction_species_index"``.
        polygon_stat (str, optional): Polygon aggregation statistic. Only ``"mean"`` is supported. Defaults to ``"mean"``.
        all_touched (bool, optional): Rasterization mode for polygon attribution. Defaults to ``True``.

    Returns:
        gpd.GeoDataFrame: Polygon GeoDataFrame with ``gC_cm3`` and ``gC_cm3_se`` columns.
        gpd.GeoDataFrame: Point GeoDataFrame with ``gC_cm3`` and ``gC_cm3_se`` columns.

    Raises:
        KeyError: If ``req_col`` is missing from either input GeoDataFrame.
        ValueError: If ``polygon_stat`` is not supported.
    """
    if req_col not in poly_gdf.columns or req_col not in points_gdf.columns:
        raise KeyError(
            "Both polygon and point GeoDataFrames must contain `prediction_species_index`."
        )

    gdf_local = poly_gdf.copy()
    pts_local = points_gdf.copy()
    if gdf_local.crs is None:
        gdf_local = gdf_local.set_crs("EPSG:4326")
    if pts_local.crs is None:
        pts_local = pts_local.set_crs("EPSG:4326")
    if str(gdf_local.crs).lower() not in {"epsg:4326", "4326"}:
        gdf_local = gdf_local.to_crs(4326)
    if str(pts_local.crs).lower() not in {"epsg:4326", "4326"}:
        pts_local = pts_local.to_crs(4326)

    gdf_local = gdf_local[
        gdf_local.geometry.notna() & ~gdf_local.geometry.is_empty
    ].copy()
    pts_local = pts_local[
        pts_local.geometry.notna() & ~pts_local.geometry.is_empty
    ].copy()
    gdf_local = gdf_local[gdf_local[req_col].notna()].copy()
    pts_local = pts_local[pts_local[req_col].notna()].copy()
    gdf_local[req_col] = gdf_local[req_col].astype(int)
    pts_local[req_col] = pts_local[req_col].astype(int)

    lon_asc = np.sort(np.asarray(pred_ds["lon"].values))
    lat_asc = np.sort(np.asarray(pred_ds["lat"].values))
    ny, nx = len(lat_asc), len(lon_asc)

    pred_all = pred_ds["prediction"].sel(lon=lon_asc, lat=lat_asc)
    se_all = pred_ds["prediction_se"].sel(lon=lon_asc, lat=lat_asc)
    fill = pred_all.attrs.get("_FillValue", -9999)
    pred_all = pred_all.where(pred_all != fill)
    se_all = se_all.where(se_all != fill)

    species_coord = np.asarray(pred_all["prediction_species_index"].values)
    species_set = set(int(v) for v in species_coord)

    lon_step = float(np.median(np.diff(lon_asc)))
    lat_step = float(np.median(np.diff(lat_asc)))
    lon_edges = np.r_[lon_asc - lon_step / 2, lon_asc[-1] + lon_step / 2]
    lat_edges = np.r_[lat_asc - lat_step / 2, lat_asc[-1] + lat_step / 2]

    pred_by_species = {
        int(sp): pred_all.sel(prediction_species_index=int(sp))
        .transpose("lat", "lon")
        .values
        for sp in species_set
    }
    se_by_species = {
        int(sp): se_all.sel(prediction_species_index=int(sp))
        .transpose("lat", "lon")
        .values
        for sp in species_set
    }

    poly_values = np.full(len(gdf_local), np.nan, dtype=float)
    poly_se_values = np.full(len(gdf_local), np.nan, dtype=float)
    for row_pos, (_, row) in tqdm(
        enumerate(gdf_local.iterrows()), total=len(gdf_local), desc="Polygons"
    ):
        sp = int(row[req_col])
        if sp not in pred_by_species or sp not in se_by_species:
            continue

        xmin, ymin, xmax, ymax = row.geometry.bounds
        ix0 = max(0, np.searchsorted(lon_edges, xmin, side="right") - 1)
        ix1 = min(nx - 1, np.searchsorted(lon_edges, xmax, side="left"))
        iy0 = max(0, np.searchsorted(lat_edges, ymin, side="right") - 1)
        iy1 = min(ny - 1, np.searchsorted(lat_edges, ymax, side="left"))
        if ix1 < ix0 or iy1 < iy0:
            continue

        win_w = ix1 - ix0 + 1
        win_h = iy1 - iy0 + 1
        win_transform = from_origin(
            west=float(lon_edges[ix0]),
            north=float(lat_edges[iy1 + 1]),
            xsize=float(lon_step),
            ysize=float(lat_step),
        )
        mask_desc = rasterize(
            [(row.geometry, 1)],
            out_shape=(win_h, win_w),
            transform=win_transform,
            fill=0,
            all_touched=all_touched,
            dtype="uint8",
        )
        mask = np.flipud(mask_desc).astype(bool)
        if not mask.any():
            continue

        pred_window = pred_by_species[sp][iy0 : iy1 + 1, ix0 : ix1 + 1]
        se_window = se_by_species[sp][iy0 : iy1 + 1, ix0 : ix1 + 1]
        vals = pred_window[mask]
        se_vals = se_window[mask]
        vals = vals[np.isfinite(vals)]
        se_vals = se_vals[np.isfinite(se_vals)]
        if vals.size == 0 or se_vals.size == 0:
            continue

        if polygon_stat == "mean":
            poly_values[row_pos] = float(vals.mean())
            poly_se_values[row_pos] = float(se_vals.mean())
        else:
            raise ValueError("Only polygon_stat='mean' is implemented.")

    point_values = np.full(len(pts_local), np.nan, dtype=float)
    point_se_values = np.full(len(pts_local), np.nan, dtype=float)
    x = pts_local.geometry.x.to_numpy()
    y = pts_local.geometry.y.to_numpy()
    sp_arr = pts_local[req_col].to_numpy(dtype=int)
    cols = np.searchsorted(lon_edges, x, side="right") - 1
    rows = np.searchsorted(lat_edges, y, side="right") - 1
    in_bounds = (cols >= 0) & (cols < nx) & (rows >= 0) & (rows < ny)
    valid_sp = np.isin(sp_arr, list(species_set))
    valid = in_bounds & valid_sp

    for sp in tqdm(np.unique(sp_arr[valid]), desc="Point species groups"):
        idx = np.where(valid & (sp_arr == sp))[0]
        if idx.size == 0:
            continue
        arr_pred = pred_by_species[int(sp)]
        arr_se = se_by_species[int(sp)]
        point_values[idx] = arr_pred[rows[idx], cols[idx]]
        point_se_values[idx] = arr_se[rows[idx], cols[idx]]

    gdf_attributed = gdf_local.copy()
    pts_attributed = pts_local.copy()
    gdf_attributed["gC_cm3"] = poly_values
    gdf_attributed["gC_cm3_se"] = poly_se_values
    pts_attributed["gC_cm3"] = point_values
    pts_attributed["gC_cm3_se"] = point_se_values
    return gdf_attributed, pts_attributed


def load_or_build_attributed_geometries(
    *,
    pred_ds: xa.Dataset,
    attributed_polygon_fp: Path,
    attributed_points_fp: Path,
    all_poly_fp: Path,
    processed_polygon_fp: Path,
    processed_points_fp: Path,
    poly_cache_fp: Path,
    points_cache_fp: Path,
    inv_species_index_map: dict,
):
    """Load cached attributed geometries or build them from source files.

    Args:
        pred_ds (xa.Dataset): Prediction dataset used for geometry attribution.
        attributed_polygon_fp (Path): Output/cache path for attributed polygons.
        attributed_points_fp (Path): Output/cache path for attributed points.
        all_poly_fp (Path): Path to raw polygon CSV.
        processed_polygon_fp (Path): Path to processed polygon CSV with species labels.
        processed_points_fp (Path): Path to processed point CSV with species labels.
        poly_cache_fp (Path): Path to cached polygon GeoJSON prepared from CSV.
        points_cache_fp (Path): Path to cached point GeoJSON prepared from CSV.
        inv_species_index_map (dict): Mapping from species label to model species index.

    Returns:
        gpd.GeoDataFrame: Attributed polygon GeoDataFrame.
        gpd.GeoDataFrame: Attributed point GeoDataFrame.
    """
    attributed_polygon_fp.parent.mkdir(parents=True, exist_ok=True)
    attributed_points_fp.parent.mkdir(parents=True, exist_ok=True)

    if attributed_polygon_fp.exists() and attributed_points_fp.exists():
        try:
            return gpd.read_parquet(attributed_polygon_fp), gpd.read_parquet(
                attributed_points_fp
            )
        except (FileNotFoundError, OSError):
            # Rebuild attributed layers if cache read fails.
            pass

    poly_gdf, points_gdf = load_or_prepare_eov_geometries(
        all_poly_fp=all_poly_fp,
        processed_polygon_fp=processed_polygon_fp,
        processed_points_fp=processed_points_fp,
        poly_cache_fp=poly_cache_fp,
        points_cache_fp=points_cache_fp,
        inv_species_index_map=inv_species_index_map,
    )
    gdf_attributed, pts_attributed = attribute_geometries_from_prediction_grid(
        poly_gdf, points_gdf, pred_ds
    )
    gdf_attributed.to_parquet(attributed_polygon_fp, index=False)
    pts_attributed.to_parquet(attributed_points_fp, index=False)
    return gdf_attributed, pts_attributed


def build_or_load_eez_joins(
    *,
    gdf_attributed: gpd.GeoDataFrame,
    pts_attributed: gpd.GeoDataFrame,
    eez_fp: Path,
    polygon_join_fp: Path,
    point_join_fp: Path,
    buffer_deg: float = 0.2,
    polygon_chunk_size: int = 3000,
    point_chunk_size: int = 10000,
):
    """Build or load cached EEZ joins for attributed polygons and points.

    Args:
        gdf_attributed (gpd.GeoDataFrame): Attributed polygon geometries.
        pts_attributed (gpd.GeoDataFrame): Attributed point geometries.
        eez_fp (Path): Path to EEZ vector data.
        polygon_join_fp (Path): Cache path for polygon-EEZ join parquet.
        point_join_fp (Path): Cache path for point-EEZ join parquet.
        buffer_deg (float, optional): Buffer in degrees around study bounds for EEZ prefiltering. Defaults to 0.2.
        polygon_chunk_size (int, optional): Chunk size for polygon spatial joins. Defaults to 3000.
        point_chunk_size (int, optional): Chunk size for point spatial joins. Defaults to 10000.

    Returns:
        gpd.GeoDataFrame: Polygon GeoDataFrame joined with EEZ attributes.
        gpd.GeoDataFrame: Point GeoDataFrame joined with EEZ attributes.
    """
    t0 = spatial_join.time.time()
    spatial_join.log_step("Reading EEZ layer...")
    eez_gdf = gpd.read_file(eez_fp, engine="pyogrio")
    spatial_join.log_step(
        f"EEZ read complete: {len(eez_gdf):,} rows in {spatial_join.time.time() - t0:.1f}s"
    )

    eez_cols = ["geometry", "GEONAME", "TERRITORY1", "SOVEREIGN1"]
    eez_gdf = eez_gdf[eez_cols].copy()
    eez_gdf = spatial_join.ensure_same_crs(gdf_attributed, eez_gdf, "eez_gdf")
    if pts_attributed.crs != gdf_attributed.crs:
        pts_attributed = pts_attributed.to_crs(gdf_attributed.crs)

    t0 = spatial_join.time.time()
    spatial_join.log_step("Prefiltering EEZ to study-area bbox...")
    eez_small = spatial_join.prefilter_eez_to_study_area(
        eez_gdf, gdf_attributed, pts_attributed, buffer_deg=buffer_deg
    )
    spatial_join.log_step(
        f"EEZ prefilter complete: {len(eez_small):,} rows in {spatial_join.time.time() - t0:.1f}s"
    )

    if polygon_join_fp.exists():
        spatial_join.log_step("Loading joined polygon data from cache...")
        gdf_attributed_joined = gpd.read_parquet(polygon_join_fp)
    else:
        t0 = spatial_join.time.time()
        spatial_join.log_step("Intersecting polygons with EEZs (chunked)...")
        poly_pts = gdf_attributed.copy()
        poly_pts["geometry"] = poly_pts.representative_point()
        joined_poly_pts = spatial_join.sjoin_chunked(
            poly_pts,
            eez_small,
            how="left",
            predicate="within",
            chunk_size=polygon_chunk_size,
            desc="polygon sjoin",
        )

        # Keep original polygon geometries while using representative points only
        # as spatial keys for EEZ matching.
        gdf_attributed_joined = gdf_attributed.loc[joined_poly_pts.index].copy()
        eez_attr_cols = [c for c in eez_cols if c != "geometry"]
        for col in eez_attr_cols:
            gdf_attributed_joined[col] = (
                joined_poly_pts[col].values
                if col in joined_poly_pts.columns
                else np.nan
            )

        gdf_attributed_joined.to_parquet(polygon_join_fp, index=False)
        spatial_join.log_step(
            f"Polygon join + save done in {spatial_join.time.time() - t0:.1f}s"
        )

    if point_join_fp.exists():
        spatial_join.log_step("Loading joined point data from cache...")
        pts_attributed_joined = gpd.read_parquet(point_join_fp)
    else:
        t0 = spatial_join.time.time()
        spatial_join.log_step("Intersecting points with EEZs (chunked)...")
        pts_attributed_joined = spatial_join.sjoin_chunked(
            pts_attributed,
            eez_small,
            how="left",
            predicate="within",
            chunk_size=point_chunk_size,
            desc="point sjoin",
        )
        drop_cols = [
            c
            for c in pts_attributed_joined.columns
            if c.startswith("index_right") or c == "geometry_right"
        ]
        pts_attributed_joined = pts_attributed_joined.drop(
            columns=drop_cols, errors="ignore"
        )
        pts_attributed_joined.to_parquet(point_join_fp, index=False)
        spatial_join.log_step(
            f"Point join + save done in {spatial_join.time.time() - t0:.1f}s"
        )

    spatial_join.log_step("All joins complete.")
    return gdf_attributed_joined, pts_attributed_joined


def prepare_european_stock_tables(
    gdf_attributed_joined: gpd.GeoDataFrame,
    pts_attributed_joined: gpd.GeoDataFrame,
    *,
    excluded_sovereign: str | list[str] = "Tunisia",
    equal_area_epsg: int = 3035,
):
    """Prepare cleaned stock-analysis GeoDataFrames for European summaries.

    Args:
        gdf_attributed_joined (gpd.GeoDataFrame): Polygon data joined to EEZ attributes.
        pts_attributed_joined (gpd.GeoDataFrame): Point data joined to EEZ attributes.
        excluded_sovereign (str, optional): Sovereign area to exclude from polygon summaries. Defaults to ``"Tunisia"``.
        equal_area_epsg (int, optional): EPSG code for equal-area reprojection. Defaults to 3035.

    Returns:
        gpd.GeoDataFrame: Filtered polygon GeoDataFrame in equal-area CRS.
        gpd.GeoDataFrame: Filtered point GeoDataFrame with join artifact columns removed.
    """
    if isinstance(excluded_sovereign, str):
        excluded_sovereign = [excluded_sovereign]
    gdf_filtered = gdf_attributed_joined[
        ~gdf_attributed_joined["SOVEREIGN1"].isin(excluded_sovereign)
    ].copy()
    pts_filtered = pts_attributed_joined.drop(
        columns=[
            c
            for c in pts_attributed_joined.columns
            if c.startswith("index_right") or c == "geometry_right"
        ],
        errors="ignore",
    )
    ea_gdf = gdf_filtered.to_crs(epsg=equal_area_epsg)
    return ea_gdf, pts_filtered


def add_stock_columns(
    ea_gdf: gpd.GeoDataFrame, *, depth_cm: float = 100.0, clip_gC_cm3=(0.0, 0.25)
):
    """Add area and carbon stock summary columns to an equal-area GeoDataFrame.

    Args:
        ea_gdf (gpd.GeoDataFrame): Equal-area polygon GeoDataFrame with carbon density columns.
        depth_cm (float, optional): Sediment depth used for stock conversion. Defaults to 100.0.
        clip_gC_cm3 (tuple, optional): Min/max bounds applied to ``gC_cm3`` before total stock calculation. Defaults to ``(0.0, 0.25)``.

    Returns:
        gpd.GeoDataFrame: Copy of input with ``area_ha``, ``stock_kg_ha``, ``stock_kg_ha_se``, ``stock_Gg``, and ``stock_Gg_se`` columns.
    """
    out = ea_gdf.copy()
    out["area_ha"] = out.geometry.area / 1e4
    # Clip density before all stock derivations so areal stock and total stock match.
    out["gC_cm3"] = out["gC_cm3"].clip(*clip_gC_cm3)
    out["stock_kg_ha"] = stocks.convert_carbon_density_to_carbon_areal_stock(
        out["gC_cm3"], depth_cm
    )
    out["stock_kg_ha_se"] = stocks.convert_carbon_density_to_carbon_areal_stock(
        out["gC_cm3_se"], depth_cm
    )
    out["stock_Gg"] = stocks.convert_carbon_density_to_carbon_stock(
        out["gC_cm3"], depth_cm, out.geometry.area
    )
    out["stock_Gg_se"] = stocks.convert_carbon_density_to_carbon_stock(
        out["gC_cm3_se"], depth_cm, out.geometry.area
    )
    return out


def load_gomis_table(gomis_fp: Path):
    """Load and clean the Gomis comparison table.

    Args:
        gomis_fp (Path): Path to the Gomis CSV file.

    Returns:
        pd.DataFrame: Cleaned Gomis table with numeric comparison columns coerced to floats.
    """
    gomis_df = pd.read_csv(gomis_fp).dropna(axis=0, how="all")
    cols_to_convert = [
        "Seagrass extent (ha)",
        "Biomass carbon stock (kg C / ha)",
        "Biomass carbon stock (Gg)",
    ]
    for col in cols_to_convert:
        if col in gomis_df:
            gomis_df[col] = (
                gomis_df[col]
                .astype(str)
                .str.replace(",", "", regex=False)
                .replace("nan", np.nan)
                .astype(float)
            )
    return gomis_df


def build_national_metrics(ea_gdf: gpd.GeoDataFrame):
    """Aggregate national stock metrics from polygon-level data.

    Args:
        ea_gdf (gpd.GeoDataFrame): Equal-area polygon GeoDataFrame with stock columns.

    Returns:
        pd.DataFrame: National summary table keyed by ``Country``.
    """
    return (
        ea_gdf.groupby("TERRITORY1")
        .apply(
            lambda g: pd.Series(
                {
                    "area_ha": g["area_ha"].sum(),
                    "stock_kg_ha": stocks.safe_weighted_average(
                        g["stock_kg_ha"], g["area_ha"]
                    ),
                    "stock_kg_ha_se": stocks.weighted_mean_standard_error(
                        g["stock_kg_ha_se"],
                        g["area_ha"],
                        g["stock_kg_ha"],
                    ),
                    "stock_Gg": g["stock_Gg"].sum(),
                    # Propagate SE for summed stock using root-sum-of-squares
                    # (assumes polygon-level errors are independent).
                    "stock_Gg_se": float(
                        np.sqrt(
                            np.nansum(np.square(g["stock_Gg_se"].to_numpy(dtype=float)))
                        )
                    ),
                }
            )
        )
        .reset_index()
        .rename(columns={"TERRITORY1": "Country"})
    )


def build_comparison_dataframe(ea_gdf: gpd.GeoDataFrame, gomis_fp: Path):
    """Build merged comparison dataframe for plotting study vs Gomis estimates.

    Args:
        ea_gdf (gpd.GeoDataFrame): Equal-area polygon GeoDataFrame with computed stock columns.
        gomis_fp (Path): Path to the Gomis CSV comparison table.

    Returns:
        pd.DataFrame: Outer-joined comparison table containing Gomis and this-study metrics.
    """
    gomis_df = load_gomis_table(gomis_fp)
    national_metrics = build_national_metrics(ea_gdf)
    return pd.merge(gomis_df, national_metrics, on="Country", how="outer")
