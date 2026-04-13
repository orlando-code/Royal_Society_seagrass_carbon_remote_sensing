from __future__ import annotations

import time
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import box
from tqdm.auto import tqdm


def log_step(msg: str) -> None:
    print(f"\n[{time.strftime('%H:%M:%S')}] {msg}")


def ensure_same_crs(left_gdf: gpd.GeoDataFrame, right_gdf: gpd.GeoDataFrame, right_name: str = "right_gdf") -> gpd.GeoDataFrame:
    if left_gdf.crs is None:
        raise ValueError("Left GeoDataFrame has no CRS.")
    if right_gdf.crs is None:
        raise ValueError(f"{right_name} has no CRS.")
    if left_gdf.crs != right_gdf.crs:
        right_gdf = right_gdf.to_crs(left_gdf.crs)
    return right_gdf


def sjoin_chunked(
    left_gdf: gpd.GeoDataFrame,
    right_gdf: gpd.GeoDataFrame,
    *,
    how: str = "left",
    predicate: str = "intersects",
    chunk_size: int = 5000,
    desc: str = "sjoin",
) -> gpd.GeoDataFrame:
    n = len(left_gdf)
    if n == 0:
        return left_gdf.copy()

    out_parts = []
    idx = np.arange(0, n, chunk_size)
    for start in tqdm(idx, desc=desc):
        chunk = left_gdf.iloc[start : start + chunk_size]
        joined = gpd.sjoin(chunk, right_gdf, how=how, predicate=predicate)
        out_parts.append(joined)

    out = pd.concat(out_parts, axis=0)
    return gpd.GeoDataFrame(out, geometry=left_gdf.geometry.name, crs=left_gdf.crs)


def prefilter_eez_to_study_area(eez_gdf: gpd.GeoDataFrame, *study_gdfs: gpd.GeoDataFrame, buffer_deg: float = 0.0) -> gpd.GeoDataFrame:
    bounds = np.array([g.total_bounds for g in study_gdfs if len(g) > 0])
    minx, miny = bounds[:, 0].min(), bounds[:, 1].min()
    maxx, maxy = bounds[:, 2].max(), bounds[:, 3].max()

    if buffer_deg > 0:
        minx -= buffer_deg
        miny -= buffer_deg
        maxx += buffer_deg
        maxy += buffer_deg

    roi = box(minx, miny, maxx, maxy)
    cand_idx = list(eez_gdf.sindex.intersection(roi.bounds))
    eez_small = eez_gdf.iloc[cand_idx]
    eez_small = eez_small[eez_small.intersects(roi)].copy()
    return eez_small

