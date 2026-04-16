from __future__ import annotations

import pickle
from pathlib import Path

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
from cartopy.mpl.ticker import LatitudeFormatter, LongitudeFormatter
from shapely import make_valid
from shapely.geometry import box


def dynamic_buffer_distance(extent, frac=0.003, min_buf=0.0005, max_buf=0.06):
    """Compute a buffer distance in degrees from map extent size.

    Args:
        extent (tuple): Map extent as ``(xmin, xmax, ymin, ymax)`` in degrees.
        frac (float, optional): Fraction of the smaller span to use. Defaults to 0.003.
        min_buf (float, optional): Minimum buffer in degrees. Defaults to 0.0005.
        max_buf (float, optional): Maximum buffer in degrees. Defaults to 0.06.

    Returns:
        float: Clipped buffer distance in degrees.
    """
    xmin, xmax, ymin, ymax = extent
    scale = min(xmax - xmin, ymax - ymin)
    return float(np.clip(scale * frac, min_buf, max_buf))


def dynamic_point_size(extent, n_points, min_s=2.0, max_s=8.0):
    """Choose scatter point size from panel extent and point density.

    Args:
        extent (tuple): Map extent as ``(xmin, xmax, ymin, ymax)``.
        n_points (int): Number of points in the panel.
        min_s (float, optional): Minimum marker size. Defaults to 2.0.
        max_s (float, optional): Maximum marker size. Defaults to 8.0.

    Returns:
        float: Marker size for ``scatter``.
    """
    xmin, xmax, ymin, ymax = extent
    area = max((xmax - xmin) * (ymax - ymin), 1e-9)
    density = n_points / area
    base = 0.18 * np.sqrt(area) + 1.2
    density_factor = 1 / (1.0 + 0.02 * np.sqrt(density))
    return float(np.clip(base * density_factor, min_s, max_s))


def auto_simplify_tolerance(extent, fraction=0.00001):
    """Default polygon simplification tolerance from map span.

    Args:
        extent (tuple): Map extent as ``(xmin, xmax, ymin, ymax)``.
        fraction (float, optional): Fraction of the smaller span. Defaults to 0.00001.

    Returns:
        float: Simplification tolerance in coordinate units.
    """
    xmin, xmax, ymin, ymax = extent
    return min(xmax - xmin, ymax - ymin) * fraction


def auto_buffer_distance(extent, fraction=0.02):
    """Default geometry buffer distance from map span.

    Args:
        extent (tuple): Map extent as ``(xmin, xmax, ymin, ymax)``.
        fraction (float, optional): Fraction of the smaller span. Defaults to 0.02.

    Returns:
        float: Buffer distance in coordinate units.
    """
    xmin, xmax, ymin, ymax = extent
    return min(xmax - xmin, ymax - ymin) * fraction


def choose_tick_step(span, target=5):
    """Pick a nice tick step for map axes.

    Args:
        span (float): Axis span (e.g. longitude or latitude range).
        target (int, optional): Approximate number of ticks. Defaults to 5.

    Returns:
        float: Chosen step from a fixed candidate set.
    """
    steps = np.array([0.25, 0.5, 1, 2, 2.5, 5, 10, 15, 20])
    ideal = max(span / max(target, 1), 1e-9)
    return float(steps[np.argmin(np.abs(steps - ideal))])


def extent_grid_ticks(lo, hi, step):
    """Regular tick positions within ``[lo, hi]``.

    Args:
        lo (float): Lower bound.
        hi (float): Upper bound.
        step (float): Tick spacing.

    Returns:
        np.ndarray: Tick positions (may be empty).
    """
    if step <= 0 or not (np.isfinite(lo) and np.isfinite(hi)):
        return np.array([])
    first = np.ceil(lo / step - 1e-12) * step
    last = np.floor(hi / step + 1e-12) * step
    if first > last:
        return np.array([])
    n = int(np.round((last - first) / step)) + 1
    return first + np.arange(n, dtype=float) * step


def clean_geoms(sub):
    """Drop empty/invalid geometries and repair with ``make_valid`` where possible.

    Args:
        sub (gpd.GeoDataFrame): Input GeoDataFrame.

    Returns:
        gpd.GeoDataFrame: Filtered copy with valid finite-bounds geometries.
    """
    sub = sub.copy()
    sub = sub[sub.geometry.notna() & ~sub.geometry.is_empty]
    sub["geometry"] = sub.geometry.map(make_valid)
    sub = sub[sub.geometry.notna() & ~sub.geometry.is_empty]
    sub = sub[sub.geometry.is_valid]
    b = sub.geometry.bounds
    finite = np.isfinite(b[["minx", "miny", "maxx", "maxy"]]).all(axis=1)
    sub = sub.loc[finite]
    return sub


def simplify_extent(gdf, extent, tolerance):
    """Crop to extent and simplify polygon geometries.

    Args:
        gdf (gpd.GeoDataFrame): Input GeoDataFrame (must support ``.cx``).
        extent (tuple): ``(xmin, xmax, ymin, ymax)`` in the same CRS as ``gdf``.
        tolerance (float): Douglas-Peucker tolerance for ``simplify``.

    Returns:
        gpd.GeoDataFrame: Cropped and simplified geometries.
    """
    xmin, xmax, ymin, ymax = extent
    sub = gdf.cx[xmin:xmax, ymin:ymax].copy()
    sub = clean_geoms(sub)
    sub["geometry"] = sub.geometry.simplify(tolerance, preserve_topology=True)
    sub = sub[sub.geometry.notna() & ~sub.geometry.is_empty]
    return sub


def get_simplified(
    name,
    gdf,
    extent,
    cache_dir: Path,
    tolerance=None,
    buffer=None,
    force_recompute=False,
):
    """Load or build simplified (and optionally buffered) polygons for one panel extent.

    Caches results as pickle files under ``cache_dir``.

    Args:
        name (str): Panel label used in cache subdirectory names.
        gdf (gpd.GeoDataFrame): Full polygon layer to crop and simplify.
        extent (tuple): ``(xmin, xmax, ymin, ymax)``.
        cache_dir (Path): Root directory for panel caches.
        tolerance (float, optional): Simplification tolerance; if ``None``, uses ``auto_simplify_tolerance``.
        buffer (float, optional): Optional buffer in degrees; if ``None``, no buffer.
        force_recompute (bool, optional): If ``True``, ignore existing cache. Defaults to ``False``.

    Returns:
        gpd.GeoDataFrame: Simplified (and optionally buffered/clipped) polygons for plotting.
    """
    tol = tolerance if tolerance is not None else auto_simplify_tolerance(extent)
    buf = buffer if buffer is not None else 0.0
    simpl_dir = cache_dir / f"{name.replace(' ', '_')}_simp_tol_{tol:.6f}"
    simpl_dir.mkdir(parents=True, exist_ok=True)
    simplified_cache_path = simpl_dir / "simplified.pkl"
    buf_dir = simpl_dir / "buffered"
    buf_dir.mkdir(parents=True, exist_ok=True)
    buffered_cache_path = buf_dir / f"buf_{buf:.6f}.pkl"
    cache_path = simplified_cache_path if buf == 0.0 else buffered_cache_path

    if cache_path.exists() and not force_recompute:
        with open(cache_path, "rb") as f:
            return pickle.load(f)

    if simplified_cache_path.exists():
        with open(simplified_cache_path, "rb") as f:
            simp = pickle.load(f)
    else:
        try:
            simp = simplify_extent(gdf, extent, tol)
        except Exception:
            xmin, xmax, ymin, ymax = extent
            sub = gdf.cx[xmin:xmax, ymin:ymax].copy()
            try:
                sub["geometry"] = sub.geometry.simplify(tol, preserve_topology=False)
            except Exception:
                sub = sub[sub.geometry.notna() & ~sub.geometry.is_empty]
                simp = sub
            else:
                sub = sub[sub.geometry.notna() & ~sub.geometry.is_empty]
                if hasattr(sub.geometry, "is_valid"):
                    sub = sub[sub.geometry.is_valid]
                simp = sub
        with open(simplified_cache_path, "wb") as f:
            pickle.dump(simp, f)

    if buf > 0.0:
        simp = simp.copy()
        simp["geometry"] = simp.geometry.buffer(buf)
        clip_box = box(extent[0], extent[2], extent[1], extent[3])
        simp["geometry"] = simp.geometry.intersection(clip_box)
        simp = simp[simp.geometry.notna() & ~simp.geometry.is_empty]
        if hasattr(simp.geometry, "is_valid"):
            simp = simp[simp.geometry.is_valid]
        with open(buffered_cache_path, "wb") as f:
            pickle.dump(simp, f)
    return simp


def style_map_axes(
    ax,
    extent,
    *,
    show_bottom=False,
    show_left=False,
    show_top=False,
    show_right=False,
    show_labels=False,
    target_ticks=4,
):
    """Apply grid, ticks, and optional degree labels on a cartopy PlateCarree axis.

    Args:
        ax: Cartopy ``GeoAxes`` with PlateCarree projection.
        extent (tuple): ``(xmin, xmax, ymin, ymax)``.
        show_bottom (bool, optional): Show x tick labels on bottom. Defaults to ``False``.
        show_left (bool, optional): Show y tick labels on left. Defaults to ``False``.
        show_top (bool, optional): Show x tick labels on top. Defaults to ``False``.
        show_right (bool, optional): Show y tick labels on right. Defaults to ``False``.
        show_labels (bool, optional): If ``True``, configure tick formatters and grid. Defaults to ``False``.
        target_ticks (int, optional): Target tick count per axis. Defaults to 4.

    Returns:
        None
    """
    xmin, xmax, ymin, ymax = extent
    xstep = choose_tick_step(xmax - xmin, target=target_ticks)
    ystep = choose_tick_step(ymax - ymin, target=target_ticks)
    xticks = extent_grid_ticks(xmin, xmax, xstep)
    yticks = extent_grid_ticks(ymin, ymax, ystep)
    if show_labels:
        ax.set_xticks(xticks, crs=ccrs.PlateCarree())
        ax.set_yticks(yticks, crs=ccrs.PlateCarree())
        ax.xaxis.set_major_formatter(
            LongitudeFormatter(
                number_format=".0f", degree_symbol="°", dateline_direction_label=False
            )
        )
        ax.yaxis.set_major_formatter(
            LatitudeFormatter(number_format=".0f", degree_symbol="°")
        )
        ax.tick_params(
            axis="x",
            labelsize=8,
            pad=1.5,
            labelbottom=show_bottom,
            labeltop=show_top,
            bottom=False,
            top=True,
        )
        ax.tick_params(
            axis="y",
            labelsize=8,
            pad=1.5,
            labelleft=show_left,
            labelright=show_right,
            left=True,
            right=False,
        )
        ax.grid(True, color="#f1f1f1", alpha=1, zorder=-100)
    else:
        ax.set_xticks([])
        ax.set_yticks([])
    ax.set_extent(extent, crs=ccrs.PlateCarree())


def add_panel_label(ax, text, fontsize=12, _corner_pad_points=(-10, -10)):
    """Draw a small boxed title in the upper-left of the axes.

    Args:
        ax: Matplotlib axes.
        text (str): Text to display.
        fontsize (int, optional): Font size. Defaults to 12.

    Returns:
        None
    """
    ax.annotate(
        text,
        xy=(1.0, 1.0),
        xycoords="axes fraction",
        xytext=_corner_pad_points,
        textcoords="offset points",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=fontsize,
        bbox=dict(facecolor="white", edgecolor="0.3", boxstyle="square,pad=0.2"),
        zorder=100,
    )


def draw_geom_panel(
    ax: plt.axes.Axes,
    name: str,
    extent: tuple[float, float, float, float],
    title: str,
    *,
    simplified: dict,
    pt_gdf: gpd.GeoDataFrame,
    poly_val_col: str,
    pt_val_col: str,
    cmap: plt.colors.Colormap,
    norm: plt.colors.Normalize,
    land_colour: str,
    border_colour: str,
    river_colour: str,
    rivers_gdf: gpd.GeoDataFrame | None = None,
    show_bottom: bool = False,
    show_left: bool = False,
    show_top: bool = False,
    show_right: bool = False,
    show_labels: bool = False,
):
    """Draw one map panel: basemap, optional river overlay, choropleth polygons, and points.

    Args:
        ax: Cartopy ``GeoAxes`` (PlateCarree).
        name (str): Key into ``simplified`` for this panel's polygon GeoDataFrame.
        extent (tuple): ``(xmin, xmax, ymin, ymax)`` in degrees.
        title (str): Panel label text (upper-left box).
        simplified (dict): Mapping from panel name to simplified polygon ``GeoDataFrame``.
        pt_gdf (gpd.GeoDataFrame): Point geometries; colored by ``pt_val_col`` when present.
        poly_val_col (str): Column name for polygon choropleth values.
        pt_val_col (str): Column name for point scatter values; if missing, points plot light grey.
        cmap: Matplotlib colormap name or ``Colormap`` instance.
        norm: Normalization for the color scale.
        land_colour (str): Fill color for land.
        border_colour (str): Color for country borders.
        river_colour (str): Color for river features.
        rivers_gdf (gpd.GeoDataFrame, optional): Optional river lines (e.g. Somme highlight). Defaults to ``None``.
        show_bottom (bool, optional): Pass-through to ``style_map_axes``. Defaults to ``False``.
        show_left (bool, optional): Pass-through to ``style_map_axes``. Defaults to ``False``.
        show_top (bool, optional): Pass-through to ``style_map_axes``. Defaults to ``False``.
        show_right (bool, optional): Pass-through to ``style_map_axes``. Defaults to ``False``.
        show_labels (bool, optional): Pass-through to ``style_map_axes``. Defaults to ``False``.

    Returns:
        None
    """
    if rivers_gdf is not None:
        rivers_gdf = rivers_gdf.to_crs(ccrs.PlateCarree())

    xmin, xmax, ymin, ymax = extent
    ax.set_extent(extent, crs=ccrs.PlateCarree())
    ax.coastlines(linewidth=0.3, color="0.35")
    ax.add_feature(cfeature.LAND, color=land_colour, alpha=1, zorder=0)
    ax.add_feature(cfeature.BORDERS, color=border_colour, alpha=1, zorder=0, lw=0.5)
    ax.add_feature(
        cfeature.RIVERS.with_scale("10m"), color=river_colour, alpha=0.1, zorder=0
    )
    if rivers_gdf is not None:
        somme_rows = rivers_gdf[rivers_gdf.name.str.contains("Somme", na=False)]
        somme_rows.plot(
            ax=ax, transform=ccrs.PlateCarree(), color=river_colour, alpha=0.1, zorder=0
        )

    if name not in simplified:
        raise KeyError(f"Panel '{name}' not found in simplified geometry dictionary.")

    poly_sub = simplified[name].copy().cx[xmin:xmax, ymin:ymax]
    pts_sub = pt_gdf.cx[xmin:xmax, ymin:ymax]
    ms = dynamic_point_size(extent, n_points=len(pts_sub))

    if poly_val_col not in poly_sub.columns:
        raise KeyError(
            f"Polygon value column '{poly_val_col}' not found for panel '{name}'. "
            f"Available columns include: {', '.join(map(str, poly_sub.columns[:12]))}"
        )

    if len(pts_sub) > 0:
        xv = pts_sub.geometry.x.to_numpy()
        yv = pts_sub.geometry.y.to_numpy()
        if pt_val_col not in pts_sub.columns:
            ax.scatter(
                xv,
                yv,
                c="grey",
                s=ms,
                alpha=0.7,
                linewidths=0,
                transform=ccrs.PlateCarree(),
                zorder=1,
            )
        else:
            vv = pts_sub[pt_val_col].to_numpy(dtype=float)
            ok = np.isfinite(vv) & (vv > 0)
            if np.any(~ok):
                ax.scatter(
                    xv[~ok],
                    yv[~ok],
                    c="grey",
                    s=ms,
                    alpha=0.7,
                    linewidths=0,
                    transform=ccrs.PlateCarree(),
                    zorder=1,
                )
            if np.any(ok):
                ax.scatter(
                    xv[ok],
                    yv[ok],
                    c=vv[ok],
                    s=ms,
                    cmap=cmap,
                    norm=norm,
                    alpha=0.7,
                    linewidths=0,
                    transform=ccrs.PlateCarree(),
                    zorder=1,
                )

    if len(poly_sub) > 0:
        poly_sub.plot(
            ax=ax,
            transform=ccrs.PlateCarree(),
            column=poly_val_col,
            cmap=cmap,
            norm=norm,
            linewidth=0.25,
            alpha=0.9,
            zorder=10,
            missing_kwds={"color": "grey"},
        )

    style_map_axes(
        ax,
        extent,
        show_bottom=show_bottom,
        show_left=show_left,
        show_top=show_top,
        show_right=show_right,
        show_labels=show_labels,
        target_ticks=4,
    )
    add_panel_label(ax, title)
    for spine in ax.spines.values():
        spine.set_linewidth(1.0)
        spine.set_edgecolor("0.15")
    ax.set_extent(extent, crs=ccrs.PlateCarree())
    ax.set_aspect("auto", adjustable="box")


def to_wgs84_clean(gdf):
    """Ensure WGS84 (EPSG:4326) and drop empty geometries.

    Args:
        gdf (gpd.GeoDataFrame): Input layer.

    Returns:
        gpd.GeoDataFrame: Copy in EPSG:4326 with only non-empty geometries.
    """
    out = gdf.copy()
    if out.crs is None:
        out = out.set_crs("EPSG:4326")
    if str(out.crs).lower() not in {"epsg:4326", "4326"}:
        out = out.to_crs(4326)
    return out[out.geometry.notna() & ~out.geometry.is_empty].copy()
