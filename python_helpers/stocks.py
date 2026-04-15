from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def convert_carbon_density_to_carbon_stock(
    carbon_density: np.ndarray, depth: float, area: np.ndarray
) -> np.ndarray:
    """Convert carbon density from gC/cm3 to total stock in GgC."""
    return carbon_density * depth * 1e4 * area / 1e9


def convert_carbon_density_to_carbon_areal_stock(
    carbon_density: np.ndarray, depth: float
) -> np.ndarray:
    """Convert carbon density from gC/cm3 to areal stock in kgC/ha."""
    return carbon_density * depth * 1e4 * 1e4 * 1e-3


def safe_weighted_average(values, weights):
    values = np.asarray(values)
    weights = np.asarray(weights)
    mask = ~np.isnan(values)
    if np.sum(mask) == 0 or np.sum(weights[mask]) == 0:
        return np.nan
    return np.average(values[mask], weights=weights[mask])


def weighted_mean_standard_error(se, weights, values_for_mask):
    """SE of area-weighted mean: sqrt(sum(w_i^2 * se_i^2)) / sum(w_i); mask matches ``safe_weighted_average``."""
    se = np.asarray(se, dtype=float)
    weights = np.asarray(weights, dtype=float)
    values_for_mask = np.asarray(values_for_mask, dtype=float)
    mask = ~np.isnan(values_for_mask)
    if np.sum(mask) == 0 or np.sum(weights[mask]) == 0:
        return np.nan
    se_m = se[mask]
    w_m = weights[mask]
    if not np.all(np.isfinite(se_m)) or not np.all(np.isfinite(w_m)):
        return np.nan
    if np.any(w_m <= 0):
        return np.nan
    return float(np.sqrt(np.sum((w_m**2) * (se_m**2))) / np.sum(w_m))


def plot_comparison_bars(df: pd.DataFrame, dpi: int = 300):
    GOMIS_COLOR = "#0C0787"
    STUDY_COLOR = "#FBA237"
    BAR_WIDTH = 0.38
    ALPHA = 0.85

    fig, axes = plt.subplots(3, 1, figsize=(12, 8), dpi=dpi)
    fig.subplots_adjust(hspace=5)

    panels = [
        {
            "ax": axes[0],
            "col_gomis": "Seagrass extent (ha)",
            "col_study": "area_ha",
            "col_study_se": None,
            "ylabel": "Seagrass extent (ha)",
            "title": "National seagrass extent",
        },
        {
            "ax": axes[1],
            "col_gomis": "Biomass carbon stock (kg C / ha)",
            "col_study": "stock_kg_ha",
            "col_study_se": "stock_kg_ha_se",
            "ylabel": "Stock density (kg C / ha)",
            "title": "Mean sediment carbon stock density",
        },
        {
            "ax": axes[2],
            "col_gomis": "Biomass carbon stock (Gg)",
            "col_study": "stock_Gg",
            "col_study_se": "stock_Gg_se",
            "ylabel": "Total carbon stock (Gg C)",
            "title": "Total sediment carbon stock",
        },
    ]

    df = df.sort_values("area_ha", ascending=False, na_position="last")
    countries = df["Country"].tolist()
    x = np.arange(len(countries))

    for idx, p in enumerate(panels):
        ax = p["ax"]
        cg, cs, cse = p["col_gomis"], p["col_study"], p["col_study_se"]
        panel_df = df[["Country", cg, cs] + ([cse] if cse is not None else [])].copy()

        gomis_vals = pd.to_numeric(panel_df[cg], errors="coerce").to_numpy(dtype=float)
        study_vals = pd.to_numeric(panel_df[cs], errors="coerce").to_numpy(dtype=float)
        study_vals_se = None
        if cse is not None:
            study_vals_se = pd.to_numeric(panel_df[cse], errors="coerce").to_numpy(
                dtype=float
            )

        ax.bar(
            x - BAR_WIDTH / 2,
            np.where(np.isfinite(gomis_vals), gomis_vals, 0),
            width=BAR_WIDTH,
            color=GOMIS_COLOR,
            alpha=ALPHA,
            label="Gomis et al. 2025 (organic biomass carbon stock)",
            zorder=3,
        )
        ax.bar(
            x + BAR_WIDTH / 2,
            np.where(np.isfinite(study_vals), study_vals, 0),
            yerr=study_vals_se,
            width=BAR_WIDTH,
            color=STUDY_COLOR,
            alpha=ALPHA,
            label="This study (sediment carbon stock)",
            zorder=3,
        )

        # Mark missing values with crosses so one-sided availability is explicit.
        missing_gomis = ~np.isfinite(gomis_vals)
        missing_study = ~np.isfinite(study_vals)

        ax.set_xticks(x)
        country_labels = [str(country).replace(" ", "\n") for country in countries]
        ax.set_xticklabels(country_labels, ha="center", fontsize=10, rotation=45)
        ax.set_ylabel(p["ylabel"], fontsize=11)
        ax.set_title(
            p["title"], fontsize=13, fontweight="bold", pad=10, x=0.0, loc="left"
        )
        ax.yaxis.grid(True, linestyle="--", alpha=0.5, zorder=0)
        ax.set_axisbelow(True)
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_yscale("log")

        # Place crosses just above lower log bound so they are visible.
        y_min, y_max = ax.get_ylim()
        if y_min > 0 and y_max > y_min:
            log_min = np.log10(y_min)
            log_max = np.log10(y_max)
            y_cross = 10 ** (log_min + 0.08 * (log_max - log_min))
        else:
            y_cross = y_min

        if np.any(missing_gomis):
            ax.plot(
                x[missing_gomis] - BAR_WIDTH / 2,
                np.full(np.sum(missing_gomis), y_cross),
                linestyle="None",
                marker="x",
                color=GOMIS_COLOR,
                markersize=6,
                markeredgewidth=1.3,
                zorder=5,
            )
        if np.any(missing_study):
            ax.plot(
                x[missing_study] + BAR_WIDTH / 2,
                np.full(np.sum(missing_study), y_cross),
                linestyle="None",
                marker="x",
                color=STUDY_COLOR,
                markersize=6,
                markeredgewidth=1.3,
                zorder=5,
            )

        if idx < 2:
            ax.set_xticklabels([])

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        fontsize=10,
        framealpha=0.92,
        ncol=2,
        loc="lower center",
        bbox_to_anchor=(0.5, 0),
    )
    plt.tight_layout(rect=[0, 0.05, 1, 1.0])
    return fig, axes


def _territory_total_stock_se(
    df: pd.DataFrame,
    territory_col: str,
    se_col: str,
) -> pd.Series:
    """Standard error of summed stock per territory (sqrt of summed variances)."""
    if se_col not in df.columns:
        return pd.Series(dtype=float)
    se = pd.to_numeric(df[se_col], errors="coerce").to_numpy(dtype=float)
    g2 = df[[territory_col]].copy()
    g2["_var"] = np.nan_to_num(np.square(se), nan=0.0)
    summed = g2.groupby(territory_col, sort=False)["_var"].sum()
    return np.sqrt(summed)


def _species_ordered_stack_pivot(
    g: pd.DataFrame,
    species_color_map: dict[str, str],
    *,
    territory_col: str = "TERRITORY1",
    species_col: str = "seagrass_species",
    value_col: str = "stock_Gg",
) -> tuple[pd.DataFrame, list[str]]:
    """Territory × species pivot with species columns ordered like ``species_color_map``."""
    stacked = (
        g.groupby([territory_col, species_col])[value_col].sum().unstack(fill_value=0)
    )
    canon = [s for s in species_color_map if s in stacked.columns]
    species_order = canon + [s for s in stacked.columns if s not in canon]
    stacked = stacked[species_order]
    colors = [species_color_map.get(s, "#777777") for s in stacked.columns]
    return stacked, colors


def _prettify_territory_xlabels(ax: plt.Axes) -> None:
    ax.set_xticklabels(
        [t.get_text().replace(" ", "\n") for t in ax.get_xticklabels()],
        ha="center",
    )
    ax.tick_params(axis="x", rotation=45)


def _plot_stacked_species_bars_on_ax(
    ax: plt.Axes,
    data: pd.DataFrame,
    colors: list[str],
    *,
    y_label: str,
    territory_se: pd.Series | None,
    errorbar_capsize: float,
    errorbar_color: str,
) -> None:
    """Draw stacked species bars and optional per-territory total error bars."""
    data.plot(
        kind="bar",
        stacked=True,
        ax=ax,
        color=colors,
        edgecolor="none",
        legend=False,
    )
    ax.set_ylabel(y_label)
    ax.set_axisbelow(True)
    ax.grid(True, which="major", axis="y", alpha=0.5)
    ax.grid(True, which="minor", axis="y", alpha=0.1)
    _prettify_territory_xlabels(ax)
    n_bars = len(data)
    if territory_se is None or territory_se.empty or n_bars == 0:
        return
    x_pos = np.arange(n_bars, dtype=float)
    totals = data.sum(axis=1).to_numpy(dtype=float)
    errs = territory_se.reindex(data.index).to_numpy(dtype=float)
    valid = np.isfinite(totals) & np.isfinite(errs)
    if np.any(valid):
        ax.errorbar(
            x_pos[valid],
            totals[valid],
            yerr=errs[valid],
            fmt="none",
            ecolor=errorbar_color,
            elinewidth=1.0,
            capsize=errorbar_capsize,
            zorder=30,
        )


def plot_territory_species_stacked_all(
    g: pd.DataFrame,
    *,
    species_color_map: dict[str, str],
    territory_col: str = "TERRITORY1",
    species_col: str = "seagrass_species",
    value_col: str = "stock_Gg",
    se_col: str = "stock_Gg_se",
    show_total_errorbars: bool = False,
    dpi: int = 300,
    figsize: tuple[float, float] = (12.0, 5.0),
    ylim: tuple[float, float] | None = None,
    xlabel: str = "National Territory",
    ylabel: str = "Carbon Stock (GgC)",
    legend_bbox: tuple[float, float] = (0.5, 1.1),
    errorbar_capsize: float = 3.0,
    errorbar_color: str = "0.15",
) -> tuple[plt.Figure, plt.Axes, pd.DataFrame]:
    """Single stacked bar chart: all territories by species, ordered by total stock.

    Args:
        g: Dataframe with territory, species, and stock (and optional SE) columns.
        species_color_map: Species-to-color mapping (column order follows this dict).
        territory_col: Territory column name.
        species_col: Species column name.
        value_col: Stock column (GgC).
        se_col: Per-row stock standard error; used only if ``show_total_errorbars``.
        show_total_errorbars: If ``True``, draw one SE bar on each stacked total
            (``sqrt(sum se^2)`` per territory).
        dpi: Figure DPI.
        figsize: ``(width, height)`` inches.
        ylim: Fixed y-axis limits, or ``None`` for matplotlib default.
        xlabel: X axis label.
        ylabel: Y axis label.
        legend_bbox: ``(x, y)`` in axes fraction terms for ``bbox_to_anchor``.
        errorbar_capsize: Cap width when error bars are shown.
        errorbar_color: Error bar color when error bars are shown.

    Returns:
        ``(fig, ax, stacked_table)`` — territory × species table actually plotted
        (territories sorted by descending total stock).
    """
    stacked, colors = _species_ordered_stack_pivot(
        g,
        species_color_map,
        territory_col=territory_col,
        species_col=species_col,
        value_col=value_col,
    )
    stacked = stacked.loc[stacked.sum(axis=1).sort_values(ascending=False).index]
    territory_se: pd.Series | None = None
    if show_total_errorbars:
        territory_se = _territory_total_stock_se(g, territory_col, se_col)

    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    _plot_stacked_species_bars_on_ax(
        ax,
        stacked,
        colors,
        y_label=ylabel,
        territory_se=territory_se,
        errorbar_capsize=errorbar_capsize,
        errorbar_color=errorbar_color,
    )
    ax.set_xlabel(xlabel)
    handles, labels = ax.get_legend_handles_labels()
    legend = ax.legend(
        handles,
        labels,
        title="Seagrass Species",
        bbox_to_anchor=legend_bbox,
        ncols=len(stacked.columns),
        loc="center",
        borderaxespad=0.0,
    )
    if legend.get_title() is not None:
        legend.get_title().set_fontweight("bold")
    if ylim is not None:
        ax.set_ylim(*ylim)
    plt.tight_layout()
    return fig, ax, stacked


def plot_territory_species_stacked_split(
    g: pd.DataFrame,
    *,
    species_color_map: dict[str, str],
    territory_col: str = "TERRITORY1",
    species_col: str = "seagrass_species",
    value_col: str = "stock_Gg",
    se_col: str = "stock_Gg_se",
    top_n: int = 5,
    dpi: int = 300,
    figsize_height: float = 5.0,
    ylim_left: tuple[float, float] | None = (0, 26e4),
    ylim_right: tuple[float, float] | None = (0, 31e2),
    left_panel_title: str = "A) Top 5 National Territories",
    right_panel_title: str = "B) Remaining National Territories",
    ylabel: str = "Carbon Stock (GgC)",
    supxlabel: str = "National Territory",
    errorbar_capsize: float = 3.0,
    errorbar_color: str = "0.15",
) -> tuple[plt.Figure, tuple[plt.Axes, plt.Axes]]:
    """Stacked species contributions by territory in two panels, with total SE error bars.

    Uncertainty for each territory total is ``sqrt(sum se_i^2)`` over rows (independence).

    Args:
        g: Dataframe with territory, species, stock, and optional per-row standard error.
        species_color_map: Species label to face color for stack segments.
        territory_col: Column name for national territory.
        species_col: Column name for seagrass species.
        value_col: Column name for carbon stock (GgC).
        se_col: Column name for per-row stock standard error; if missing, no error bars.
        top_n: Number of largest territories (by total stock) on the left panel.
        dpi: Figure resolution.
        figsize_height: Figure height in inches; width follows panel count ratio.
        ylim_left: Y limits for the left (top territories) axis; ``None`` to autoscale.
        ylim_right: Y limits for the right axis; ``None`` to autoscale.
        left_panel_title: Annotation on the top-right of the left axes.
        right_panel_title: Annotation on the top-right of the right axes.
        ylabel: Left axis y-label (right panel y-label is left blank).
        supxlabel: Figure-level x label below both panels.
        errorbar_capsize: Cap width for total-stock error bars.
        errorbar_color: Error bar line color.

    Returns:
        ``(fig, (ax_left, ax_right))``
    """
    stacked, colors = _species_ordered_stack_pivot(
        g,
        species_color_map,
        territory_col=territory_col,
        species_col=species_col,
        value_col=value_col,
    )

    territory_totals = stacked.sum(axis=1).sort_values(ascending=False)
    territory_se = _territory_total_stock_se(g, territory_col, se_col)
    top = territory_totals.index[:top_n]
    rest = territory_totals.index[top_n:]
    stacked_top, stacked_rest = stacked.loc[top], stacked.loc[rest]

    n_top = len(stacked_top)
    n_rest = len(stacked_rest)
    width_ratios = [1, max(1, n_rest / max(n_top, 1))]

    fig, (ax_l, ax_r) = plt.subplots(
        1,
        2,
        figsize=(14.0, figsize_height),
        dpi=dpi,
        gridspec_kw={"width_ratios": width_ratios},
    )

    _plot_stacked_species_bars_on_ax(
        ax_l,
        stacked_top,
        colors,
        y_label=ylabel,
        territory_se=territory_se,
        errorbar_capsize=errorbar_capsize,
        errorbar_color=errorbar_color,
    )
    ax_l.set_xlabel("")
    # Inset from each panel's top-right in points (not axes fraction) so padding matches
    # when subplot widths differ.
    _corner_pad_pts = (-10.0, -10.0)
    ax_l.annotate(
        left_panel_title,
        xy=(1.0, 1.0),
        xycoords="axes fraction",
        xytext=_corner_pad_pts,
        textcoords="offset points",
        ha="right",
        va="top",
        fontsize=11,
        bbox=dict(
            boxstyle="round,pad=0.3", facecolor="white", edgecolor="#cccccc", alpha=0.85
        ),
    )

    _plot_stacked_species_bars_on_ax(
        ax_r,
        stacked_rest,
        colors,
        y_label="",
        territory_se=territory_se,
        errorbar_capsize=errorbar_capsize,
        errorbar_color=errorbar_color,
    )
    ax_r.set_xlabel("")
    ax_r.annotate(
        right_panel_title,
        xy=(1.0, 1.0),
        xycoords="axes fraction",
        xytext=_corner_pad_pts,
        textcoords="offset points",
        ha="right",
        va="top",
        fontsize=11,
        bbox=dict(
            boxstyle="round,pad=0.3", facecolor="white", edgecolor="#cccccc", alpha=0.85
        ),
    )

    fig.supxlabel(supxlabel)
    handles, labels = ax_l.get_legend_handles_labels()
    legend = fig.legend(
        handles,
        labels,
        title="Seagrass Species",
        bbox_to_anchor=(0.5, 1.0),
        ncol=len(stacked.columns),
        loc="center",
        borderaxespad=0.0,
    )
    if legend.get_title():
        legend.get_title().set_fontweight("bold")

    if ylim_left is not None:
        ax_l.set_ylim(*ylim_left)
    if ylim_right is not None:
        ax_r.set_ylim(*ylim_right)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    return fig, (ax_l, ax_r)
