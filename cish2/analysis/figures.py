"""Publication figures for a cis-heritability benchmarking campaign.

Each function makes one figure that answers one question, and returns the
paths it wrote (PNG for review at 400 dpi, PDF as the vector master with real
embedded fonts).  Figures adapt to the campaign: an axis that does not vary in
a given run is simply not drawn.
"""

from __future__ import annotations

from typing import Dict, List, Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from .load import per_cell_timings
from .metrics import summarize_by, win_rate
from .style import (
    BLUE_RAMP,
    DIVERGING,
    ESTIMATOR_ORDER,
    INK,
    INK_MUTED,
    INK_SECONDARY,
    REFERENCE,
    SEQUENTIAL,
    SERIES,
    apply_style,
    finish,
    ordinal_colors,
    panel_label,
    save,
)

AXIS_LABELS = {
    "ncausal": "causal SNPs",
    "h2_true": "true cis-$h^2$",
    "n": "sample size $n$",
    "effect_dist": "effect-size prior",
    "maf_bin": "MAF bin",
    "ld_arch": "LD architecture",
}


def present(frame: pd.DataFrame) -> List[str]:
    """Estimators in this frame, in the fixed palette order.

    Figures iterate this rather than ESTIMATOR_ORDER so a campaign run before
    HEELS existed still plots, with two series instead of three.
    """
    if "estimator" in frame.columns:
        have = set(map(str, frame["estimator"].dropna().unique()))
    else:
        have = {e for e in ESTIMATOR_ORDER if f"abs_err_{e}" in frame.columns}
    return [e for e in ESTIMATOR_ORDER if e in have]


def _estimator_legend(ax, loc: str = "best", extra: Optional[List] = None,
                      estimators: Optional[List[str]] = None) -> None:
    handles = [
        Line2D([], [], color=SERIES[e], lw=2.2, marker="o", ms=5, label=e)
        for e in (estimators or ESTIMATOR_ORDER)
    ]
    ax.legend(handles=(extra or []) + handles, loc=loc, handlelength=1.6)


def _estimator_patch_legend(ax, loc: str = "upper right",
                            estimators: Optional[List[str]] = None) -> None:
    """Legend for bar panels - swatches, not line+marker handles."""
    ests = estimators or ESTIMATOR_ORDER
    ax.legend(handles=[Patch(facecolor=SERIES[e], edgecolor="white", label=e)
                       for e in ests],
              loc=loc, ncol=len(ests), handlelength=1.1, columnspacing=1.2)


def _best_estimator_panel(ax, fig, wins: pd.DataFrame, rows, cols,
                          row_key: str, col_key: str, estimators: List[str],
                          row_label: str, col_label: str) -> None:
    """Which estimator wins each design cell, as a categorical grid.

    A win *rate* heatmap only works for two methods; with three the readable
    encoding is the identity of the winner, coloured by the same fixed hue each
    estimator carries everywhere else, annotated with its winning share.
    """
    from matplotlib.colors import BoundaryNorm, ListedColormap

    share = np.full((len(rows), len(cols), len(estimators)), np.nan)
    for i, r in enumerate(rows):
        for j, c in enumerate(cols):
            hit = wins[(wins[row_key] == r) & (wins[col_key] == c)]
            if hit.empty:
                continue
            for k, est in enumerate(estimators):
                share[i, j, k] = float(hit.iloc[0].get(f"win_{est}", np.nan))

    best = np.full((len(rows), len(cols)), np.nan)
    with np.errstate(invalid="ignore"):
        ok = np.isfinite(share).all(axis=2)
        best[ok] = np.argmax(share[ok], axis=1)

    cmap = ListedColormap([SERIES[e] for e in estimators])
    norm = BoundaryNorm(np.arange(-0.5, len(estimators)), cmap.N)
    im = ax.imshow(best, cmap=cmap, norm=norm, aspect="auto", origin="lower")
    for i in range(len(rows)):
        for j in range(len(cols)):
            if not np.isfinite(best[i, j]):
                continue
            k = int(best[i, j])
            r, g, b, _ = cmap(k)
            lum = 0.2126 * r + 0.7152 * g + 0.0722 * b
            ax.text(j, i, f"{share[i, j, k]:.0%}", ha="center", va="center",
                    fontsize=7.5, color="white" if lum < 0.5 else INK)
    _heat_axes(ax, rows, cols)
    finish(ax, xlabel=col_label, ylabel=row_label, title="Which estimator wins")
    ax.legend(handles=[Patch(facecolor=SERIES[e], edgecolor="white", label=e)
                       for e in estimators],
              loc="upper center", bbox_to_anchor=(0.5, -0.24),
              ncol=len(estimators), handlelength=1.1, columnspacing=1.0,
              fontsize=8)


def _headroom(ax, frac: float = 0.28) -> None:
    """Grow the axis upward so a legend never lands on the data."""
    lo, hi = ax.get_ylim()
    ax.set_ylim(lo, hi + frac * (hi - lo))


def _facet_grid(n_panels: int, width: float = 2.6, height: float = 2.5, max_cols: int = 5):
    ncols = min(max_cols, n_panels)
    nrows = int(np.ceil(n_panels / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(width * ncols + 1.0, height * nrows + 0.6),
                             squeeze=False, sharex=True, sharey=True)
    return fig, axes.ravel(), nrows, ncols


# ---------------------------------------------------------------------------
# Figure 1 - accuracy across the heritability grid
# ---------------------------------------------------------------------------
def fig_accuracy_grid(long_df: pd.DataFrame, out_dir: str, facet: str = "ncausal") -> List[str]:
    """Mean estimated h2 against true h2, one panel per architecture level."""
    apply_style()
    levels = sorted(long_df[facet].dropna().unique())
    stats = summarize_by(long_df, [facet, "h2_true"])

    ests = present(long_df)
    fig, axes, nrows, ncols = _facet_grid(len(levels))
    for k, (ax, level) in enumerate(zip(axes, levels)):
        lims = [0, max(long_df["h2_true"].max() * 1.08, 0.9)]
        ax.plot(lims, lims, ls=(0, (4, 3)), lw=1.1, color=REFERENCE, zorder=1)
        for est in ests:
            part = stats[(stats["estimator"] == est) & (stats[facet] == level)].sort_values("h2_true")
            if part.empty:
                continue
            ax.errorbar(part["h2_true"], part["mean_h2_hat"],
                        yerr=1.96 * part["bias_mcse"], color=SERIES[est],
                        marker="o", ms=5, lw=2.0, elinewidth=1.1, capsize=2.4,
                        markeredgecolor="white", markeredgewidth=0.7, zorder=3)
        ax.set_xlim(lims)
        ax.set_ylim(lims)
        finish(ax, title=f"{AXIS_LABELS.get(facet, facet)} = {_fmt(level)}")

    for ax in axes[len(levels):]:
        ax.set_visible(False)
    for i, ax in enumerate(axes[:len(levels)]):
        if i % ncols == 0:
            ax.set_ylabel("mean estimated $h^2$")
        if i >= len(levels) - ncols:
            ax.set_xlabel("true $h^2$")

    _estimator_legend(
        axes[0], loc="upper left", estimators=ests,
        extra=[Line2D([], [], color=REFERENCE, ls=(0, (4, 3)), lw=1.1, label="unbiased")],
    )
    fig.suptitle("Estimated vs. true cis-heritability across genetic architectures",
                 fontsize=11.5, fontweight="bold", color=INK, y=1.0)
    fig.text(0.5, -0.02,
             "Points are means over replicate loci; bars are 95% Monte-Carlo intervals.",
             ha="center", fontsize=8.5, color=INK_SECONDARY)
    fig.tight_layout()
    return save(fig, out_dir, "fig1_accuracy_grid")


# ---------------------------------------------------------------------------
# Figure 2 - bias as a function of polygenicity
# ---------------------------------------------------------------------------
def fig_bias_by_polygenicity(long_df: pd.DataFrame, out_dir: str) -> List[str]:
    """Bias vs. number of causal SNPs, one panel per true h2."""
    apply_style()
    levels = sorted(long_df["h2_true"].dropna().unique())
    stats = summarize_by(long_df, ["h2_true", "ncausal"])

    ests = present(long_df)
    fig, axes, nrows, ncols = _facet_grid(len(levels), width=2.5, height=2.3, max_cols=6)
    for ax, level in zip(axes, levels):
        ax.axhline(0, color=REFERENCE, lw=1.1, ls=(0, (4, 3)), zorder=1)
        for est in ests:
            part = stats[(stats["estimator"] == est) & (stats["h2_true"] == level)].sort_values("ncausal")
            if part.empty:
                continue
            ax.errorbar(part["ncausal"], part["bias"], yerr=1.96 * part["bias_mcse"],
                        color=SERIES[est], marker="o", ms=5, lw=2.0, elinewidth=1.1,
                        capsize=2.4, markeredgecolor="white", markeredgewidth=0.7, zorder=3)
        ax.set_xscale("log")
        ax.set_xticks(sorted(long_df["ncausal"].unique()))
        ax.get_xaxis().set_major_formatter(plt.ScalarFormatter())
        finish(ax, title=f"true $h^2$ = {_fmt(level)}")

    for ax in axes[len(levels):]:
        ax.set_visible(False)
    for i, ax in enumerate(axes[:len(levels)]):
        if i % ncols == 0:
            ax.set_ylabel("bias  ($\\hat{h}^2 - h^2$)")
        if i >= len(levels) - ncols:
            ax.set_xlabel("causal SNPs")
    _estimator_legend(axes[0], loc="lower left", estimators=ests)
    fig.suptitle("Bias grows with polygenicity for the fine-mapping estimator",
                 fontsize=11.5, fontweight="bold", color=INK, y=1.0)
    fig.tight_layout()
    return save(fig, out_dir, "fig2_bias_by_polygenicity")


# ---------------------------------------------------------------------------
# Figure 3 - RMSE landscape and where each estimator wins
# ---------------------------------------------------------------------------
def fig_rmse_landscape(long_df: pd.DataFrame, wide_df: pd.DataFrame, out_dir: str) -> List[str]:
    """Two RMSE heatmaps plus the paired win-rate over the same grid."""
    apply_style()
    stats = summarize_by(long_df, ["ncausal", "h2_true"])
    rows = sorted(long_df["ncausal"].unique())
    cols = sorted(long_df["h2_true"].unique())

    def grid_for(est: str) -> np.ndarray:
        part = stats[stats["estimator"] == est]
        return (part.pivot(index="ncausal", columns="h2_true", values="rmse")
                .reindex(index=rows, columns=cols).to_numpy(dtype=float))

    ests = present(long_df)
    grids = {e: grid_for(e) for e in ests}
    wins = win_rate(wide_df, ["ncausal", "h2_true"])

    fig, axes = plt.subplots(1, len(ests) + 1, figsize=(4.1 * (len(ests) + 1), 3.9))
    vmax = float(np.nanmax([g for g in grids.values()]))

    # every RMSE panel shares one scale, so they share one colorbar
    labels = {"GCTA": "GCTA / REML", "SuSiE": "SuSiE fine-mapping", "HEELS": "HEELS"}
    for k, est in enumerate(ests):
        ax = axes[k]
        im = ax.imshow(grids[est], cmap=SEQUENTIAL, vmin=0, vmax=vmax,
                       aspect="auto", origin="lower")
        _annotate_cells(ax, grids[est], im, fmt="{:.3f}")
        _heat_axes(ax, rows, cols)
        finish(ax, xlabel="true $h^2$", ylabel="causal SNPs",
               title=f"RMSE - {labels.get(est, est)}")
        panel_label(ax, chr(ord("A") + k))
    cb = fig.colorbar(im, ax=axes[len(ests) - 1], fraction=0.042, pad=0.03)
    cb.outline.set_visible(False)
    cb.set_label("RMSE", fontsize=8.5, color=INK_SECONDARY)

    _best_estimator_panel(axes[-1], fig, wins, rows, cols, "ncausal", "h2_true",
                          ests, "causal SNPs", "true $h^2$")
    panel_label(axes[-1], chr(ord("A") + len(ests)))

    fig.suptitle("Where each estimator is more accurate", fontsize=11.5,
                 fontweight="bold", color=INK, y=1.03)
    fig.tight_layout()
    return save(fig, out_dir, "fig3_rmse_landscape")


def _heat_axes(ax, rows, cols) -> None:
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels([_fmt(c) for c in cols])
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([_fmt(r) for r in rows])
    ax.grid(False)


def _annotate_cells(ax, grid, im, fmt: str) -> None:
    """Label every cell, picking ink or white from the cell's own luminance."""
    for i in range(grid.shape[0]):
        for j in range(grid.shape[1]):
            v = grid[i, j]
            if not np.isfinite(v):
                continue
            r, g, b, _ = im.cmap(im.norm(v))
            luminance = 0.2126 * r + 0.7152 * g + 0.0722 * b
            ax.text(j, i, fmt.format(v), ha="center", va="center", fontsize=7.5,
                    color="white" if luminance < 0.5 else INK)


# ---------------------------------------------------------------------------
# Figure 4 - uncertainty calibration
# ---------------------------------------------------------------------------
def fig_calibration(long_df: pd.DataFrame, out_dir: str) -> List[str]:
    """Are the reported standard errors honest, and do the intervals cover?"""
    apply_style()
    stats = summarize_by(long_df, ["ncausal", "h2_true"])
    se_mean = (long_df.groupby(["estimator", "ncausal", "h2_true"], observed=True)["h2_se"]
               .mean().reset_index().rename(columns={"h2_se": "mean_se"}))
    stats = stats.merge(se_mean, on=["estimator", "ncausal", "h2_true"], how="left")

    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.0))

    ests = present(long_df)
    ax = axes[0]
    hi = float(np.nanmax([stats["mean_se"].max(), stats["sd_h2_hat"].max()])) * 1.08
    ax.plot([0, hi], [0, hi], ls=(0, (4, 3)), lw=1.1, color=REFERENCE, zorder=1)
    for est in ests:
        part = stats[stats["estimator"] == est]
        ax.scatter(part["sd_h2_hat"], part["mean_se"], s=34, color=SERIES[est],
                   edgecolor="white", linewidth=0.7, alpha=0.9, zorder=3)
    ax.set_xlim(0, hi)
    ax.set_ylim(0, hi)
    finish(ax, xlabel="empirical SD across replicate loci",
           ylabel="mean reported SE", title="Reported vs. actual uncertainty")
    panel_label(ax, "A")
    _estimator_legend(ax, loc="upper left", estimators=ests,
                      extra=[Line2D([], [], color=REFERENCE, ls=(0, (4, 3)), lw=1.1,
                                    label="honest SE")])

    ax = axes[1]
    ax.axhline(0.95, color=REFERENCE, lw=1.1, ls=(0, (4, 3)), zorder=1)
    cov = summarize_by(long_df, ["h2_true"])
    for est in ests:
        part = cov[cov["estimator"] == est].sort_values("h2_true")
        ax.errorbar(part["h2_true"], part["coverage95"],
                    yerr=1.96 * part["coverage95_mcse"], color=SERIES[est],
                    marker="o", ms=5, lw=2.0, elinewidth=1.1, capsize=2.4,
                    markeredgecolor="white", markeredgewidth=0.7, zorder=3)
    ax.set_ylim(0, 1.02)
    ax.text(long_df["h2_true"].max(), 0.955, "nominal 95%", ha="right", va="bottom",
            fontsize=8, color=INK_SECONDARY)
    finish(ax, xlabel="true $h^2$", ylabel="coverage of the 95% interval",
           title="Interval coverage")
    panel_label(ax, "B")
    _estimator_legend(ax, loc="lower right", estimators=ests)

    fig.suptitle("Calibration of the reported standard errors", fontsize=11.5,
                 fontweight="bold", color=INK, y=1.02)
    fig.tight_layout()
    return save(fig, out_dir, "fig4_calibration")


# ---------------------------------------------------------------------------
# Figure 5 - head-to-head agreement
# ---------------------------------------------------------------------------
def fig_head_to_head(wide_df: pd.DataFrame, out_dir: str) -> List[str]:
    """Paired comparison on identical phenotypes, GCTA as the reference axis.

    GCTA is the individual-level benchmark, so the summary-statistic estimators
    are each plotted against it rather than against one another.
    """
    apply_style()
    ests = present(wide_df)
    others = [e for e in ests if e != "GCTA"]
    col = {"GCTA": "gcta_h2", "SuSiE": "susie_h2", "HEELS": "heels_h2"}

    fig, axes = plt.subplots(1, 3, figsize=(12.6, 4.0))

    # A: each summary-statistic estimator against the individual-level benchmark
    ax = axes[0]
    vals = pd.concat([wide_df[col[e]] for e in ests])
    lo = float(min(vals.quantile(0.001), 0.0))
    hi = float(max(vals.quantile(0.999), 1.0))
    ax.plot([lo, hi], [lo, hi], ls=(0, (4, 3)), lw=1.1, color=REFERENCE, zorder=1)
    for est in others:
        ax.scatter(wide_df["gcta_h2"], wide_df[col[est]], s=7, color=SERIES[est],
                   alpha=0.4, linewidths=0, zorder=3)
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    finish(ax, xlabel="GCTA $\\hat{h}^2$ (individual-level benchmark)",
           ylabel="summary-statistic $\\hat{h}^2$", title="Per-locus agreement")
    panel_label(ax, "A")
    ax.legend(handles=[Line2D([], [], color=SERIES[e], marker="o", ls="none", ms=5,
                              label=e) for e in others]
                     + [Line2D([], [], color=REFERENCE, ls=(0, (4, 3)), lw=1.1,
                               label="agreement")],
              loc="upper left", handlelength=1.8)

    # B: Bland-Altman for each, against the same benchmark
    ax = axes[1]
    ax.axhline(0, color=REFERENCE, lw=1.1, ls=(0, (4, 3)), zorder=1)
    spreads = []
    for est in others:
        diff = wide_df[col[est]] - wide_df["gcta_h2"]
        mean_est = (wide_df[col[est]] + wide_df["gcta_h2"]) / 2
        md, sd = float(diff.mean()), float(diff.std(ddof=1))
        spreads.append((md, sd))
        ax.scatter(mean_est, diff, s=7, color=SERIES[est], alpha=0.4, linewidths=0, zorder=3)
        ax.axhline(md, color=SERIES[est], lw=1.3, zorder=4)
        ax.text(0.02, md, f"{est} mean {md:+.3f}", transform=ax.get_yaxis_transform(),
                ha="left", va="bottom", fontsize=8, color=INK,
                bbox=dict(facecolor="white", edgecolor="none", pad=1.2, alpha=0.85))
    widest = max((abs(m) + 3 * s for m, s in spreads), default=0.4)
    ax.set_ylim(-widest, widest)
    finish(ax, xlabel="mean of the pair", ylabel="estimator $-$ GCTA",
           title="Bland-Altman vs. GCTA")
    panel_label(ax, "B")

    # C: who is closest to the truth, as a share of loci
    ax = axes[2]
    wins = win_rate(wide_df, ["ncausal"]).sort_values("ncausal")
    x = np.arange(len(wins))
    bottom = np.zeros(len(wins))
    for est in ests:
        share = wins[f"win_{est}"].to_numpy(dtype=float)
        ax.bar(x, share, width=0.62, bottom=bottom, color=SERIES[est], zorder=3,
               edgecolor="white", linewidth=1.8)
        for xi, v, b in zip(x, share, bottom):
            if v > 0.08:
                ax.text(xi, b + v / 2, f"{v:.0%}", ha="center", va="center",
                        fontsize=7.5, color="white")
        bottom = bottom + share
    ax.set_xticks(x)
    ax.set_xticklabels([_fmt(v) for v in wins["ncausal"]])
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_ylim(0, 1.3)  # headroom so the legend never sits on the bars
    ax.grid(axis="x", visible=False)
    finish(ax, xlabel="causal SNPs", ylabel="share of loci",
           title="Which estimator is closest")
    panel_label(ax, "C")
    _estimator_patch_legend(ax, loc="upper left", estimators=ests)

    fig.suptitle("Head-to-head on identical simulated phenotypes", fontsize=11.5,
                 fontweight="bold", color=INK, y=1.03)
    fig.tight_layout()
    return save(fig, out_dir, "fig5_head_to_head")


# ---------------------------------------------------------------------------
# Figure 6 - an extra design axis (n, effect prior, MAF bin, LD architecture)
# ---------------------------------------------------------------------------
def fig_axis_effect(long_df: pd.DataFrame, axis: str, out_dir: str,
                    wide_df: Optional[pd.DataFrame] = None,
                    name: Optional[str] = None) -> List[str]:
    """Bias, RMSE and the paired win-rate across one extra design axis."""
    apply_style()
    levels = _ordered_levels(long_df, axis)
    causal_levels = sorted(long_df["ncausal"].unique())

    fig, axes = plt.subplots(1, 3, figsize=(12.6, 3.9))
    x = np.arange(len(levels))
    width = 0.8 / max(len(ESTIMATOR_ORDER), 1)

    ests = present(long_df)
    width = 0.8 / max(len(ests), 1)
    ax = axes[0]
    ax.axhline(0, color=REFERENCE, lw=1.1, ls=(0, (4, 3)), zorder=1)
    agg = summarize_by(long_df, [axis])
    for k, est in enumerate(ests):
        part = agg[agg["estimator"] == est].set_index(axis).reindex(levels).reset_index()
        off = (k - (len(ests) - 1) / 2) * width
        ax.bar(x + off, part["bias"], width=width * 0.9,
               color=SERIES[est], zorder=3, edgecolor="white", linewidth=1.8)
        ax.errorbar(x + off, part["bias"], yerr=1.96 * part["bias_mcse"],
                    fmt="none", ecolor=INK, elinewidth=1.0, capsize=2.6, zorder=4)
    ax.set_xticks(x)
    ax.set_xticklabels([_fmt(v) for v in levels])
    ax.grid(axis="x", visible=False)
    finish(ax, xlabel=AXIS_LABELS.get(axis, axis), ylabel="bias  ($\\hat{h}^2 - h^2$)",
           title=f"Bias by {AXIS_LABELS.get(axis, axis)}")
    panel_label(ax, "A")
    _headroom(ax)
    _estimator_patch_legend(ax, estimators=ests)

    ax = axes[1]
    for k, est in enumerate(ests):
        part = agg[agg["estimator"] == est].set_index(axis).reindex(levels).reset_index()
        off = (k - (len(ests) - 1) / 2) * width
        ax.bar(x + off, part["rmse"], width=width * 0.9,
               color=SERIES[est], zorder=3, edgecolor="white", linewidth=1.8)
    ax.set_xticks(x)
    ax.set_xticklabels([_fmt(v) for v in levels])
    ax.grid(axis="x", visible=False)
    finish(ax, xlabel=AXIS_LABELS.get(axis, axis), ylabel="RMSE",
           title=f"RMSE by {AXIS_LABELS.get(axis, axis)}")
    panel_label(ax, "B")
    _headroom(ax)
    _estimator_patch_legend(ax, estimators=ests)

    # C: the interaction that matters - does the axis change *who* wins?
    ax = axes[2]
    if wide_df is not None and axis in wide_df.columns:
        wins = win_rate(wide_df, ["ncausal", axis])
        _best_estimator_panel(ax, fig, wins, causal_levels, levels, "ncausal", axis,
                              ests, "causal SNPs", AXIS_LABELS.get(axis, axis))
    else:
        by_causal = summarize_by(long_df, [axis, "ncausal"])
        ax.axhline(0, color=REFERENCE, lw=1.1, ls=(0, (4, 3)), zorder=1)
        for est in ests:
            for j, level in enumerate(levels):
                part = (by_causal[(by_causal["estimator"] == est) & (by_causal[axis] == level)]
                        .sort_values("ncausal"))
                ax.plot(part["ncausal"], part["bias"], color=SERIES[est], lw=1.8,
                        alpha=0.4 + 0.6 * j / max(len(levels) - 1, 1), zorder=3)
        ax.set_xscale("log")
        ax.set_xticks(causal_levels)
        ax.get_xaxis().set_major_formatter(plt.ScalarFormatter())
        finish(ax, xlabel="causal SNPs", ylabel="bias  ($\\hat{h}^2 - h^2$)",
               title="Bias by polygenicity (faint = low level)")
        _estimator_legend(ax, loc="best", estimators=ests)
    panel_label(ax, "C")

    fig.suptitle(f"Effect of {AXIS_LABELS.get(axis, axis)} on estimator accuracy",
                 fontsize=11.5, fontweight="bold", color=INK, y=1.03)
    fig.tight_layout()
    return save(fig, out_dir, name or f"fig6_axis_{axis}")


# ---------------------------------------------------------------------------
# Figure 7 - compute profile
# ---------------------------------------------------------------------------
def fig_compute_profile(wide_df: pd.DataFrame, out_dir: str) -> List[str]:
    """Per-cell cost against region size, and where the time actually goes."""
    apply_style()
    stages = [("t_susie_s", "SuSiE fine-mapping", SERIES["SuSiE"]),
              ("t_gcta_s", "GCTA REML", SERIES["GCTA"]),
              ("t_ld_s", "LD build + Cholesky", SERIES["third"]),
              ("t_sim_s", "simulation", "#8a8880"),
              ("t_ols_s", "marginal sumstats", "#c9c7c0")]
    stages = [s for s in stages if s[0] in wide_df.columns]
    if not stages:
        return []

    # region-level timings are amortised so every bar is "seconds per cell"
    timings = per_cell_timings(wide_df)
    frame = pd.concat([wide_df[["region_id", "p_snps"]], timings], axis=1)
    per_region = (frame.groupby("region_id", observed=True)
                  .agg(p_snps=("p_snps", "first"),
                       **{k: (k, "mean") for k, _, _ in stages})
                  .reset_index())
    per_region["total"] = per_region[[k for k, _, _ in stages]].sum(axis=1)

    fig, axes = plt.subplots(1, 2, figsize=(9.8, 4.0))

    ax = axes[0]
    ax.scatter(per_region["p_snps"], per_region["total"], s=22, color=SERIES["GCTA"],
               edgecolor="white", linewidth=0.6, alpha=0.85, zorder=3)
    ok = (per_region["p_snps"] > 0) & (per_region["total"] > 0)
    if ok.sum() > 3:
        slope, intercept = np.polyfit(np.log(per_region.loc[ok, "p_snps"]),
                                      np.log(per_region.loc[ok, "total"]), 1)
        xs = np.linspace(per_region["p_snps"].min(), per_region["p_snps"].max(), 60)
        ax.plot(xs, np.exp(intercept) * xs**slope, color=INK, lw=1.4, zorder=4)
        ax.text(0.04, 0.93, f"time $\\propto p^{{{slope:.2f}}}$", transform=ax.transAxes,
                fontsize=9, color=INK)
    ax.set_xscale("log")
    ax.set_yscale("log")
    finish(ax, xlabel="SNPs in the cis-region ($p$)", ylabel="seconds per grid cell",
           title="Cost is driven by region size")
    panel_label(ax, "A")

    ax = axes[1]
    totals = [timings[k].mean() for k, _, _ in stages]
    order = np.argsort(totals)[::-1]
    y = np.arange(len(stages))
    for rank, idx in enumerate(order):
        key, label, color = stages[idx]
        ax.barh(rank, totals[idx], height=0.62, color=color, zorder=3,
                edgecolor="white", linewidth=2.0)
        ax.text(totals[idx], rank, f"  {totals[idx]:.2f}s", va="center", fontsize=8.5,
                color=INK_SECONDARY)
    ax.set_yticks(y)
    ax.set_yticklabels([stages[i][1] for i in order])
    ax.invert_yaxis()
    ax.grid(axis="y", visible=False)
    ax.set_xlim(0, max(totals) * 1.25)
    finish(ax, xlabel="mean seconds per grid cell", title="Where the time goes")
    panel_label(ax, "B")

    fig.suptitle("Compute profile of the simulation engine", fontsize=11.5,
                 fontweight="bold", color=INK, y=1.02)
    fig.tight_layout()
    return save(fig, out_dir, "fig7_compute_profile")


# ---------------------------------------------------------------------------
# Cross-campaign comparison
# ---------------------------------------------------------------------------
def fig_condition_grid(conditions, out_dir: str, max_cols: int = 5) -> List[str]:
    """Bias vs. polygenicity, **one panel per level** of every axis swept.

    Panels are never averaged over a campaign's internal axis: a sample-size
    sweep gets one panel per n, an effect-prior sweep one panel per prior, each
    titled with its direction relative to the baseline, so "sample size" can
    never be read without knowing whether it went up or down.
    """
    apply_style()
    n = len(conditions)
    fig, axes, nrows, ncols = _facet_grid(n, width=2.55, height=2.55, max_cols=max_cols)

    ests = present(conditions[0].long) if conditions else list(ESTIMATOR_ORDER)
    for ax, cond in zip(axes, conditions):
        stats = summarize_by(cond.long, ["ncausal"])
        ax.axhline(0, color=REFERENCE, lw=1.1, ls=(0, (4, 3)), zorder=1)
        for est in ests:
            part = stats[stats["estimator"] == est].sort_values("ncausal")
            ax.errorbar(part["ncausal"], part["bias"], yerr=1.96 * part["bias_mcse"],
                        color=SERIES[est], marker="o", ms=4.5, lw=1.9, elinewidth=1.0,
                        capsize=2.2, markeredgecolor="white", markeredgewidth=0.6, zorder=3)
        ax.set_xscale("log")
        ax.set_xticks(sorted(cond.long["ncausal"].unique()))
        ax.get_xaxis().set_major_formatter(plt.ScalarFormatter())
        ax.set_title(cond.label, color=INK, pad=13, fontsize=9.5)
        ax.text(0.5, 1.015, cond.direction, transform=ax.transAxes, ha="center",
                va="bottom", fontsize=8, color=INK_SECONDARY)
        if cond.is_baseline:  # the reference every other panel is read against
            for spine in ("left", "bottom"):
                ax.spines[spine].set_color(INK)
                ax.spines[spine].set_linewidth(1.6)
        finish(ax)

    for ax in axes[n:]:
        ax.set_visible(False)
    for i, ax in enumerate(axes[:n]):
        if i % ncols == 0:
            ax.set_ylabel("bias  ($\\hat{h}^2 - h^2$)")
        # a column can end before the last row; sharex hides those tick labels
        if i + ncols >= n:
            ax.set_xlabel("causal SNPs")
            ax.tick_params(labelbottom=True)
    _estimator_legend(axes[0], loc="lower left", estimators=ests)

    fig.suptitle("Bias by polygenicity — every swept level shown separately",
                 fontsize=12, fontweight="bold", color=INK, y=1.005)
    fig.text(0.5, -0.015,
             "Bold axes mark the baseline panel. Bars are 95% Monte-Carlo intervals.",
             ha="center", fontsize=8.5, color=INK_SECONDARY)
    fig.tight_layout()
    return save(fig, out_dir, "fig8_condition_grid")


def fig_condition_effects(conditions, baseline, out_dir: str) -> List[str]:
    """What each perturbation *does*, stated as a change from the baseline.

    Panel A is the shift in bias relative to the baseline condition; panel B is
    the absolute paired win rate, which is the decision-relevant summary.
    """
    apply_style()
    others = [c for c in conditions if c is not baseline]
    labels = [f"{c.label}  ·  {c.direction}" for c in others]
    y = np.arange(len(others))
    ests = present(baseline.long)
    height = 0.8 / max(len(ests), 1)

    base_stats = summarize_by(baseline.long, []).set_index("estimator")

    fig, axes = plt.subplots(1, 2, figsize=(12.2, 0.46 * len(others) + 2.4),
                             gridspec_kw={"width_ratios": [1.35, 1]})
    letter_y = 1.0 + 0.5 / (0.46 * len(others) + 2.4)  # keep letters beside the titles

    ax = axes[0]
    ax.axvline(0, color=REFERENCE, lw=1.2, ls=(0, (4, 3)), zorder=2)
    for k, est in enumerate(ests):
        deltas, errs = [], []
        for cond in others:
            s = summarize_by(cond.long, []).set_index("estimator")
            deltas.append(float(s.loc[est, "bias"] - base_stats.loc[est, "bias"]))
            errs.append(1.96 * float(np.hypot(s.loc[est, "bias_mcse"],
                                              base_stats.loc[est, "bias_mcse"])))
        offset = (k - (len(ests) - 1) / 2) * height
        ax.barh(y + offset, deltas, height=height * 0.9, color=SERIES[est], zorder=3,
                edgecolor="white", linewidth=1.6)
        ax.errorbar(deltas, y + offset, xerr=errs, fmt="none", ecolor=INK,
                    elinewidth=1.0, capsize=2.4, zorder=4)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.grid(axis="y", visible=False)
    finish(ax, xlabel="change in bias vs. baseline", title="How each change moves the bias")
    ax.text(-0.42, letter_y, "A", transform=ax.transAxes, fontsize=11,
            fontweight="bold", color=INK, ha="left", va="bottom")
    ax.text(0.99, 0.005, "left = more downward bias", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=8, color=INK_SECONDARY)
    _estimator_patch_legend(ax, loc="lower left", estimators=ests)

    # B: a dot plot, not bars - the whole story lives in a narrow window, and
    # bars anchored at zero would compress it into invisibility
    ax = axes[1]
    base_wins = win_rate(baseline.wide, [])
    rates = {est: [] for est in ests}
    errs = {est: [] for est in ests}
    for cond in others:
        w = win_rate(cond.wide, [])
        for est in ests:
            rates[est].append(float(w[f"win_{est}"].iloc[0]))
            errs[est].append(1.96 * float(w[f"win_{est}_mcse"].iloc[0]))

    all_vals = [v for est in ests for v in rates[est]]
    all_err = [v for est in ests for v in errs[est]]
    lo = max(0.0, min(all_vals) - max(all_err) - 0.03)
    hi = min(1.0, max(all_vals) + max(all_err) + 0.05)
    even = 1.0 / len(ests)
    ax.axvline(even, color=REFERENCE, lw=1.2, ls=(0, (4, 3)), zorder=2)
    offsets = np.linspace(-0.22, 0.22, len(ests)) if len(ests) > 1 else [0.0]
    for est, off in zip(ests, offsets):
        ax.errorbar(rates[est], y + off, xerr=errs[est], fmt="o", ms=5.5,
                    color=SERIES[est], ecolor=INK, elinewidth=0.9, capsize=2.4,
                    markeredgecolor="white", markeredgewidth=0.8, zorder=4)
        ax.axvline(float(base_wins[f"win_{est}"].iloc[0]), color=SERIES[est],
                   lw=1.0, alpha=0.55, zorder=1)
    ax.set_yticks(y)
    ax.set_yticklabels([])
    ax.set_xlim(lo, hi)
    ax.set_ylim(len(others) - 0.4, -1.15)
    ax.grid(axis="y", visible=False)
    ax.text(even, -1.05, "even split", ha="center", va="top", fontsize=8,
            color=INK_SECONDARY)
    ax.text(lo + 0.005, -0.5, "faint lines = baseline", ha="left", va="top",
            fontsize=8, color=INK_SECONDARY)
    finish(ax, xlabel="share of loci where this estimator is closest",
           title="Which estimator to prefer")
    ax.text(-0.05, letter_y, "B", transform=ax.transAxes, fontsize=11,
            fontweight="bold", color=INK, ha="left", va="bottom")

    fig.suptitle("Effect of each experimental change, relative to the baseline",
                 fontsize=12, fontweight="bold", color=INK, y=1.0)
    fig.tight_layout()
    return save(fig, out_dir, "fig9_condition_effects")


# ---------------------------------------------------------------------------
def _fmt(v) -> str:
    if isinstance(v, (int, np.integer)):
        return str(int(v))
    if isinstance(v, (float, np.floating)):
        return f"{v:g}"
    return str(v)


def _ordered_levels(df: pd.DataFrame, axis: str) -> List:
    vals = df[axis].dropna().unique().tolist()
    preferred = {
        "effect_dist": ["normal", "fixed", "pareto"],
        "maf_bin": ["all", "common", "rare"],
        "ld_arch": ["low", "regular", "high"],
    }
    if axis in preferred:
        return [v for v in preferred[axis] if v in vals] + [v for v in vals if v not in preferred[axis]]
    return sorted(vals)
