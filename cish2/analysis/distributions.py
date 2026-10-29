"""Distribution figures: the full spread of every estimate, not just its mean.

The metric figures elsewhere in this package report means with Monte-Carlo
intervals, which is the right summary for *bias* but hides everything about
*spread* -- and spread is most of what separates these estimators at small
sample sizes.  These are box-and-whisker views in the style the project used
before: one group per true heritability, one box per estimator, the median as a
solid line and the **mean as a dashed line** so the two can be compared at a
glance (they separate exactly when the distribution is skewed, which is the
interesting case).

Three views of the same grid, in increasing order of diagnostic usefulness:

``estimates``       what each estimator returns, against a marker for the truth
``error``           estimate - truth, so every panel shares a zero reference
``relative error``  (estimate - truth) / truth, which is the scale that matters
                    when the true h2 is 0.05 and an absolute error of 0.03 is
                    a 60% miss rather than a rounding difference
"""

from __future__ import annotations

from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from .figures import present
from .style import (
    INK,
    INK_MUTED,
    INK_SECONDARY,
    REFERENCE,
    SERIES,
    apply_style,
    finish,
    save,
)

# the mean is an annotation on the box, not a data series, so it wears ink
MEAN_STYLE = dict(linestyle=(0, (3, 2)), linewidth=1.6, color=INK)
# white reads on every fill colour, but needs an outline to survive both the
# lightest box and a white legend background
MEDIAN_STYLE = dict(linewidth=1.5, color="white",
                    path_effects=[pe.withStroke(linewidth=3.0, foreground=INK_SECONDARY)])
FLIER_STYLE = dict(marker="o", markersize=2.2, markerfacecolor=INK_MUTED,
                   markeredgecolor="none", alpha=0.45)


def _grouped_boxes(ax, groups: List[List[np.ndarray]], estimators: Sequence[str],
                   width_frac: float = 0.78, showfliers: bool = True) -> None:
    """Draw one box per estimator inside each group position.

    ``groups[i][k]`` is the sample for group ``i``, estimator ``k``.
    """
    n_est = len(estimators)
    box_w = width_frac / n_est
    for k, est in enumerate(estimators):
        offsets = (k - (n_est - 1) / 2) * box_w
        data, pos = [], []
        for i, group in enumerate(groups):
            values = np.asarray(group[k], dtype=float)
            values = values[np.isfinite(values)]
            if values.size == 0:
                continue
            data.append(values)
            pos.append(i + offsets)
        if not data:
            continue
        bp = ax.boxplot(
            data, positions=pos, widths=box_w * 0.86, patch_artist=True,
            showmeans=True, meanline=True, showfliers=showfliers,
            meanprops=MEAN_STYLE, medianprops=MEDIAN_STYLE,
            flierprops=FLIER_STYLE,
            boxprops=dict(linewidth=0.7, edgecolor=INK_SECONDARY),
            whiskerprops=dict(linewidth=0.8, color=INK_SECONDARY),
            capprops=dict(linewidth=0.8, color=INK_SECONDARY),
            manage_ticks=False, zorder=3,
        )
        for box in bp["boxes"]:
            box.set_facecolor(SERIES[est])


def _box_legend(ax, estimators: Sequence[str], loc: str = "upper left",
                extra: Optional[List] = None) -> None:
    handles = [Patch(facecolor=SERIES[e], edgecolor=INK_SECONDARY, label=e)
               for e in estimators]
    handles += [
        Line2D([], [], color="white", lw=1.6, label="median",
               path_effects=[pe.withStroke(linewidth=3.0, foreground=INK_SECONDARY)]),
        Line2D([], [], color=INK, lw=1.6, ls=(0, (3, 2)), label="mean"),
    ]
    leg = ax.legend(handles=(extra or []) + handles, loc=loc, fontsize=7.5,
                    handlelength=1.6, framealpha=0.9, frameon=True)
    leg.get_frame().set_edgecolor(REFERENCE)
    leg.get_frame().set_linewidth(0.6)


def _panels(n: int, width: float = 3.4, height: float = 3.2, max_cols: int = 5,
            sharey: bool = True):
    ncols = min(max_cols, n)
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(width * ncols + 0.8, height * nrows + 0.8),
                             squeeze=False, sharey=sharey)
    return fig, axes.ravel(), nrows, ncols


def _collect(long_df: pd.DataFrame, facet_value, facet: str,
             groups: Sequence, group_col: str, estimators: Sequence[str],
             value_col: str) -> List[List[np.ndarray]]:
    part = long_df[long_df[facet] == facet_value]
    out = []
    for g in groups:
        sub = part[part[group_col] == g]
        out.append([sub.loc[sub["estimator"] == e, value_col].to_numpy(dtype=float)
                    for e in estimators])
    return out


# ---------------------------------------------------------------------------
def fig_estimate_boxplots(long_df: pd.DataFrame, out_dir: str,
                          name: str = "figD1_estimate_distributions") -> List[str]:
    """What each estimator actually returns, across the whole grid."""
    apply_style()
    ests = present(long_df)
    causal = sorted(long_df["ncausal"].dropna().unique())
    h2 = sorted(long_df["h2_true"].dropna().unique())

    fig, axes, nrows, ncols = _panels(len(causal))
    for idx, (ax, nc) in enumerate(zip(axes, causal)):
        groups = _collect(long_df, nc, "ncausal", h2, "h2_true", ests, "h2_hat")
        # the truth traced across the panel: a rule inside each group would be
        # hidden by the boxes it is meant to be compared against
        ax.plot(range(len(h2)), h2, color=REFERENCE, ls=(0, (4, 3)), lw=1.3,
                marker="_", ms=26, mew=1.6, zorder=1)
        _grouped_boxes(ax, groups, ests)
        ax.set_xticks(range(len(h2)))
        ax.set_xticklabels([f"{v:g}" for v in h2])
        ax.set_xlim(-0.6, len(h2) - 0.4)
        ax.grid(axis="x", visible=False)
        finish(ax, title=f"{nc} causal SNP{'s' if nc != 1 else ''}")
        if idx % ncols == 0:
            ax.set_ylabel("estimated $h^2$")
        if idx >= len(causal) - ncols:
            ax.set_xlabel("true $h^2$")

    for ax in axes[len(causal):]:
        ax.set_visible(False)
    _box_legend(axes[0], ests, loc="upper left",
                extra=[Line2D([], [], color=REFERENCE, ls=(0, (4, 3)), lw=1.2,
                              label="true $h^2$")])
    fig.suptitle("Distribution of estimates across replicate loci",
                 fontsize=12, fontweight="bold", color=INK, y=1.0)
    fig.text(0.5, -0.015,
             "Boxes span the interquartile range; whiskers 1.5x IQR; dots are "
             "individual loci beyond that.",
             ha="center", fontsize=8.5, color=INK_SECONDARY)
    fig.tight_layout()
    return save(fig, out_dir, name)


def fig_error_boxplots(long_df: pd.DataFrame, out_dir: str,
                       name: str = "figD2_error_distributions") -> List[str]:
    """Estimate minus truth, so every panel shares one zero reference."""
    apply_style()
    ests = present(long_df)
    causal = sorted(long_df["ncausal"].dropna().unique())
    h2 = sorted(long_df["h2_true"].dropna().unique())
    frame = long_df.assign(err=long_df["h2_hat"] - long_df["h2_true"])

    fig, axes, nrows, ncols = _panels(len(causal))
    for idx, (ax, nc) in enumerate(zip(axes, causal)):
        ax.axhline(0, color=REFERENCE, lw=1.2, ls=(0, (4, 3)), zorder=1)
        groups = _collect(frame, nc, "ncausal", h2, "h2_true", ests, "err")
        _grouped_boxes(ax, groups, ests)
        ax.set_xticks(range(len(h2)))
        ax.set_xticklabels([f"{v:g}" for v in h2])
        ax.set_xlim(-0.6, len(h2) - 0.4)
        ax.grid(axis="x", visible=False)
        finish(ax, title=f"{nc} causal SNP{'s' if nc != 1 else ''}")
        if idx % ncols == 0:
            ax.set_ylabel("error  ($\\hat{h}^2 - h^2$)")
        if idx >= len(causal) - ncols:
            ax.set_xlabel("true $h^2$")

    for ax in axes[len(causal):]:
        ax.set_visible(False)
    # a few runaway fits would otherwise flatten every box into a line
    lo, hi = np.nanpercentile(frame["err"].to_numpy(dtype=float), [0.5, 99.5])
    pad = 0.12 * (hi - lo)
    axes[0].set_ylim(lo - pad, hi + pad)
    _box_legend(axes[0], ests, loc="lower left")
    fig.suptitle("Distribution of errors — where each estimator sits relative to the truth",
                 fontsize=12, fontweight="bold", color=INK, y=1.0)
    fig.text(0.5, -0.015,
             "Axis trimmed to the 0.5-99.5 percentile of the error; a handful of "
             "runaway fits sit outside it.",
             ha="center", fontsize=8.5, color=INK_SECONDARY)
    fig.tight_layout()
    return save(fig, out_dir, name)


def fig_relative_error_boxplots(long_df: pd.DataFrame, out_dir: str,
                                name: str = "figD3_relative_error") -> List[str]:
    """Percent error, which is the scale that matters at low true h2."""
    apply_style()
    ests = present(long_df)
    causal = sorted(long_df["ncausal"].dropna().unique())
    h2 = sorted(long_df["h2_true"].dropna().unique())
    denom = long_df["h2_true"].abs().clip(lower=1e-12)
    frame = long_df.assign(pe=100.0 * (long_df["h2_hat"] - long_df["h2_true"]) / denom)

    fig, axes, nrows, ncols = _panels(len(causal))
    for idx, (ax, nc) in enumerate(zip(axes, causal)):
        ax.axhline(0, color=REFERENCE, lw=1.2, ls=(0, (4, 3)), zorder=1)
        ax.axhline(-100, color=INK_MUTED, lw=0.9, ls=(0, (2, 3)), zorder=1)
        groups = _collect(frame, nc, "ncausal", h2, "h2_true", ests, "pe")
        _grouped_boxes(ax, groups, ests)
        ax.set_xticks(range(len(h2)))
        ax.set_xticklabels([f"{v:g}" for v in h2])
        ax.set_xlim(-0.6, len(h2) - 0.4)
        ax.grid(axis="x", visible=False)
        finish(ax, title=f"{nc} causal SNP{'s' if nc != 1 else ''}")
        if idx % ncols == 0:
            ax.set_ylabel("percent error (%)")
        if idx >= len(causal) - ncols:
            ax.set_xlabel("true $h^2$")

    for ax in axes[len(causal):]:
        ax.set_visible(False)
    axes[0].set_ylim(-160, 200)
    axes[0].text(len(h2) - 0.5, -100, " estimate = 0", ha="right", va="bottom",
                 fontsize=7.5, color=INK_MUTED)
    _box_legend(axes[0], ests, loc="upper left")
    fig.suptitle("Relative error — the same errors as a fraction of the truth",
                 fontsize=12, fontweight="bold", color=INK, y=1.0)
    fig.text(0.5, -0.015,
             "The lower dotted rule at -100% is what an estimate of exactly zero "
             "looks like. Axis clipped to [-160, 200]%.",
             ha="center", fontsize=8.5, color=INK_SECONDARY)
    fig.tight_layout()
    return save(fig, out_dir, name)


def fig_axis_boxplots(long_df: pd.DataFrame, axis: str, out_dir: str,
                      name: Optional[str] = None) -> List[str]:
    """Error distributions across the levels of one extra design axis."""
    from .figures import AXIS_LABELS, _ordered_levels

    apply_style()
    ests = present(long_df)
    levels = _ordered_levels(long_df, axis)
    h2 = sorted(long_df["h2_true"].dropna().unique())
    frame = long_df.assign(err=long_df["h2_hat"] - long_df["h2_true"])

    fig, axes, nrows, ncols = _panels(len(levels), max_cols=5)
    for idx, (ax, level) in enumerate(zip(axes, levels)):
        ax.axhline(0, color=REFERENCE, lw=1.2, ls=(0, (4, 3)), zorder=1)
        groups = _collect(frame, level, axis, h2, "h2_true", ests, "err")
        _grouped_boxes(ax, groups, ests)
        ax.set_xticks(range(len(h2)))
        ax.set_xticklabels([f"{v:g}" for v in h2])
        ax.set_xlim(-0.6, len(h2) - 0.4)
        ax.grid(axis="x", visible=False)
        label = AXIS_LABELS.get(axis, axis)
        finish(ax, title=f"{label} = {level}")
        if idx % ncols == 0:
            ax.set_ylabel("error  ($\\hat{h}^2 - h^2$)")
        if idx >= len(levels) - ncols:
            ax.set_xlabel("true $h^2$")

    for ax in axes[len(levels):]:
        ax.set_visible(False)
    lo, hi = np.nanpercentile(frame["err"].to_numpy(dtype=float), [0.5, 99.5])
    pad = 0.12 * (hi - lo)
    axes[0].set_ylim(lo - pad, hi + pad)
    _box_legend(axes[0], ests, loc="lower left")
    fig.suptitle(f"Error distributions across {AXIS_LABELS.get(axis, axis)}",
                 fontsize=12, fontweight="bold", color=INK, y=1.0)
    fig.tight_layout()
    return save(fig, out_dir, name or f"figD4_boxplots_{axis}")


# ---------------------------------------------------------------------------
def fig_gtex_boxplots(df: pd.DataFrame, out_dir: str,
                      name: str = "gtex_figD1_distributions") -> List[str]:
    """Real-data distributions: estimates per sample size, and the paired gaps.

    Truth is unavailable, so the second panel shows each summary-statistic
    estimator minus REML on the *same gene* -- a paired difference, which is the
    real-data analogue of an error distribution.
    """
    apply_style()
    order = ["GCTA", "SuSiE", "HEELS"]
    col = {"GCTA": "gcta_h2", "SuSiE": "susie_h2", "HEELS": "heels_h2"}
    ests = [e for e in order if col[e] in df.columns]
    levels = sorted(df["n"].dropna().unique())

    fig, axes = plt.subplots(1, 2, figsize=(12.6, 4.6))

    ax = axes[0]
    groups = [[df.loc[df["n"] == n, col[e]].to_numpy(dtype=float) for e in ests]
              for n in levels]
    _grouped_boxes(ax, groups, ests, showfliers=False)
    ax.set_xticks(range(len(levels)))
    ax.set_xticklabels([f"{int(v)}" for v in levels])
    ax.set_xlim(-0.6, len(levels) - 0.4)
    ax.set_ylim(-0.02, 0.30)
    ax.grid(axis="x", visible=False)
    finish(ax, xlabel="individuals used", ylabel="estimated $h^2$",
           title="What each estimator returns")
    _box_legend(ax, ests, loc="upper right")

    ax = axes[1]
    others = [e for e in ests if e != "GCTA"]
    ax.axhline(0, color=REFERENCE, lw=1.2, ls=(0, (4, 3)), zorder=1)
    groups = [[(df.loc[df["n"] == n, col[e]] - df.loc[df["n"] == n, "gcta_h2"]).to_numpy(dtype=float)
               for e in others] for n in levels]
    _grouped_boxes(ax, groups, others, showfliers=False)
    ax.set_xticks(range(len(levels)))
    ax.set_xticklabels([f"{int(v)}" for v in levels])
    ax.set_xlim(-0.6, len(levels) - 0.4)
    ax.set_ylim(-0.10, 0.06)
    ax.grid(axis="x", visible=False)
    finish(ax, xlabel="individuals used", ylabel="estimator $-$ GCTA (same gene)",
           title="Paired disagreement with the individual-level benchmark")
    _box_legend(ax, others, loc="lower right")

    fig.suptitle("GTEx Whole Blood — distributions over ~19,000 genes",
                 fontsize=12, fontweight="bold", color=INK, y=1.0)
    fig.text(0.5, -0.02,
             "Whiskers span 1.5x IQR; individual outlier genes are omitted "
             "(19,000 per box would obscure everything). Dashed line inside each "
             "box is the mean, solid line the median. Axes trimmed.",
             ha="center", fontsize=8.5, color=INK_SECONDARY)
    fig.tight_layout()
    return save(fig, out_dir, name)
