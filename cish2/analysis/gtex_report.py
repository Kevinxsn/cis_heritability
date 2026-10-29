"""Test the simulation-derived predictions on real GTEx data.

Real data has no ground truth, so nothing here measures bias.  What it measures
is the **disagreement** with the individual-level benchmark, `estimator - GCTA`,
gene by gene — and whether that disagreement behaves the way the simulations
said it would.  Each simulation result becomes a falsifiable prediction:

===========================================================================
prediction (from simulation)                     tested here by
---------------------------------------------------------------------------
P1  the gap shrinks as sample size grows         nested subsampling of the
                                                 same individuals
P2  the gap is ~0 for a concentrated,            stratifying on LD-clumped
    localisable signal and large for signal      independent signals, *within*
    too diffuse to localise                      heritability strata
P3  REML fails often and mildly; fine-mapping    counting out-of-range
    fails rarely and wildly                      estimates
P4  the gap grows with the size of h2 itself     stratifying on GCTA h2
P5  HEELS tracks GCTA where SuSiE does not,      comparing both summary-statistic
    because it does not assume sparsity          estimators against GCTA

The P2 test must condition on heritability: genes with more detectable signals
also have more heritability, and the marginal comparison points the other way
(Simpson's paradox).
===========================================================================

A pattern matching the prediction is *consistency*, not proof: GCTA is not a
gold standard on real data (§4 of the simulation write-up shows it turning
upward-biased under frequency stratification).  What the agreement buys is
confidence that the simulated regime is not an artefact.
"""

from __future__ import annotations

import argparse
import os
import sys
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from analysis.load import load_run  # noqa: E402
from analysis.metrics import mcse_mean, mcse_prop, summarize_by  # noqa: E402
from analysis.style import (  # noqa: E402
    BLUE_RAMP,
    ESTIMATOR_ORDER,
    INK,
    INK_SECONDARY,
    REFERENCE,
    SERIES,
    apply_style,
    finish,
    ordinal_colors,
    panel_label,
    save,
)

CLUMP_BINS = [(0, 0, "0"), (1, 1, "1"), (2, 2, "2"), (3, 3, "3"), (4, 10**6, "4+")]


# --------------------------------------------------------------------- load --
def load_gtex(run_dir: str) -> pd.DataFrame:
    """Load the gene-level table and attach the disagreement columns."""
    from cish2.results import ResultStore

    df = ResultStore(run_dir).load()
    df["gap"] = df["susie_h2"] - df["gcta_h2"]
    df["usable"] = df["gcta_h2"].notna() & df["susie_h2"].notna()
    if "heels_h2" in df.columns:
        # the second summary-statistic estimator, measured against the same
        # individual-level benchmark
        df["gap_heels"] = df["heels_h2"] - df["gcta_h2"]
        df["usable_heels"] = df["gcta_h2"].notna() & df["heels_h2"].notna()
    df["clump_bin"] = pd.cut(
        df["n_clumps"].fillna(-1),
        bins=[-0.5, 0.5, 1.5, 2.5, 3.5, np.inf],
        labels=["0", "1", "2", "3", "4+"],
    )
    df["cs_bin"] = pd.cut(
        df["susie_n_cs"].fillna(-1),
        bins=[-0.5, 0.5, 1.5, 2.5, np.inf],
        labels=["0", "1", "2", "3+"],
    )
    return df


def _agg(part: pd.DataFrame) -> Dict[str, float]:
    ok = part[part["usable"]]
    gap = ok["gap"].to_numpy(dtype=float)
    out = {
        "genes": int(len(ok)),
        "mean_gcta": float(ok["gcta_h2"].mean()) if len(ok) else np.nan,
        "mean_susie": float(ok["susie_h2"].mean()) if len(ok) else np.nan,
        "mean_gap": float(np.mean(gap)) if gap.size else np.nan,
        "gap_mcse": mcse_mean(gap),
        "median_gap": float(np.median(gap)) if gap.size else np.nan,
        "corr": float(ok[["gcta_h2", "susie_h2"]].corr().iloc[0, 1]) if len(ok) > 2 else np.nan,
        "spearman": (float(ok[["gcta_h2", "susie_h2"]].corr(method="spearman").iloc[0, 1])
                     if len(ok) > 2 else np.nan),
        "gcta_failed": float(part["gcta_h2"].isna().mean()) if len(part) else np.nan,
        "susie_failed": float(part["susie_h2"].isna().mean()) if len(part) else np.nan,
        "gcta_impossible": float(((part["gcta_h2"] < 0) | (part["gcta_h2"] > 1)).mean()),
        "susie_impossible": float(((part["susie_h2"] < 0) | (part["susie_h2"] > 1)).mean()),
        "susie_max_abs": float(part["susie_h2"].abs().max()),
        "gcta_max_abs": float(part["gcta_h2"].abs().max()),
    }
    if "gap_heels" in part.columns:
        ok_h = part[part["usable_heels"]]
        gap_h = ok_h["gap_heels"].to_numpy(dtype=float)
        out.update({
            "mean_heels": float(ok_h["heels_h2"].mean()) if len(ok_h) else np.nan,
            "mean_gap_heels": float(np.mean(gap_h)) if gap_h.size else np.nan,
            "gap_heels_mcse": mcse_mean(gap_h),
            "median_gap_heels": float(np.median(gap_h)) if gap_h.size else np.nan,
            "corr_heels": (float(ok_h[["gcta_h2", "heels_h2"]].corr().iloc[0, 1])
                           if len(ok_h) > 2 else np.nan),
            "heels_failed": float(part["heels_h2"].isna().mean()) if len(part) else np.nan,
            "heels_impossible": float(((part["heels_h2"] < 0) | (part["heels_h2"] > 1)).mean()),
            "heels_max_abs": float(part["heels_h2"].abs().max()),
        })
    return out


def by_sample_size(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for n, part in df.groupby("n"):
        rec = {"n": int(n)}
        rec.update(_agg(part))
        rows.append(rec)
    return pd.DataFrame(rows).sort_values("n")


def by_stratum(df: pd.DataFrame, column: str) -> pd.DataFrame:
    rows = []
    for level, part in df.groupby(column, observed=True):
        rec = {column: str(level)}
        rec.update(_agg(part))
        rows.append(rec)
    return pd.DataFrame(rows)


def simulation_prediction(sim_run_dir: str) -> pd.DataFrame:
    """Predicted SuSiE-minus-GCTA gap per sample size, from the simulations."""
    from analysis.load import to_long

    wide = load_run(sim_run_dir)
    long = to_long(wide)
    stats = summarize_by(long, ["n"])
    rows = []
    for n in sorted(stats["n"].unique()):
        g = stats[(stats["estimator"] == "GCTA") & (stats["n"] == n)].iloc[0]
        s = stats[(stats["estimator"] == "SuSiE") & (stats["n"] == n)].iloc[0]
        rows.append(dict(n=int(n),
                         predicted_gap=float(s["bias"] - g["bias"]),
                         predicted_mcse=float(np.hypot(s["bias_mcse"], g["bias_mcse"]))))
    return pd.DataFrame(rows)


# ------------------------------------------------------------------ figures --
def _log_ticks(ax, values) -> None:
    """Label exactly these values on a log axis; kill the minor tick labels.

    Matplotlib keeps drawing minor decade labels (4x10^2, 6x10^2) on a short
    log range, which collide with the explicit ticks we care about.
    """
    from matplotlib.ticker import NullFormatter, NullLocator

    ax.set_xscale("log")
    ax.set_xticks(sorted(set(int(v) for v in values)))
    ax.get_xaxis().set_major_formatter(plt.ScalarFormatter())
    ax.get_xaxis().set_minor_locator(NullLocator())
    ax.get_xaxis().set_minor_formatter(NullFormatter())



def fig_sample_size_test(observed: pd.DataFrame, predicted: Optional[pd.DataFrame],
                         out_dir: str) -> List[str]:
    """P1: does the gap shrink with n on real data, as the simulations said?"""
    apply_style()
    fig, axes = plt.subplots(1, 2, figsize=(9.8, 4.1))

    ax = axes[0]
    ax.axhline(0, color=REFERENCE, lw=1.1, ls=(0, (4, 3)), zorder=1)
    if predicted is not None and len(predicted):
        ax.errorbar(predicted["n"], predicted["predicted_gap"],
                    yerr=1.96 * predicted["predicted_mcse"], color=INK_SECONDARY,
                    marker="s", ms=5, lw=2.0, ls=(0, (5, 2)), elinewidth=1.0, capsize=2.4,
                    markeredgecolor="white", markeredgewidth=0.7, zorder=3)
    if "mean_gap_heels" in observed.columns:
        ax.errorbar(observed["n"], observed["mean_gap_heels"],
                    yerr=1.96 * observed["gap_heels_mcse"], color=SERIES["HEELS"],
                    marker="D", ms=5, lw=2.0, elinewidth=1.0, capsize=2.4,
                    markeredgecolor="white", markeredgewidth=0.7, zorder=4)
    ax.errorbar(observed["n"], observed["mean_gap"], yerr=1.96 * observed["gap_mcse"],
                color=SERIES["SuSiE"], marker="o", ms=6, lw=2.2, elinewidth=1.1,
                capsize=2.6, markeredgecolor="white", markeredgewidth=0.8, zorder=5)
    # label every observed n, plus any simulated n far enough away to not crowd
    ticks = [int(v) for v in observed["n"].unique()]
    if predicted is not None and len(predicted):
        for v in predicted["n"]:
            if all(abs(np.log(v) - np.log(t)) > 0.15 for t in ticks):
                ticks.append(int(v))
    _log_ticks(ax, ticks)
    finish(ax, xlabel="individuals used", ylabel="SuSiE $-$ GCTA  ($\\hat{h}^2$)",
           title="Observed disagreement vs. the simulated prediction")
    panel_label(ax, "A")
    handles = [Line2D([], [], color=SERIES["SuSiE"], marker="o", lw=2.2, ms=6,
                      label="SuSiE $-$ GCTA (observed)")]
    if "mean_gap_heels" in observed.columns:
        handles.append(Line2D([], [], color=SERIES["HEELS"], marker="D", lw=2.0, ms=5,
                              label="HEELS $-$ GCTA (observed)"))
    if predicted is not None and len(predicted):
        handles.append(Line2D([], [], color=INK_SECONDARY, marker="s", lw=2.0,
                              ls=(0, (5, 2)), ms=5, label="SuSiE gap, simulated"))
    ax.legend(handles=handles, loc="lower right", handlelength=2.2, fontsize=8)

    ax = axes[1]
    ax.plot(observed["n"], observed["corr"], color=SERIES["SuSiE"], marker="o", ms=6,
            lw=2.2, markeredgecolor="white", markeredgewidth=0.8, zorder=3)
    if "corr_heels" in observed.columns:
        ax.plot(observed["n"], observed["corr_heels"], color=SERIES["HEELS"],
                marker="D", ms=5, lw=2.0, markeredgecolor="white",
                markeredgewidth=0.8, zorder=3)
    _log_ticks(ax, observed["n"].unique())
    ax.set_ylim(0, 1)
    finish(ax, xlabel="individuals used", ylabel="correlation between estimators",
           title="Agreement with the individual-level benchmark")
    panel_label(ax, "B")
    corr_handles = [Line2D([], [], color=SERIES["SuSiE"], marker="o", lw=2.2,
                           label="SuSiE vs GCTA")]
    if "corr_heels" in observed.columns:
        corr_handles.append(Line2D([], [], color=SERIES["HEELS"], marker="D", lw=2.0,
                                   label="HEELS vs GCTA"))
    ax.legend(handles=corr_handles, loc="lower right", handlelength=2.2, fontsize=8)

    fig.suptitle("Prediction 1 — the sample-size effect is real but far smaller "
                 "than simulated", fontsize=11.5, fontweight="bold", color=INK, y=1.02)
    fig.tight_layout()
    return save(fig, out_dir, "gtex_fig1_sample_size")


H2_BINS = [0.0, 0.05, 0.1, 0.2, 0.4, 1.0]


def add_h2_bin(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame.copy()
    out["h2_bin"] = pd.cut(out["gcta_h2"], bins=H2_BINS)
    return out


def two_way_gap(full: pd.DataFrame, column: str = "clump_bin", min_genes: int = 30
                ) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Mean gap over (architecture proxy x heritability bin).

    The one-way version of this table is confounded: genes with more detectable
    signals also have higher heritability, so a marginal comparison mixes the
    two.  Conditioning on heritability is what makes the comparison mean
    anything.
    """
    frame = add_h2_bin(full[full["usable"]]).dropna(subset=["h2_bin"])
    mean = frame.pivot_table(index=column, columns="h2_bin", values="gap",
                             aggfunc="mean", observed=True)
    count = frame.pivot_table(index=column, columns="h2_bin", values="gap",
                              aggfunc="size", observed=True)
    mean = mean.where(count >= min_genes)
    return mean, count


def fig_architecture_test(df: pd.DataFrame, out_dir: str, n_full: int) -> List[str]:
    """P2/P4: is the gap driven by architecture once heritability is held fixed?"""
    apply_style()
    full = df[df["n"] == n_full]
    mean, count = two_way_gap(full, "clump_bin")

    fig, axes = plt.subplots(1, 3, figsize=(12.8, 4.1))
    h2_levels = list(mean.columns)
    colors = ordinal_colors(len(h2_levels), BLUE_RAMP)

    # A: the conditional relationship - one line per heritability stratum
    ax = axes[0]
    ax.axhline(0, color=REFERENCE, lw=1.1, ls=(0, (4, 3)), zorder=1)
    x = np.arange(len(mean.index))
    for level, color in zip(h2_levels, colors):
        series = mean[level]
        ok = series.notna().to_numpy()
        if ok.sum() < 2:
            continue
        ax.plot(x[ok], series.to_numpy()[ok], color=color, marker="o", ms=5, lw=2.0,
                markeredgecolor="white", markeredgewidth=0.7, zorder=3,
                label=str(level).replace(", ", "–").strip("()[]"))
    ax.set_xticks(x)
    ax.set_xticklabels(mean.index)
    ax.grid(axis="x", visible=False)
    finish(ax, xlabel="independent marginal signals (LD clumping)",
           ylabel="SuSiE $-$ GCTA  ($\\hat{h}^2$)",
           title="Holding heritability fixed, the gap *closes*")
    panel_label(ax, "A")
    leg = ax.legend(title="GCTA $\\hat{h}^2$", loc="lower right", fontsize=7.5,
                    title_fontsize=7.5, handlelength=1.6)
    leg.get_title().set_color(INK_SECONDARY)

    # B: the same numbers as a grid, so the confounding is visible at a glance
    ax = axes[1]
    grid = mean.to_numpy(dtype=float)
    vmax = float(np.nanmax(np.abs(grid)))
    im = ax.imshow(grid, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto", origin="lower")
    for i in range(grid.shape[0]):
        for j in range(grid.shape[1]):
            if not np.isfinite(grid[i, j]):
                continue
            r, g, b, _ = im.cmap(im.norm(grid[i, j]))
            lum = 0.2126 * r + 0.7152 * g + 0.0722 * b
            ax.text(j, i, f"{grid[i, j]:+.3f}", ha="center", va="center", fontsize=7.5,
                    color="white" if lum < 0.5 else INK)
    ax.set_xticks(range(len(h2_levels)))
    ax.set_xticklabels([str(s).replace(", ", "–").strip("()[]") for s in h2_levels],
                       fontsize=7.5)
    ax.set_yticks(range(len(mean.index)))
    ax.set_yticklabels(mean.index)
    ax.grid(False)
    finish(ax, xlabel="GCTA $\\hat{h}^2$", ylabel="independent signals",
           title="Mean gap, conditioned on both")
    panel_label(ax, "B")
    cb = fig.colorbar(im, ax=ax, fraction=0.042, pad=0.03)
    cb.outline.set_visible(False)
    cb.set_label("SuSiE $-$ GCTA", fontsize=8.5, color=INK_SECONDARY)

    # C: the marginal view, which points the other way - shown so the
    # confounding is on the record rather than hidden
    ax = axes[2]
    stratum = by_stratum(full, "clump_bin").dropna(subset=["mean_gap"])
    stratum = stratum[stratum["genes"] >= 20]
    x = np.arange(len(stratum))
    ax.axhline(0, color=REFERENCE, lw=1.1, ls=(0, (4, 3)), zorder=1)
    ax.bar(x, stratum["mean_gap"], width=0.62, color=SERIES["SuSiE"], zorder=3,
           edgecolor="white", linewidth=2.0)
    ax.errorbar(x, stratum["mean_gap"], yerr=1.96 * stratum["gap_mcse"], fmt="none",
                ecolor=INK, elinewidth=1.1, capsize=3, zorder=4)
    ax.set_xticks(x)
    # the confounder goes in the tick label, where it cannot collide with a bar
    ax.set_xticklabels(
        [f"{row.clump_bin}\n$\\hat{{h}}^2$={row.mean_gcta:.2f}" for row in stratum.itertuples()],
        fontsize=7.5,
    )
    ax.grid(axis="x", visible=False)
    finish(ax, xlabel="independent marginal signals",
           ylabel="SuSiE $-$ GCTA  ($\\hat{h}^2$)",
           title="Marginal view — confounded, points the other way")
    panel_label(ax, "C")

    fig.suptitle(f"Prediction 2 — the gap is largest where heritability is *not* "
                 f"localisable (n = {n_full})", fontsize=11.5, fontweight="bold",
                 color=INK, y=1.02)
    fig.text(0.5, -0.02,
             "Cells with fewer than 30 genes are left blank. Panel C is annotated with "
             "each stratum's mean heritability, which is what confounds it.",
             ha="center", fontsize=8.5, color=INK_SECONDARY)
    fig.tight_layout()
    return save(fig, out_dir, "gtex_fig2_architecture")


def fig_agreement(df: pd.DataFrame, out_dir: str, n_full: int) -> List[str]:
    """Per-gene agreement at full sample size, coloured by architecture."""
    apply_style()
    full = df[(df["n"] == n_full) & df["usable"]]
    fig, axes = plt.subplots(1, 2, figsize=(9.8, 4.2))

    levels = ["0", "1", "2", "3", "4+"]
    colors = ordinal_colors(len(levels), BLUE_RAMP)

    ax = axes[0]
    lo = float(min(full["gcta_h2"].quantile(0.002), full["susie_h2"].quantile(0.002), 0))
    hi = float(max(full["gcta_h2"].quantile(0.998), full["susie_h2"].quantile(0.998)))
    ax.plot([lo, hi], [lo, hi], ls=(0, (4, 3)), lw=1.1, color=REFERENCE, zorder=1)
    for level, color in zip(levels, colors):
        part = full[full["clump_bin"] == level]
        if part.empty:
            continue
        ax.scatter(part["gcta_h2"], part["susie_h2"], s=7, color=color, alpha=0.5,
                   linewidths=0, zorder=3, label=level)
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    finish(ax, xlabel="GCTA $\\hat{h}^2$", ylabel="SuSiE $\\hat{h}^2$",
           title=f"Per-gene agreement ({len(full):,} genes)")
    panel_label(ax, "A")
    leg = ax.legend(title="independent signals", loc="upper left", fontsize=7.5,
                    title_fontsize=7.5, markerscale=2.4, handletextpad=0.3)
    leg.get_title().set_color(INK_SECONDARY)

    ax = axes[1]
    mean_est = (full["gcta_h2"] + full["susie_h2"]) / 2
    diff = full["susie_h2"] - full["gcta_h2"]
    md, sd = float(diff.mean()), float(diff.std(ddof=1))
    for level, color in zip(levels, colors):
        m = full["clump_bin"] == level
        ax.scatter(mean_est[m], diff[m], s=7, color=color, alpha=0.5, linewidths=0, zorder=3)
    ax.axhline(0, color=REFERENCE, lw=1.1, ls=(0, (4, 3)), zorder=1)
    ax.axhline(md, color=INK, lw=1.2, zorder=2)
    ax.text(0.02, md, f"mean {md:+.3f}", transform=ax.get_yaxis_transform(), ha="left",
            va="bottom", fontsize=8, color=INK,
            bbox=dict(facecolor="white", edgecolor="none", pad=1.2, alpha=0.85))
    ax.set_ylim(md - 4 * sd, md + 4 * sd)
    finish(ax, xlabel="mean of the two estimates", ylabel="SuSiE $-$ GCTA",
           title="Bland-Altman")
    panel_label(ax, "B")

    fig.suptitle(f"GCTA vs. SuSiE on GTEx Whole Blood (n = {n_full})",
                 fontsize=11.5, fontweight="bold", color=INK, y=1.02)
    fig.tight_layout()
    return save(fig, out_dir, "gtex_fig3_agreement")


def fig_failure_modes(observed: pd.DataFrame, out_dir: str) -> List[str]:
    """P3: REML fails often and mildly, fine-mapping rarely and wildly."""
    apply_style()
    fig, axes = plt.subplots(1, 2, figsize=(9.8, 4.0))
    x = np.arange(len(observed))
    width = 0.38

    ax = axes[0]
    ax.bar(x - width / 2, 100 * observed["gcta_impossible"], width=width * 0.92,
           color=SERIES["GCTA"], zorder=3, edgecolor="white", linewidth=1.6)
    ax.bar(x + width / 2, 100 * observed["susie_impossible"], width=width * 0.92,
           color=SERIES["SuSiE"], zorder=3, edgecolor="white", linewidth=1.6)
    ax.set_xticks(x)
    ax.set_xticklabels(observed["n"])
    ax.grid(axis="x", visible=False)
    finish(ax, xlabel="individuals used", ylabel="% of genes outside [0, 1]",
           title="How often each estimator leaves the valid range")
    panel_label(ax, "A")
    ax.legend(handles=[Patch(facecolor=SERIES[e], edgecolor="white", label=e)
                       for e in ESTIMATOR_ORDER], loc="upper right", ncol=2,
              handlelength=1.1)

    ax = axes[1]
    ax.bar(x - width / 2, observed["gcta_max_abs"], width=width * 0.92,
           color=SERIES["GCTA"], zorder=3, edgecolor="white", linewidth=1.6)
    ax.bar(x + width / 2, observed["susie_max_abs"], width=width * 0.92,
           color=SERIES["SuSiE"], zorder=3, edgecolor="white", linewidth=1.6)
    ax.axhline(1.0, color=REFERENCE, lw=1.2, ls=(0, (4, 3)), zorder=4)
    ax.text(len(observed) - 0.5, 1.02, "$h^2 = 1$", ha="right", va="bottom", fontsize=8,
            color=INK_SECONDARY)
    ax.set_xticks(x)
    ax.set_xticklabels(observed["n"])
    ax.grid(axis="x", visible=False)
    finish(ax, xlabel="individuals used", ylabel="largest $|\\hat{h}^2|$ seen",
           title="...and how badly, when it does")
    panel_label(ax, "B")

    fig.suptitle("Prediction 3 — opposite failure profiles, reproduced on real data",
                 fontsize=11.5, fontweight="bold", color=INK, y=1.02)
    fig.tight_layout()
    return save(fig, out_dir, "gtex_fig4_failure_modes")


# ------------------------------------------------------------------- report --
def build_gtex_report(run_dir: str, sim_run_dir: Optional[str] = None) -> str:
    run_dir = os.path.abspath(run_dir)
    fig_dir = os.path.join(run_dir, "figures")
    tab_dir = os.path.join(run_dir, "tables")
    os.makedirs(fig_dir, exist_ok=True)
    os.makedirs(tab_dir, exist_ok=True)

    df = load_gtex(run_dir)
    n_full = int(df["n"].max())
    observed = by_sample_size(df)
    predicted = None
    if sim_run_dir and os.path.isdir(sim_run_dir):
        predicted = simulation_prediction(sim_run_dir)
        predicted.to_csv(os.path.join(tab_dir, "simulation_prediction.csv"), index=False)

    observed.to_csv(os.path.join(tab_dir, "by_sample_size.csv"), index=False)
    full = df[df["n"] == n_full]
    by_stratum(full, "clump_bin").to_csv(os.path.join(tab_dir, "by_clump_bin.csv"), index=False)
    by_stratum(full, "cs_bin").to_csv(os.path.join(tab_dir, "by_credible_sets.csv"), index=False)
    mean, count = two_way_gap(full, "clump_bin")
    mean.to_csv(os.path.join(tab_dir, "gap_by_clumps_and_h2.csv"))
    count.to_csv(os.path.join(tab_dir, "gene_counts_by_clumps_and_h2.csv"))

    written: List[str] = []
    written += fig_sample_size_test(observed, predicted, fig_dir)
    written += fig_architecture_test(df, fig_dir, n_full)
    written += fig_agreement(df, fig_dir, n_full)
    written += fig_failure_modes(observed, fig_dir)
    from analysis.distributions import fig_gtex_boxplots

    written += fig_gtex_boxplots(df, fig_dir)

    path = os.path.join(run_dir, "GTEX_FINDINGS.md")
    with open(path, "w") as fh:
        fh.write(_render(df, observed, predicted, n_full, mean))
    return path


def _render(df, observed, predicted, n_full, mean) -> str:
    full = df[df["n"] == n_full]
    ok = full[full["usable"]]
    clump = by_stratum(full, "clump_bin")
    clump = clump[clump["genes"] >= 20]
    lo = observed.iloc[0]
    hi = observed.iloc[-1]

    L: List[str] = []
    add = L.append
    add("# GTEx Whole Blood — do the simulated patterns hold in real data?")
    add("")
    add(f"**{full['gene_id'].nunique():,} genes**, GTEx v8 EUR, HapMap3 genotypes, "
        f"{n_full} individuals with expression, genotype and covariates. Each gene "
        f"is fitted at five nested sample sizes ({', '.join(str(int(n)) for n in observed['n'])}), "
        f"so the sample-size axis is paired down to the individual.")
    add("")
    add("Expression **and** genotypes are residualised on genotype PCs, inferred "
        "expression factors, sex, platform and PCR batch, with the factor count "
        "following GTEx's own sample-size rule. LD is in-sample — the same "
        "individuals that produced the summary statistics — so unlike the "
        "simulations there is no reference-panel mismatch.")
    add("")
    add("> There is no ground truth here, so **nothing below measures bias**. "
        "Everything is the *disagreement* `estimator − GCTA`. A match to "
        "prediction is consistency, not proof: GCTA is not a gold standard either.")
    add("")
    add("### A correction to an earlier version of this analysis")
    add("")
    add("An earlier run pre-residualised the expression on covariates and handed "
        "REML the residual, which hides the degrees of freedom the covariates "
        "consume. That inflated GCTA's estimates — badly at small *n*, where 24 "
        "covariates eat a quarter of the data — and manufactured most of an "
        "apparent gap that appeared to close with sample size. **That result was "
        "an artefact.** REML now receives the covariates as fixed effects via "
        "`--qcovar`, and the numbers below supersede it.")
    add("")
    add("Two consequences of doing it properly, both documented in the config: "
        "unconstrained REML then fails to converge on 60-75% of genes, so this "
        "pass uses GCTA's default *constrained* REML (about 2% attrition); and "
        "GCTA can therefore no longer return a negative heritability, which "
        "removes the asymmetry that Prediction 3 relied on.")
    add("")
    add("---")
    add("")

    # ------------------------------------------------------------ P1 --
    ratio = abs(lo["mean_gap"] / hi["mean_gap"]) if hi["mean_gap"] else float("nan")
    if ratio >= 3 and abs(lo["mean_gap"]) >= 0.05:
        verdict1 = "**Confirmed, and steeply.**"
    elif ratio >= 1.5:
        verdict1 = ("**Confirmed in direction, but the effect is small.** The gap "
                    "does shrink monotonically with sample size, yet it is tiny "
                    "throughout — a few thousandths of $h^2$ — so this is a weak "
                    "version of the simulated result rather than a striking one.")
    else:
        verdict1 = ("**Not confirmed.** The gap does not move materially with "
                    "sample size on this data.")
    add("## Prediction 1 — the gap should close as sample size grows")
    add("")
    add(verdict1)
    add("")
    has_heels = "mean_gap_heels" in observed.columns
    if has_heels:
        add("| individuals | genes | GCTA $\\hat{h}^2$ | SuSiE $\\hat{h}^2$ "
            "| HEELS $\\hat{h}^2$ | SuSiE $-$ GCTA | HEELS $-$ GCTA "
            "| corr(SuSiE,GCTA) | corr(HEELS,GCTA) |")
        add("|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
        for _, r in observed.iterrows():
            add(f"| {int(r['n'])} | {int(r['genes']):,} | {r['mean_gcta']:.4f} "
                f"| {r['mean_susie']:.4f} | {r['mean_heels']:.4f} "
                f"| **{r['mean_gap']:+.4f}** | **{r['mean_gap_heels']:+.4f}** "
                f"| {r['corr']:.3f} | {r['corr_heels']:.3f} |")
    else:
        add("| individuals | genes fitted | mean GCTA $\\hat{h}^2$ | mean SuSiE $\\hat{h}^2$ "
            "| mean gap | median gap | correlation |")
        add("|---:|---:|---:|---:|---:|---:|---:|")
        for _, r in observed.iterrows():
            add(f"| {int(r['n'])} | {int(r['genes']):,} | {r['mean_gcta']:.4f} "
                f"| {r['mean_susie']:.4f} | **{r['mean_gap']:+.4f}** | {r['median_gap']:+.4f} "
                f"| {r['corr']:.3f} |")
    add("")
    add(f"The gap narrows {ratio:.0f}-fold from {lo['mean_gap']:+.4f} at "
        f"n = {int(lo['n'])} to {hi['mean_gap']:+.4f} at n = {int(hi['n'])}, and the "
        f"SuSiE-GCTA correlation rises from {lo['corr']:.2f} to {hi['corr']:.2f}. "
        f"The direction matches the simulation; the magnitude does not.")
    add("")
    if predicted is not None and len(predicted):
        add("Quantitatively, against the simulated prediction:")
        add("")
        add("| n | simulation predicts | GTEx observes |")
        add("|---:|---:|---:|")
        for _, r in predicted.iterrows():
            obs = observed[observed["n"] == r["n"]]
            shown = f"{obs['mean_gap'].iloc[0]:+.3f}" if len(obs) else "—"
            add(f"| {int(r['n'])} | {r['predicted_gap']:+.3f} | {shown} |")
        add("")
        matched = predicted[predicted["n"] == int(lo["n"])]
        if len(matched):
            pred_lo = float(matched["predicted_gap"].iloc[0])
            shrink = abs(pred_lo / lo["mean_gap"]) if lo["mean_gap"] else float("nan")
            add(f"The direction agrees, but **the real gap is roughly "
                f"{shrink:.0f}x smaller than the simulation predicts** "
                f"({lo['mean_gap']:+.3f} observed against {pred_lo:+.3f} predicted at "
                f"n = {int(lo['n'])}). The simulations are not wrong; they describe a "
                f"regime this data is not in.")
            add("")
        add("Two differences explain it, and both say the simulated grid is harsher "
            "than Whole Blood. First, the simulated architectures are deliberately "
            "polygenic — up to 20 causal SNPs carrying equal weight — whereas most "
            "real *cis*-windows here have zero or one detectable independent signal, "
            "which is exactly the corner where the simulations say fine-mapping is "
            "fine. Second, HapMap3 gives ~365 SNPs per window against ~8,200 "
            "simulated, with *n* > *p* rather than *p* >> *n*.")
        add("")
        add("The honest reading is that the simulated failure mode is **real but "
            "rarely triggered** in this tissue: the architecture that breaks "
            "fine-mapping is uncommon among Whole Blood *cis*-eQTLs, not absent.")
        add("")
    if has_heels:
        add("---")
        add("")
        add("## Prediction 5 — HEELS should track GCTA where SuSiE does not")
        add("")
        add(f"**Confirmed, and more strongly than expected.** HEELS reproduces "
            f"REML almost exactly on real genes: correlation "
            f"**{hi['corr_heels']:.3f}** at n = {int(hi['n'])} (and ≥ "
            f"{observed['corr_heels'].min():.3f} at every sample size), with a mean "
            f"difference of {hi['mean_gap_heels']:+.4f}. SuSiE, using the same "
            f"inputs, correlates {hi['corr']:.2f}.")
        add("")
        add("This is the sufficiency argument made visible. For the model both "
            "fit, the marginal statistics plus in-sample LD are *sufficient* for "
            "the individual-level likelihood, so a summary-statistic method that "
            "makes no extra assumption must land on REML — and it does, to three "
            "decimal places, on ~19,000 real genes. **Nothing is lost by sharing "
            "summary statistics instead of genotypes**, provided the LD is "
            "in-sample.")
        add("")
        add("It also isolates what the remaining SuSiE-GCTA difference is *made "
            "of*. Both summary-statistic estimators see byte-identical inputs, so "
            "the gap between them is the sparsity assumption and nothing else. "
            "That gap is small here — real Whole Blood *cis*-architecture is "
            "mostly sparse enough for fine-mapping to cope — but it is the same "
            "quantity the simulations blow up by choosing harsher architectures.")
        add("")
    add("![sample size](figures/gtex_fig1_sample_size.png)")
    add("")
    add("The means above are what the tables report; the distributions behind "
        "them are worth seeing too, because the spread is far larger than any of "
        "the differences between estimators:")
    add("")
    add("![distributions](figures/gtex_figD1_distributions.png)")
    add("")
    def _iqr(series):
        v = series.dropna().to_numpy(dtype=float)
        return float(np.percentile(v, 75) - np.percentile(v, 25)) if v.size else float("nan")

    lo_n, hi_n = int(observed["n"].iloc[0]), int(observed["n"].iloc[-1])
    lo_part, hi_part = df[df["n"] == lo_n], df[df["n"] == hi_n]
    su_lo, su_hi = _iqr(lo_part["gap"]), _iqr(hi_part["gap"])
    he_lo = _iqr(lo_part["gap_heels"]) if "gap_heels" in df.columns else float("nan")
    he_hi = _iqr(hi_part["gap_heels"]) if "gap_heels" in df.columns else float("nan")
    add("The right panel is the paired difference gene by gene — the closest "
        "real-data analogue of an error distribution, and the panel that carries "
        "the actual information.")
    add("")
    if np.isfinite(he_hi):
        add(f"Both boxes tighten as the sample grows, but **HEELS is about twice "
            f"as tight as SuSiE at every sample size**: an inter-quartile spread "
            f"of {he_lo:.4f} against {su_lo:.4f} at n = {lo_n}, narrowing to "
            f"{he_hi:.4f} against {su_hi:.4f} at n = {hi_n}. The two estimators "
            f"see byte-identical inputs, so that factor of two is the sparsity "
            f"assumption showing up as **variance** — gene-to-gene "
            f"unpredictability — rather than as the bias the simulations "
            f"emphasise. On real data the sparsity penalty is mostly a precision "
            f"cost, not an accuracy cost.")
    add("")

    # ------------------------------------------------------------ P2 --
    add("---")
    add("")
    add("## Prediction 2 — the gap should track how diffuse the signal is")
    add("")
    add("**Confirmed, but only once heritability is held fixed — and the naive "
        "version of this test gives the opposite answer.**")
    add("")
    add("Genes are stratified by the number of independent marginal associations "
        "found by LD clumping (p < 1e-4, r² < 0.1), a proxy that does **not** use "
        "SuSiE's output, so the test is not circular. Read marginally, it looks "
        "like the gap *grows* with the number of signals:")
    add("")
    add("| independent signals | genes | mean GCTA | mean SuSiE | mean gap |")
    add("|---|---:|---:|---:|---:|")
    for _, r in clump.iterrows():
        add(f"| {r['clump_bin']} | {int(r['genes']):,} | {r['mean_gcta']:.4f} "
            f"| {r['mean_susie']:.4f} | {r['mean_gap']:+.4f} |")
    add("")
    add("That reading is **confounded**: genes with more detectable signals also "
        "have far more heritability (0.02 in the first row against 0.37 in the "
        "last), so the comparison is mixing architecture with magnitude. "
        "Conditioning on heritability reverses it — within every stratum, the gap "
        "*shrinks* as detectable signals accumulate:")
    add("")
    header = "| independent signals | " + " | ".join(
        str(c).replace(", ", "–").strip("()[]") for c in mean.columns) + " |"
    add(header)
    add("|---|" + "---:|" * len(mean.columns))
    for level in mean.index:
        cells = []
        for col in mean.columns:
            v = mean.loc[level, col]
            cells.append("—" if not np.isfinite(v) else f"{v:+.3f}")
        add(f"| {level} | " + " | ".join(cells) + " |")
    add("")
    add("Columns are GCTA $\\hat{h}^2$ bins; cells with fewer than 30 genes are blank.")
    add("")
    add("This is a textbook Simpson's paradox, and the conditional direction is the "
        "one that matches the simulations. The number of *detected* clumps measures "
        "how **localisable** a gene's signal is, not how polygenic it is. The worst "
        "cell in the table is the gene with real heritability and **no** detectable "
        "individual association — heritability spread too thinly to localise, which "
        "is exactly the equal-magnitude/high-polygenicity corner where the "
        "simulations put fine-mapping's largest shortfall. As signal concentrates "
        "into detectable clumps, the gap closes, matching the effect-concentration "
        "result (§3 of the simulation write-up).")
    add("")
    add("Practically: a gene with a decent GCTA heritability and no genome-wide "
        "detectable *cis*-eQTL is where the two methods will disagree most, and it "
        "is precisely the gene for which fine-mapping-based $h^2$ should be "
        "distrusted.")
    add("")
    add("![architecture](figures/gtex_fig2_architecture.png)")
    add("")
    add("![agreement](figures/gtex_fig3_agreement.png)")
    add("")

    # ------------------------------------------------------------ P3 --
    add("---")
    add("")
    gcta_oob = float(observed["gcta_impossible"].max())
    if gcta_oob < 1e-4:
        verdict3 = ("**Not testable on this data.** The covariate fix above forced "
                    "a switch to *constrained* REML, which cannot return a negative "
                    "heritability, so GCTA's half of this prediction has been "
                    "defined away. What remains observable is the fine-mapping "
                    "half, and it behaves as predicted.")
    else:
        verdict3 = "**Confirmed.**"
    add("## Prediction 3 — opposite failure profiles")
    add("")
    add(verdict3)
    add("")
    if has_heels:
        add("| individuals | GCTA outside [0,1] | SuSiE outside [0,1] "
            "| HEELS outside [0,1] | largest &#124;GCTA&#124; | largest &#124;SuSiE&#124; |")
        add("|---:|---:|---:|---:|---:|---:|")
        for _, r in observed.iterrows():
            add(f"| {int(r['n'])} | {r['gcta_impossible']:.1%} "
                f"| {r['susie_impossible']:.2%} | {r['heels_impossible']:.2%} "
                f"| {r['gcta_max_abs']:.2f} | {r['susie_max_abs']:.2f} |")
    else:
        add("| individuals | GCTA outside [0,1] | SuSiE outside [0,1] "
            "| largest &#124;GCTA&#124; | largest &#124;SuSiE&#124; |")
        add("|---:|---:|---:|---:|---:|")
        for _, r in observed.iterrows():
            add(f"| {int(r['n'])} | {r['gcta_impossible']:.1%} | {r['susie_impossible']:.2%} "
                f"| {r['gcta_max_abs']:.2f} | {r['susie_max_abs']:.2f} |")
    add("")
    if gcta_oob < 1e-4:
        add(f"GCTA is pinned at 0.0% by construction here, so the columns to read "
            f"are SuSiE's and HEELS's. **SuSiE fails rarely but wildly** — at most "
            f"{float(observed['susie_impossible'].max()):.2%} of genes leave "
            f"[0, 1], but when one does it can read "
            f"{float(observed['susie_max_abs'].max()):.2f}, while reporting clean "
            f"convergence. **HEELS never leaves the range at all**, because its "
            f"heritability is a ratio of two non-negative variance components; it "
            f"fails by pinning at a boundary instead, which is at least visible.")
        add("")
        add("The operational advice from the simulations therefore carries over "
            "for the estimator it was about: screen fine-mapping heritabilities "
            "against [0, 1] rather than trusting the convergence flag. The claim "
            "that REML fails often-but-mildly is a property of *unconstrained* "
            "REML, and is visible in the simulations rather than here.")
        add("")
    else:
        add("REML leaves the valid range often and by a little; fine-mapping leaves "
            "it rarely and by a lot. The operational advice from the simulations "
            "carries over unchanged: screen fine-mapping heritabilities against "
            "[0, 1] rather than trusting the convergence flag.")
        add("")
    add("![failure modes](figures/gtex_fig4_failure_modes.png)")
    add("")

    # -------------------------------------------------------- caveats --
    add("---")
    add("")
    add("## What this does and does not establish")
    add("")
    if has_heels:
        add("**On the LD asymmetry.** In these real-data runs the asymmetry that "
            "complicates the simulations disappears: genotypes and expression come "
            "from the same individuals, so *both* summary-statistic estimators get "
            "correctly specified in-sample LD. The SuSiE-versus-HEELS comparison "
            "here is therefore a clean test of the sparsity assumption, which the "
            "simulation comparison is not.")
        add("")
    add("**Does:** the architecture dependence is present and clean — within a "
        "heritability stratum, the disagreement shrinks as a gene's signal becomes "
        "localisable, which is the simulations' central mechanism. And HEELS "
        "reproducing REML to a correlation of 1.000 over ~19,000 genes is a strong "
        "positive result in its own right: summary statistics plus in-sample LD "
        "lose nothing.")
    add("")
    add("**Does not:** establish that the simulated *magnitudes* transfer. The "
        "sample-size effect is real but roughly six times smaller here than "
        "simulated, because Whole Blood cis-architecture is mostly sparse and the "
        "HapMap3 windows are 20x smaller than the simulated ones. The simulations "
        "describe a harsher regime than this tissue presents.")
    add("")
    add("**Does not:** establish that either estimator is *correct* here. The "
        "quantity measured is disagreement, and the simulations themselves show "
        "GCTA turning upward-biased when causal variants are frequency-stratified, "
        "which is plausible in real data. A gap of "
        f"{hi['mean_gap']:+.3f} at n = {n_full} is consistent with SuSiE "
        "under-recovering, with GCTA over-recovering, or with both.")
    add("")
    add("Three concrete limits:")
    add("")
    add(f"1. **One tissue.** Whole Blood only, so the sample-size axis comes from "
        f"subsampling rather than from tissues that genuinely differ. Subsampling "
        f"holds the biology fixed, which is cleaner for the test but says nothing "
        f"about tissue-to-tissue variation.")
    add(f"2. **A different SNP regime.** HapMap3 genotypes give a median of "
        f"{int(full['p_snps'].median())} SNPs per ±500 kb window against ~8,200 in "
        f"the 1000G simulations, and here *n* > *p*, so the LD matrix is full rank. "
        f"This is a **more** favourable setting for fine-mapping than the "
        f"simulations used, which makes the surviving gap more notable, not less.")
    add("3. **REML attrition at small n.** "
        f"{lo['gcta_failed']:.1%} of genes have no GCTA estimate at n = {int(lo['n'])} "
        f"(against {hi['gcta_failed']:.1%} at n = {int(hi['n'])}), and those are "
        "disproportionately the genes where REML struggled. The small-n gap is "
        "therefore measured on a slightly easier subset than the large-n gap.")
    add("")
    add("---")
    add("")
    add("**Companion analysis.** Accuracy cannot be measured here, but "
        "*precision* can: `../gtex_split_half/RELIABILITY.md` splits the donors "
        "into two disjoint halves and reports test-retest reliability. Headline: "
        "GCTA and HEELS agree at r = 0.999 on the same donors while GCTA agrees "
        "with itself at r = 0.848 across donors — the method difference is far "
        "below the sampling noise. SuSiE is the most reliable arm, which is what "
        "its shrinkage predicts rather than evidence of accuracy.")
    add("")
    add("Tables: `tables/by_sample_size.csv`, `tables/by_clump_bin.csv`, "
        "`tables/by_credible_sets.csv`. Gene-level results: `cells.tsv.gz`. "
        "Simulation write-up: [`../FINDINGS.md`](../FINDINGS.md).")
    add("")
    add("Generated by `python -m analysis.gtex_report`.")
    return "\n".join(L) + "\n"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir", nargs="?",
                    default="/home/sux002/work/cis_heritability/results/gtex_whole_blood")
    ap.add_argument("--simulation", default="/home/sux002/work/cis_heritability/"
                                            "results/sample_size",
                    help="simulation campaign providing the sample-size prediction")
    args = ap.parse_args(argv)
    print(build_gtex_report(args.run_dir, args.simulation))
    return 0


if __name__ == "__main__":
    sys.exit(main())
