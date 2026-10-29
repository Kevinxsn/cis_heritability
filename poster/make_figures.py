"""Poster figures for IGVF 2026.

These are *not* the paper figures.  A poster is read from ~1.5 m, so every panel
here is redrawn with fewer facets, thicker strokes and much larger type than the
manuscript versions in ../results/*/figures/.  Numbers are read from the result
tables at build time, exactly as the manuscript figures are, so the poster
cannot drift away from the analysis.

Estimator naming: the fine-mapping arm is called SuSiE in the pipeline and
CHEWS (Cis-Heritability Estimation With SuSiE) on the poster; the rename happens
in DISPLAY only.

    python make_figures.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent / "results"
OUT = HERE / "figures"
OUT.mkdir(parents=True, exist_ok=True)

# --- palette: unchanged from analysis/style.py (validated slots 1-3) --------
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
SERIES = {"GCTA": BLUE, "CHEWS": ORANGE, "HEELS": AQUA}
COL = {"GCTA": "gcta", "CHEWS": "susie", "HEELS": "heels"}
ORDER = ["GCTA", "CHEWS", "HEELS"]

INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#8a8880"
GRID = "#e2e1dd"

# Poster type scale: base 20 pt in figure units, rendered at 300 dpi onto a
# ~11 in column, i.e. roughly 30 pt on the printed sheet.
mpl.rcParams.update({
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
    # Liberation Sans is metrically identical to Arial, which the poster body
    # text uses - so figure labels and prose share one typeface.
    "font.family": "Liberation Sans",
    "font.size": 22,
    "axes.titlesize": 24,
    "axes.titleweight": "bold",
    "axes.labelsize": 23,
    "axes.labelcolor": INK,
    "axes.edgecolor": INK2,
    "axes.linewidth": 1.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.color": GRID,
    "grid.linewidth": 1.4,
    "xtick.color": INK2,
    "ytick.color": INK2,
    "xtick.labelsize": 21,
    "ytick.labelsize": 21,
    "xtick.major.width": 1.6,
    "ytick.major.width": 1.6,
    "legend.frameon": False,
    "legend.fontsize": 22,
    "lines.linewidth": 4.0,
    "lines.markersize": 12,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "pdf.fonttype": 42,
})


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"{name}.{ext}")
    plt.close(fig)
    print(f"  wrote figures/{name}.png")


def legend_handles(names=ORDER):
    return [Line2D([0], [0], color=SERIES[k], lw=4.5, marker="o", ms=11, label=k)
            for k in names]


TITLE_PAD = 30

# A figure drawn W inches wide but placed in an 11.07 in column has all of its
# type multiplied by 11.07/W.  Drawing at the column width makes every point
# size in rcParams the size it will actually print at, which is the only way
# the >=20 pt floor for annotations is meaningful.
COLW_IN = 11.0
FULLW_IN = 33.0


def region_bootstrap_ci(df, value_fn, by, n_boot=2000, seed=20260908):
    """95% CI resampling REGIONS, not rows.

    The same cis-region supplies a cell at every causal count and every h2, so
    rows are clustered within region and an ordinary SE understates the
    uncertainty.  Resampling whole regions with replacement respects that.
    """
    rng = np.random.default_rng(seed)
    regions = df["region_id"].unique()
    idx_by_region = {r: g.index.to_numpy() for r, g in df.groupby("region_id")}
    keys = sorted(df[by].unique())
    draws = {k: [] for k in keys}
    for _ in range(n_boot):
        pick = rng.choice(regions, size=len(regions), replace=True)
        rows = np.concatenate([idx_by_region[r] for r in pick])
        boot = df.loc[rows]
        g = value_fn(boot).groupby(boot[by]).mean() if False else None
        vals = value_fn(boot)
        agg = vals.groupby(boot[by]).mean()
        for k in keys:
            draws[k].append(agg.get(k, np.nan))
    lo = {k: float(np.nanpercentile(v, 2.5)) for k, v in draws.items()}
    hi = {k: float(np.nanpercentile(v, 97.5)) for k, v in draws.items()}
    return keys, lo, hi


def pairwise_win(cells, by):
    """% of paired fits where CHEWS lands closer to the truth than GCTA.

    Head-to-head, not the pipeline's three-way `susie_win_rate`: the poster
    claims "beats GCTA", and that is the number which supports it (it is also
    the 76% quoted in the submitted abstract at one causal SNP).
    """
    d = cells.dropna(subset=["susie_h2", "gcta_h2"]).copy()
    d["_s"] = (d["susie_h2"] - d["h2_true"]).abs()
    d["_g"] = (d["gcta_h2"] - d["h2_true"]).abs()
    g = d.groupby(by).apply(lambda x: pd.Series({
        "rate": float((x["_s"] < x["_g"]).mean()), "n": len(x)}))
    g["se"] = np.sqrt(g["rate"] * (1 - g["rate"]) / g["n"])
    return g.reset_index()


def panel_tag(ax, letter, dx=-0.145, dy=1.26):
    """Tag sits *above* the title, so the two never collide."""
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=30,
            fontweight="bold", color=INK, ha="left", va="top")


def panel_tag_in(ax, letter):
    ax.text(0.028, 0.965, letter, transform=ax.transAxes, fontsize=26,
            fontweight="bold", color=INK, ha="left", va="top", zorder=6,
            path_effects=[pe.withStroke(linewidth=4.0, foreground="white")])


def titled(ax, text, size=23.5):
    ax.set_title(text, pad=TITLE_PAD, fontsize=size)


def cell_text(ax, j, i, label, norm_v):
    """Dark ink on pale cells, white on saturated ones - a white-with-stroke
    label is unreadable on the near-white middle of a diverging ramp."""
    dark = abs(norm_v) > 0.45
    ax.text(j, i, label, ha="center", va="center", fontsize=19,
            fontweight="bold", color="white" if dark else INK,
            path_effects=[pe.withStroke(linewidth=2.6,
                                        foreground="#00000055" if dark else "#ffffffcc")])


def _ticks(ax):
    """Log-n axis; 800 is drawn but not labelled - it collides with 1000."""
    ax.set_xscale("log")
    ax.set_xticks([100, 300, 500, 800, 1000])
    ax.set_xticklabels(["100", "300", "500", "", "1000"])
    ax.minorticks_off()


# ---------------------------------------------------------------------------
# P1. Polygenicity is the switch, and it is specific to fine-mapping
# ---------------------------------------------------------------------------
def fig_polygenicity():
    cells = pd.read_csv(RESULTS / "baseline_n300" / "cells.tsv.gz", sep="\t")

    fig, axes = plt.subplots(1, 2, figsize=(COLW_IN, 4.7))

    # A: bias vs causal SNPs, pooled over the h2 grid
    ax = axes[0]
    for name in ORDER:
        c = COL[name]
        d = cells.dropna(subset=[f"{c}_h2"]).copy()
        d["err"] = d[f"{c}_h2"] - d["h2_true"]
        g = d.groupby("ncausal")["err"].agg(["mean", "sem"]).reset_index()
        ax.errorbar(g["ncausal"], g["mean"], yerr=1.96 * g["sem"],
                    color=SERIES[name], marker="o", capsize=6, capthick=2.4,
                    elinewidth=2.4, label=name, zorder=3)
    ax.axhline(0, color=MUTED, ls="--", lw=2.2, zorder=1)
    ax.set_xscale("log")
    ax.set_xticks([1, 3, 5, 10, 20])
    ax.set_xticklabels(["1", "3", "5", "10", "20"])
    ax.minorticks_off()
    ax.set_xlabel("causal SNPs in the $\\it{cis}$-window")
    ax.set_ylabel("bias  ($\\hat{h}^2 - h^2$)")
    titled(ax, "Bias grows with polygenicity")
    ax.legend(handles=legend_handles(), loc="lower left")
    panel_tag(ax, "A")

    # B: how often CHEWS is the closer estimate
    ax = axes[1]
    g = pairwise_win(cells, "ncausal")
    bars = ax.bar(np.arange(len(g)), 100 * g["rate"],
                  yerr=196 * g["se"], color=ORANGE, width=0.66,
                  error_kw=dict(ecolor=INK2, elinewidth=2.4, capsize=7, capthick=2.2),
                  zorder=3)
    ax.axhline(50, color=MUTED, ls="--", lw=2.2, zorder=1)
    ax.text(len(g) - 0.42, 51.5, "coin flip", color=INK2, fontsize=17,
            ha="right", va="bottom", style="italic")
    for i, (b, r) in enumerate(zip(bars, g["rate"])):
        ax.text(b.get_x() + b.get_width() / 2, 100 * r + 5.5, f"{100*r:.0f}%",
                ha="center", va="bottom", fontsize=18, fontweight="bold", color=INK)
    ax.set_xticks(np.arange(len(g)))
    ax.set_xticklabels([str(int(v)) for v in g["ncausal"]])
    ax.set_xlabel("causal SNPs in the $\\it{cis}$-window")
    ax.set_ylabel("loci where CHEWS beats\nGCTA (%)")
    ax.set_ylim(0, 88)
    titled(ax, "CHEWS beats GCTA when signal is sparse")
    panel_tag(ax, "B")

    fig.tight_layout(w_pad=1.6)
    save(fig, "P1_polygenicity")


# ---------------------------------------------------------------------------
# P2. Where each estimator wins, across the whole simulated grid
# ---------------------------------------------------------------------------
# One encoding for "how often CHEWS is closer to truth", shared by P2 and the
# round-4 panels (P3 B, P6 B): a given colour means the same win rate in every
# section.  Centred on 50%.
WIN_CMAP = "RdBu_r"
WIN_VMIN, WIN_VMAX = 15, 85


def fig_landscape():
    cells = pd.read_csv(RESULTS / "baseline_n300" / "cells.tsv.gz", sep="\t")
    g = pairwise_win(cells, ["h2_true", "ncausal"])
    piv = g.pivot(index="h2_true", columns="ncausal", values="rate")
    piv = piv.sort_index(ascending=False)

    fig, ax = plt.subplots(figsize=(COLW_IN, 5.2))
    im = ax.imshow(100 * piv.values, cmap=WIN_CMAP, vmin=WIN_VMIN, vmax=WIN_VMAX, aspect="auto")
    ax.set_xticks(range(piv.shape[1]))
    ax.set_xticklabels([str(int(c)) for c in piv.columns])
    ax.set_yticks(range(piv.shape[0]))
    ax.set_yticklabels([f"{v:g}" for v in piv.index])
    ax.set_xlabel("causal SNPs")
    ax.set_ylabel("true $h^2$")
    ax.grid(False)
    for i in range(piv.shape[0]):
        for j in range(piv.shape[1]):
            v = 100 * piv.values[i, j]
            cell_text(ax, j, i, f"{v:.0f}", (v - 50) / 35)
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    cb.set_label("% of loci, CHEWS closer to truth", fontsize=18)
    cb.ax.tick_params(labelsize=16)
    titled(ax, "Fraction of simulations with lower\nCHEWS absolute error", size=21)
    fig.tight_layout()
    save(fig, "P2_win_landscape")


def _panel_dispersion(ax, m):
    """§5 panel A, shared by the round-3 and round-4 P3 figures so it is drawn
    by the same code in both (PI round 4: "keep as is")."""
    sd = (m.groupby(["estimator", "n"])
            .apply(lambda d: np.average(d["sd_h2_hat"], weights=d["n_cells"]))
            .rename("sd").reset_index())
    for name in ORDER:
        d = sd[sd["estimator"] == name].sort_values("n")
        ax.plot(d["n"], d["sd"], color=SERIES[name], marker="o", label=name, zorder=3)
    _ticks(ax)
    ax.set_xlabel("individuals ($n$)")
    ax.set_ylabel("SD of $\\hat{h}^2$ across loci")
    titled(ax, "Dispersion of $\\hat{h}^2$", size=21)
    for name, dy in (("HEELS", 10), ("GCTA", -14), ("CHEWS", -12)):
        v = sd[(sd["estimator"] == name) & (sd["n"] == 100)]["sd"].iloc[0]
        ax.annotate(name, xy=(100, v), xytext=(12, dy), textcoords="offset points",
                    color=SERIES[name], fontsize=19, fontweight="bold",
                    ha="left", va="center", zorder=6,
                    path_effects=[pe.withStroke(linewidth=3.5, foreground="white")])
    ax.set_xlim(left=88)
    panel_tag(ax, "A", dx=-0.30, dy=1.30)


def _panel_mae_diff(fig, ax, m):
    """Round-3 §5 panel B (MAE difference).  Also used by the round-4 figure
    to reproduce the round-3 layout exactly before swapping the panel."""
    piv = (m.pivot_table(index="ncausal", columns=["estimator", "n"],
                         values="mae", aggfunc=lambda v: np.average(v)))
    diff = piv["CHEWS"] - piv["GCTA"]
    diff = diff.sort_index(ascending=False)
    v = diff.values.astype(float)
    lim = np.nanmax(np.abs(v))
    im = ax.imshow(v, cmap="RdBu_r", vmin=-lim, vmax=lim, aspect="auto")
    ax.set_xticks(range(diff.shape[1]))
    ax.set_xticklabels([str(int(c)) for c in diff.columns])
    ax.set_yticks(range(diff.shape[0]))
    ax.set_yticklabels([str(int(i)) for i in diff.index], fontsize=22, fontweight="bold")
    ax.set_xlabel("individuals ($n$)")
    ax.set_ylabel("causal SNPs")
    ax.grid(False)
    for i in range(v.shape[0]):
        for j in range(v.shape[1]):
            lab = f"{v[i, j]:+.02f}".replace("+0.00", "0.00")
            dark = abs(v[i, j] / lim) > 0.45
            ax.text(j, i, lab, ha="center", va="center", fontsize=15.5,
                    fontweight="bold", color="white" if dark else INK,
                    path_effects=[pe.withStroke(
                        linewidth=2.4,
                        foreground="#00000055" if dark else "#ffffffcc")])
    cb = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.03)
    cb.set_label("MAE$_{CHEWS}$ $-$ MAE$_{GCTA}$\n(negative favors CHEWS)", fontsize=17)
    cb.ax.tick_params(labelsize=16)
    titled(ax, "Mean absolute error difference", size=21)
    panel_tag(ax, "B", dx=-0.26, dy=1.34)
    return cb


# ---------------------------------------------------------------------------
# P3. Sample size: bias falls, but the advantage falls with it
# ---------------------------------------------------------------------------
def fig_sample_size():
    """Precision directly, and the n-by-architecture interaction.

    Replaces the old bias+MAE pair.  Bias against n is already carried by the
    text; what was missing is (a) an estimator-level SD, which is what
    "precision" actually means and what the title claims, and (b) the
    interaction that the decision rule depends on - small n AND concentrated
    effects, not small n alone.
    """
    m = pd.read_csv(RESULTS / "sample_size" / "tables" / "metrics_by_n.csv")
    m["estimator"] = m["estimator"].replace({"SuSiE": "CHEWS"})

    fig, axes = plt.subplots(1, 2, figsize=(COLW_IN, 4.9),
                             gridspec_kw=dict(width_ratios=[1.0, 1.30]))
    _panel_dispersion(axes[0], m)

    # B: where CHEWS actually lowers error - n x causal count
    _panel_mae_diff(fig, axes[1], m)

    fig.tight_layout(w_pad=1.6)
    save(fig, "P3_sample_size")


# ---------------------------------------------------------------------------
# P4. HEELS separates the data requirement from the modelling assumption
# ---------------------------------------------------------------------------
def fig_sufficiency():
    cells = pd.read_csv(RESULTS / "baseline_n300" / "cells.tsv.gz", sep="\t")
    eff = pd.read_csv(RESULTS / "effect_size" / "cells.tsv.gz", sep="\t")

    fig, axes = plt.subplots(1, 2, figsize=(COLW_IN, 4.7))

    # A: bias by effect concentration - the second axis of the same mechanism
    ax = axes[0]
    labels = {"fixed": "equal\nmagnitude", "normal": "normal", "pareto": "heavy\ntailed"}
    keys = ["fixed", "normal", "pareto"]
    xs = np.arange(len(keys))
    w = 0.26
    for k, name in enumerate(ORDER):
        c = COL[name]
        vals, errs = [], []
        for key in keys:
            d = eff[(eff["effect_dist"] == key)].dropna(subset=[f"{c}_h2"])
            e = d[f"{c}_h2"] - d["h2_true"]
            vals.append(e.mean())
            errs.append(1.96 * e.sem())
        ax.bar(xs + (k - 1) * w, vals, width=w, yerr=errs, color=SERIES[name],
               label=name, zorder=3,
               error_kw=dict(ecolor=INK2, elinewidth=2.0, capsize=5, capthick=1.8))
    ax.axhline(0, color=INK2, lw=1.8, zorder=2)
    ax.set_xticks(xs)
    ax.set_xticklabels([labels[k] for k in keys])
    ax.set_xlabel("how the fixed total effect is spread")
    ax.set_ylabel("bias  ($\\hat{h}^2 - h^2$)")
    titled(ax, "Concentration drives the gap")
    ax.legend(handles=legend_handles(), loc="lower right", ncol=1)
    panel_tag(ax, "A")

    # B: causal-MAF stratification - the one place GCTA breaks
    ax = axes[1]
    maf = pd.read_csv(RESULTS / "maf_arch" / "cells.tsv.gz", sep="\t")
    mlab = {"common": "common causal\n(MAF 0.05-0.5)", "rare": "rare causal\n(MAF 0.01-0.05)"}
    mkeys = [k for k in ["common", "rare"] if (maf["maf_bin"] == k).any()]
    xs = np.arange(len(mkeys))
    for k, name in enumerate(ORDER):
        c = COL[name]
        vals, errs = [], []
        for key in mkeys:
            d = maf[maf["maf_bin"] == key].dropna(subset=[f"{c}_h2"])
            e = d[f"{c}_h2"] - d["h2_true"]
            vals.append(e.mean())
            errs.append(1.96 * e.sem())
        ax.bar(xs + (k - 1) * w, vals, width=w, yerr=errs, color=SERIES[name],
               label=name, zorder=3,
               error_kw=dict(ecolor=INK2, elinewidth=2.0, capsize=5, capthick=1.8))
    ax.axhline(0, color=INK2, lw=1.8, zorder=2)
    ax.set_xticks(xs)
    ax.set_xticklabels([mlab[k] for k in mkeys])
    ax.set_xlabel("frequency of the causal variants")
    ax.set_ylabel("bias  ($\\hat{h}^2 - h^2$)")
    titled(ax, "The one place GCTA breaks")
    panel_tag(ax, "B")

    fig.tight_layout(w_pad=1.6)
    save(fig, "P4_architecture")


# ---------------------------------------------------------------------------
# P5. GTEx Whole Blood - the simulated pattern, in real data
# ---------------------------------------------------------------------------
def fig_gtex():
    t = pd.read_csv(RESULTS / "gtex_whole_blood" / "tables" / "by_sample_size.csv")

    fig, axes = plt.subplots(1, 2, figsize=(COLW_IN, 4.7))

    ax = axes[0]
    ax.errorbar(t["n"], t["mean_gap"], yerr=1.96 * t["gap_mcse"], color=ORANGE,
                marker="o", capsize=6, capthick=2.4, elinewidth=2.4,
                label="CHEWS $-$ GCTA", zorder=3)
    ax.errorbar(t["n"], t["mean_gap_heels"], yerr=1.96 * t["gap_heels_mcse"],
                color=AQUA, marker="o", capsize=6, capthick=2.4, elinewidth=2.4,
                label="HEELS $-$ GCTA", zorder=3)
    ax.axhline(0, color=MUTED, ls="--", lw=2.2, zorder=1)
    ax.set_xticks(t["n"])
    ax.set_xlabel("GTEx donors analyzed")
    ax.set_ylabel("mean difference from GCTA")
    titled(ax, "Gap closes with sample size")
    ax.legend(loc="lower right")
    panel_tag(ax, "A")

    ax = axes[1]
    ax.plot(t["n"], t["corr_heels"], color=AQUA, marker="o", label="HEELS vs GCTA", zorder=3)
    ax.plot(t["n"], t["corr"], color=ORANGE, marker="o", label="CHEWS vs GCTA", zorder=3)
    ax.set_xticks(t["n"])
    ax.set_ylim(0.75, 1.02)
    ax.set_xlabel("GTEx donors analyzed")
    ax.set_ylabel("correlation with GCTA $\\hat{h}^2$")
    titled(ax, "HEELS reproduces REML exactly")
    ax.annotate(f"$r$ = {t['corr_heels'].iloc[-1]:.3f}",
                xy=(t["n"].iloc[-1], t["corr_heels"].iloc[-1]),
                xytext=(-14, -46), textcoords="offset points", fontsize=19,
                fontweight="bold", color=AQUA, ha="right")
    ax.legend(loc="lower right")
    panel_tag(ax, "B")

    fig.tight_layout(w_pad=1.6)
    save(fig, "P5_gtex_sample_size")


# ---------------------------------------------------------------------------
# P6. GTEx: the disagreement is architectural, once h2 is held fixed
# ---------------------------------------------------------------------------
def fig_gtex_architecture():
    gap = pd.read_csv(RESULTS / "gtex_whole_blood" / "tables" / "gap_by_clumps_and_h2.csv")
    gap = gap.set_index("clump_bin")
    vals = gap.values.astype(float)

    fig, ax = plt.subplots(figsize=(COLW_IN, 5.2))
    cmap = plt.get_cmap("RdBu").copy()
    cmap.set_bad("#eceae5")          # NaN must not read as "gap = 0"
    im = ax.imshow(np.ma.masked_invalid(vals), cmap=cmap,
                   vmin=-0.055, vmax=0.055, aspect="auto")
    ax.set_xticks(range(gap.shape[1]))
    ax.set_xticklabels(["0–.05", ".05–.1", ".1–.2", ".2–.4", ".4–1"])
    ax.set_yticks(range(gap.shape[0]))
    ax.set_yticklabels(list(gap.index))
    ax.set_xlabel("GCTA $\\hat{h}^2$ stratum")
    ax.set_ylabel("independent $\\it{cis}$-eQTL signals\n(LD clumping, $p<10^{-4}$)")
    ax.grid(False)
    for i in range(vals.shape[0]):
        for j in range(vals.shape[1]):
            v = vals[i, j]
            if np.isnan(v):
                ax.text(j, i, "n/a", ha="center", va="center", color=MUTED,
                        fontsize=16, style="italic")
            else:
                cell_text(ax, j, i, f"{v:+.3f}", v / 0.055)
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    cb.set_label("CHEWS $-$ GCTA", fontsize=18)
    cb.ax.tick_params(labelsize=16)
    titled(ax, "Worst gap: real $h^2$, no clear signal")
    fig.tight_layout()
    save(fig, "P6_gtex_architecture")


# ---------------------------------------------------------------------------
# P7. Test-retest reliability on disjoint donor halves
# ---------------------------------------------------------------------------
def fig_reliability():
    rel = pd.read_csv(RESULTS / "gtex_split_half" / "tables" / "reliability_signal_genes.csv")
    cross = pd.read_csv(RESULTS / "gtex_split_half" / "tables" / "cross_estimator_agreement.csv")

    fig, axes = plt.subplots(1, 2, figsize=(COLW_IN, 4.7))

    ax = axes[0]
    rel = rel.copy()
    rel["label"] = rel["estimator"].replace({"SuSiE": "CHEWS"})
    rel = rel.set_index("label").loc[[e for e in ORDER if e in set(rel["label"])]]
    bars = ax.bar(np.arange(len(rel)), rel["pearson"],
                  color=[SERIES[k] for k in rel.index], width=0.62, zorder=3)
    for b, v in zip(bars, rel["pearson"]):
        ax.text(b.get_x() + b.get_width() / 2, v + 0.012, f"{v:.3f}", ha="center",
                va="bottom", fontsize=20, fontweight="bold", color=INK)
    ax.set_xticks(np.arange(len(rel)))
    ax.set_xticklabels(list(rel.index))
    ax.set_ylim(0.6, 0.99)
    ax.set_ylabel("half A vs half B,\nPearson $r$")
    titled(ax, "Test-retest on disjoint halves")
    panel_tag(ax, "A")

    # B: the reframing - method difference vs donor sampling noise
    ax = axes[1]
    items = [
        ("GCTA vs HEELS\nsame donors", 0.999, AQUA),
        ("GCTA vs itself\ndifferent donors", 0.848, BLUE),
        ("GCTA vs CHEWS\nsame donors", 0.835, ORANGE),
    ]
    ys = np.arange(len(items))[::-1]
    for y, (lab, v, c) in zip(ys, items):
        ax.barh(y, v, color=c, height=0.58, zorder=3)
        ax.text(v - 0.008, y, f"{v:.3f}", ha="right", va="center", fontsize=19,
                fontweight="bold", color="white")
    ax.set_yticks(ys)
    ax.set_yticklabels([lab for lab, _, _ in items], fontsize=17.5)
    ax.set_xlim(0.7, 1.02)
    ax.set_xlabel("Pearson $r$")
    titled(ax, "Method gap vs sampling noise")
    ax.grid(axis="y", visible=False)
    panel_tag(ax, "B", dx=-0.34)

    fig.tight_layout(w_pad=1.6)
    save(fig, "P7_reliability")


# ===========================================================================
# Figures added for the revised (post-review) poster
# ===========================================================================
TBL = HERE / "tables"


def fig_workflow():
    """What CHEWS actually computes, as a strip for the header band."""
    fig, ax = plt.subplots(figsize=(FULLW_IN, 2.4))
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 22)
    ax.axis("off")
    ax.grid(False)

    boxes = [
        (1.0, "marginal\nsummary stats", "$\\hat\\beta$, SE"),
        (25.0, "in-sample or\npanel LD", "$R$"),
        (49.0, "SuSiE posterior\nper SNP", "$\\mu_j,\\ \\sigma_j$"),
        (73.0, "posterior genetic\nvariance", "$\\hat h^2$"),
    ]
    w = 22.0
    for i, (x, label, sym) in enumerate(boxes):
        face = "#eef4fc" if i < 2 else ("#fdf1ea" if i == 2 else "#e8f7f1")
        edge = BLUE if i < 2 else (ORANGE if i == 2 else AQUA)
        ax.add_patch(plt.Rectangle((x, 3.0), w, 15.0, facecolor=face,
                                   edgecolor=edge, linewidth=3.0, zorder=2))
        ax.text(x + w / 2, 13.4, label, ha="center", va="center",
                fontsize=25, fontweight="bold", color=INK, zorder=3)
        ax.text(x + w / 2, 6.6, sym, ha="center", va="center",
                fontsize=26, color=edge, zorder=3)
        if i < len(boxes) - 1:
            ax.annotate("", xy=(x + w + 1.7, 10.5), xytext=(x + w + 0.3, 10.5),
                        arrowprops=dict(arrowstyle="-|>", lw=3.5, color=INK2,
                                        mutation_scale=32))
    fig.tight_layout()
    save(fig, "P0_workflow")


def fig_polygenicity_single():
    """Bias by causal count for all three estimators, plus the quadratic form.

    Merges what were two figures: the win-rate panel duplicated P2, and a
    separate trace-vs-quadratic panel duplicated this x-axis.  One panel carries
    the whole mechanism argument - GCTA and HEELS flat, the trace form falling
    away, and how much of that fall the omitted cross terms account for.
    """
    cells = pd.read_csv(RESULTS / "baseline_n300" / "cells.tsv.gz", sep="\t")
    fig, ax = plt.subplots(figsize=(COLW_IN, 4.7))
    for name in ORDER:
        c = COL[name]
        d = cells.dropna(subset=[f"{c}_h2"]).copy().reset_index(drop=True)
        d["err"] = d[f"{c}_h2"] - d["h2_true"]
        mean = d.groupby("ncausal")["err"].mean()
        keys, lo, hi = region_bootstrap_ci(d, lambda b: b["err"], "ncausal")
        ax.errorbar(keys, [mean[k] for k in keys],
                    yerr=[[mean[k] - lo[k] for k in keys],
                          [hi[k] - mean[k] for k in keys]],
                    color=SERIES[name], marker="o", capsize=6, capthick=2.4,
                    elinewidth=2.4, label=name, zorder=3)

    run = RESULTS / "estimand_check" / "cells.tsv.gz"
    handles = legend_handles()
    if run.exists():
        q = pd.read_csv(run, sep="\t").dropna(subset=["susie_h2_quad"]).reset_index(drop=True)
        q["err"] = q["susie_h2_quad"] - q["h2_true"]
        qm = q.groupby("ncausal")["err"].mean()
        qk, qlo, qhi = region_bootstrap_ci(q, lambda b: b["err"], "ncausal")
        ax.errorbar(qk, [qm[k] for k in qk],
                    yerr=[[qm[k] - qlo[k] for k in qk],
                          [qhi[k] - qm[k] for k in qk]],
                    color="#7d4bd4", marker="s", linestyle="--", capsize=6,
                    capthick=2.4, elinewidth=2.4, zorder=3)
        handles.append(Line2D([0], [0], color="#7d4bd4", lw=4.5, marker="s",
                              ms=11, linestyle="--",
                              label="CHEWS, exact LD-aware"))
    ax.axhline(0, color=MUTED, ls="--", lw=2.2, zorder=1)
    ax.set_xscale("log")
    ax.set_xticks([1, 3, 5, 10, 20])
    ax.set_xticklabels(["1", "3", "5", "10", "20"])
    ax.minorticks_off()
    ax.set_xlabel("causal SNPs in the $\\it{cis}$-window")
    ax.set_ylabel("bias  ($\\hat{h}^2 - h^2$)")
    titled(ax, "CHEWS bias by causal-effect count", size=22)
    ax.legend(handles=handles, loc="lower left", fontsize=17.5)
    fig.tight_layout()
    save(fig, "P1_bias_only")


def fig_agreement():
    """Agreement, not correlation: slope and spread, per sample size."""
    t = pd.read_csv(TBL / "gtex_agreement.csv")
    fig, axes = plt.subplots(1, 2, figsize=(COLW_IN, 4.7))

    for ax, (metric, ylab, ttl, ref) in zip(axes, [
        ("mean_diff", "mean difference\nfrom GCTA", "Systematic offset", 0.0),
        ("slope", "Deming slope vs GCTA", "Calibration slope", 1.0),
    ]):
        for name, col in (("CHEWS", ORANGE), ("HEELS", AQUA)):
            d = t[t["estimator"] == name].sort_values("n")
            ax.plot(d["n"], d[metric], color=col, marker="o", label=name, zorder=3)
        if ref is not None:
            ax.axhline(ref, color=MUTED, ls="--", lw=2.4, zorder=1)
        ax.set_xticks([100, 300, 574])
        ax.set_xlabel("GTEx donors")
        ax.set_ylabel(ylab)
        titled(ax, ttl)
    axes[1].legend(loc="lower right", fontsize=20)
    for ax, L in zip(axes, "AB"):
        panel_tag_in(ax, L)
    fig.tight_layout(w_pad=1.6)
    save(fig, "P5_agreement")


def fig_disjoint_strata():
    """The conditional pattern, with strata defined on donors NOT used for the gap.

    Stacked rather than side by side: two 5x5 grids will not fit across an 11 in
    column without the annotations dropping below the legibility floor.
    """
    mean = pd.read_csv(TBL / "gap_disjoint_strata.csv", index_col=0)
    cnt = pd.read_csv(TBL / "gap_disjoint_counts.csv", index_col=0)
    coup = pd.read_csv(TBL / "gap_coupled_strata.csv", index_col=0)

    fig, axes = plt.subplots(2, 1, figsize=(COLW_IN, 7.4))
    for ax, (M, ttl, tag) in zip(axes, [
        (coup, "Strata from the same donors as the gap", "A"),
        (mean, "Strata from a disjoint half of donors", "B"),
    ]):
        v = M.values.astype(float)
        cmap = plt.get_cmap("RdBu").copy()
        cmap.set_bad("#eceae5")
        im = ax.imshow(np.ma.masked_invalid(v), cmap=cmap, vmin=-0.06, vmax=0.06,
                       aspect="auto")
        ax.set_xticks(range(M.shape[1]))
        ax.set_xticklabels(list(M.columns))
        ax.set_yticks(range(M.shape[0]))
        ax.set_yticklabels([str(i) for i in M.index])
        ax.set_ylabel("independent\n$\\it{cis}$-eQTL signals")
        ax.grid(False)
        for i in range(v.shape[0]):
            for j in range(v.shape[1]):
                if np.isnan(v[i, j]):
                    ax.text(j, i, "n/a", ha="center", va="center",
                            color=MUTED, fontsize=18, style="italic")
                else:
                    dark = abs(v[i, j] / 0.06) > 0.45
                    ax.text(j, i, f"{v[i,j]:+.3f}", ha="center", va="center",
                            fontsize=19, fontweight="bold",
                            color="white" if dark else INK,
                            path_effects=[pe.withStroke(
                                linewidth=2.6,
                                foreground="#00000055" if dark else "#ffffffcc")])
        titled(ax, ttl, size=22)
        panel_tag(ax, tag, dx=-0.155, dy=1.30)
    axes[0].set_xticklabels([])          # shared with B; its title needs the room
    axes[1].set_xlabel("GCTA $\\hat{h}^2$ stratum")
    fig.subplots_adjust(hspace=0.42)
    cb = fig.colorbar(im, ax=axes, fraction=0.038, pad=0.02)
    cb.set_label("mean CHEWS $-$ GCTA", fontsize=21)
    cb.ax.tick_params(labelsize=18)
    save(fig, "P6_strata_compare")


def fig_reliability_sets():
    """Reliability is stable across gene-selection rules; the edge is small."""
    t = pd.read_csv(TBL / "reliability_by_gene_set.csv")
    short = {"all genes (no filter)": "all genes\n(no filter)",
             "prespecified: >=1 cis-eQTL in the full 574": "prespecified\n(≥1 cis-eQTL)",
             "post-hoc: h2>0.05 in all six estimates": "post-hoc\n($\\hat h^2$>0.05 ×6)"}
    sets = list(short)
    fig, ax = plt.subplots(figsize=(COLW_IN, 4.6))
    xs = np.arange(len(sets))
    w = 0.26
    for k, name in enumerate(["GCTA", "CHEWS", "HEELS"]):
        vals = [float(t[(t.gene_set == g) & (t.estimator == name)]["pearson"].iloc[0])
                for g in sets]
        b = ax.bar(xs + (k - 1) * w, vals, width=w, color=SERIES[name],
                   label=name, zorder=3)
        for r, v in zip(b, vals):
            ax.text(r.get_x() + r.get_width() / 2, v + 0.004, f"{v:.3f}",
                    ha="center", va="bottom", fontsize=16, fontweight="bold")
    ax.set_xticks(xs)
    ax.set_xticklabels([short[g] for g in sets])
    ax.set_ylim(0.78, 0.95)
    ax.set_ylabel("split-half Pearson $r$")
    titled(ax, "Reliability ordering is stable across gene sets")
    ax.legend(loc="upper right", ncol=3, fontsize=19)
    fig.tight_layout()
    save(fig, "P7_reliability_sets")


def fig_estimand():
    """Diagonal (trace) vs full quadratic form on the same 3,000 paired fits.

    Same seed as baseline_n300, so the phenotypes are identical and the
    diagonal column reproduces the original run exactly; the only difference is
    the extra mu'R mu term.
    """
    run = RESULTS / "estimand_check" / "cells.tsv.gz"
    if not run.exists():
        print("  (skipping P8: estimand_check has not finished)")
        return
    d = pd.read_csv(run, sep="\t").dropna(subset=["susie_h2", "susie_h2_quad"])
    if not len(d):
        return

    fig, ax = plt.subplots(figsize=(COLW_IN, 4.3))
    series = [("trace  Σ(μ²+σ²)  — as implemented", "susie_h2", ORANGE, "o", "-"),
              ("quadratic  μ'Rμ+Σσ²", "susie_h2_quad", "#7d4bd4", "s", "--")]
    for label, col, c, mk, ls in series:
        g = d.assign(e=d[col] - d["h2_true"]).groupby("ncausal")["e"].agg(["mean", "sem"])
        ax.errorbar(g.index, g["mean"], yerr=1.96 * g["sem"], color=c, marker=mk,
                    linestyle=ls, capsize=6, capthick=2.4, elinewidth=2.4,
                    label=label, zorder=3)
    ax.axhline(0, color=MUTED, ls="--", lw=2.2, zorder=1)
    ax.set_xscale("log")
    ax.set_xticks([1, 3, 5, 10, 20])
    ax.set_xticklabels(["1", "3", "5", "10", "20"])
    ax.minorticks_off()
    ax.set_xlabel("causal SNPs in the $\\it{cis}$-window")
    ax.set_ylabel("bias  ($\\hat{h}^2 - h^2$)")
    titled(ax, "The omitted term is not small", size=22)
    ax.legend(loc="lower left", fontsize=18)
    fig.tight_layout()
    save(fig, "P8_estimand")


# ===========================================================================
# Round 4 (PI feedback) - docs/pi_feedback_round4_design.md.  New files only;
# every round-3 figure above is still built, unchanged.
# ===========================================================================
def _win_cell(ax, j, i, rate_pct, sub=None, sub_size=16.5, main_size=19):
    """Win-rate label, with an optional second line (relative MAE).  Same
    dark/white switching rule as cell_text, on the shared 15-85% ramp."""
    norm_v = (rate_pct - 50) / ((WIN_VMAX - WIN_VMIN) / 2)
    dark = abs(norm_v) > 0.45
    kw = dict(ha="center", fontweight="bold", color="white" if dark else INK,
              path_effects=[pe.withStroke(linewidth=2.6,
                                          foreground="#00000055" if dark else "#ffffffcc")])
    if sub is None:
        ax.text(j, i, f"{rate_pct:.0f}%", va="center", fontsize=main_size, **kw)
        return
    ax.text(j, i - 0.04, f"{rate_pct:.0f}%", va="bottom", fontsize=main_size, **kw)
    kw["fontweight"] = "normal"
    ax.text(j, i + 0.06, sub, va="top", fontsize=sub_size, **kw)


def fig_sample_size_winrate(two_line=False, name="P3_sample_size_r4"):
    """§5 with panel B redesigned: how OFTEN CHEWS wins, not by how much.

    Panel A is drawn by the same helper as round 3.  Panel B uses the P2
    encoding (WIN_CMAP, WIN_VMIN-WIN_VMAX, centred at 50%).  With
    ``two_line`` the relative MAE change rides under each win rate (design
    option A3); without it, only the win rate (A2).
    """
    m = pd.read_csv(RESULTS / "sample_size" / "tables" / "metrics_by_n.csv")
    m["estimator"] = m["estimator"].replace({"SuSiE": "CHEWS"})
    t = pd.read_csv(TBL / "win_by_n_causal.csv")

    fig, axes = plt.subplots(1, 2, figsize=(COLW_IN, 4.9),
                             gridspec_kw=dict(width_ratios=[1.0, 1.30]))
    _panel_dispersion(axes[0], m)
    # Lay out the ROUND-3 figure, freeze its geometry, then swap panel B in
    # place - so panel A sits exactly where it did (PI: "keep as is").
    # tight_layout is not called again.
    cb_old = _panel_mae_diff(fig, axes[1], m)
    fig.tight_layout(w_pad=1.6)
    pos_b, pos_cb = axes[1].get_position(), cb_old.ax.get_position()
    cb_old.remove()
    axes[1].remove()
    ax = fig.add_axes(pos_b)
    cax = fig.add_axes(pos_cb)
    cax.grid(False)          # poster rcParams turn grids on for new axes
    rate = (100 * t.pivot(index="ncausal", columns="n", values="rate")).sort_index(ascending=False)
    rel = t.pivot(index="ncausal", columns="n", values="rel_mae").reindex(rate.index)
    lo = t.pivot(index="ncausal", columns="n", values="lo").reindex(rate.index)
    hi = t.pivot(index="ncausal", columns="n", values="hi").reindex(rate.index)
    v = rate.values.astype(float)
    im = ax.imshow(v, cmap=WIN_CMAP, vmin=WIN_VMIN, vmax=WIN_VMAX, aspect="auto")
    ax.set_xticks(range(rate.shape[1]))
    ax.set_xticklabels([str(int(c)) for c in rate.columns])
    ax.set_yticks(range(rate.shape[0]))
    ax.set_yticklabels([str(int(i)) for i in rate.index], fontsize=22, fontweight="bold")
    ax.set_xlabel("individuals ($n$)")
    ax.set_ylabel("causal SNPs")
    ax.grid(False)
    for i in range(v.shape[0]):
        for j in range(v.shape[1]):
            sub = None
            if two_line:
                sub = f"{100 * rel.values[i, j]:+.0f}%".replace("-", "\u2212")
            _win_cell(ax, j, i, v[i, j], sub)
            # cells whose region-bootstrap CI includes 50% are drawn with a
            # dotted frame: "not distinguishable from a coin flip"
            if lo.values[i, j] <= 0.5 <= hi.values[i, j]:
                ax.add_patch(plt.Rectangle((j - 0.47, i - 0.47), 0.94, 0.94, fill=False,
                                           ec=INK2, lw=1.8, ls=(0, (2, 2)), zorder=4))
    cb = fig.colorbar(im, cax=cax)
    cb.set_label("% of loci where CHEWS\nis closer to truth", fontsize=17)
    cb.ax.tick_params(labelsize=16)
    # one line, so it stays within the height round 3's title used
    titled(ax, "How often CHEWS has lower error", size=21)
    panel_tag(ax, "B", dx=-0.33, dy=1.34)
    save(fig, name)


def fig_better_or_different():
    """§7, round 4: is CHEWS *better* than GCTA in each stratum, not just different?

    A = the observed disjoint-donor gap (round-3 panel B, same scale).
    B = the semi-synthetic GTEx run (same donors, halves, covariates, strata):
        % of genes where CHEWS is closer to truth under the best-fit prior,
        on the shared win-rate encoding.  Hatched = the acceptable priors do not
        agree on the side of 50% (part_g decision rule 4); grey = ESS < 30.
    """
    mean = pd.read_csv(TBL / "gap_disjoint_strata.csv", index_col=0)
    win = pd.read_csv(TBL / "sim_strata_win.csv")
    win["clump_bin"] = win["clump_bin"].astype(str)
    order = [str(i) for i in mean.index]
    W = (100 * win.pivot(index="clump_bin", columns="h2_bin", values="win")).reindex(order)[list(mean.columns)]
    masked = win.pivot(index="clump_bin", columns="h2_bin", values="masked").reindex(order)[list(mean.columns)]
    rule = win.pivot(index="clump_bin", columns="h2_bin", values="rule").reindex(order)[list(mean.columns)]

    fig, axes = plt.subplots(2, 1, figsize=(COLW_IN, 7.4))
    cmap = plt.get_cmap("RdBu").copy()
    cmap.set_bad("#eceae5")
    ax = axes[0]
    v = mean.values.astype(float)
    imA = ax.imshow(np.ma.masked_invalid(v), cmap=cmap, vmin=-0.06, vmax=0.06, aspect="auto")
    for i in range(v.shape[0]):
        for j in range(v.shape[1]):
            if np.isnan(v[i, j]):
                ax.text(j, i, "n/a", ha="center", va="center", color=MUTED, fontsize=18,
                        style="italic")
            else:
                dark = abs(v[i, j] / 0.06) > 0.45
                ax.text(j, i, f"{v[i,j]:+.3f}", ha="center", va="center", fontsize=19,
                        fontweight="bold", color="white" if dark else INK,
                        path_effects=[pe.withStroke(linewidth=2.6,
                                                    foreground="#00000055" if dark else "#ffffffcc")])
    titled(ax, "Observed: CHEWS $-$ GCTA (strata from a disjoint half)", size=21)
    cb = fig.colorbar(imA, ax=ax, fraction=0.046, pad=0.02)
    cb.set_label("mean CHEWS $-$ GCTA", fontsize=18)
    cb.ax.tick_params(labelsize=16)

    ax = axes[1]
    w = W.values.astype(float)
    wm = np.where(masked.values.astype(bool), np.nan, w)
    cmapW = plt.get_cmap(WIN_CMAP).copy()
    cmapW.set_bad("#eceae5")
    imB = ax.imshow(np.ma.masked_invalid(wm), cmap=cmapW, vmin=WIN_VMIN, vmax=WIN_VMAX,
                    aspect="auto")
    for i in range(w.shape[0]):
        for j in range(w.shape[1]):
            if np.isnan(wm[i, j]):
                ax.text(j, i, "n/a", ha="center", va="center", color=MUTED, fontsize=18,
                        style="italic")
                continue
            if str(rule.values[i, j]).startswith("4"):
                ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, ec=INK2,
                                           lw=0, hatch="//", alpha=0.55, zorder=3))
            _win_cell(ax, j, i, wm[i, j])
    titled(ax, "Simulated on GTEx genotypes: % genes\nwhere CHEWS is closer to truth", size=21)
    import json
    if not json.load(open(TBL / "calibration_summary.json"))["calibration_passes"]:
        # pre-registered rule 1 failed: this panel must not be read as a result
        fig.text(0.5, -0.02, "EXPLORATORY: simulation calibration failed (rule 1) - not for the poster",
                 ha="center", va="top", fontsize=17, color="#b3261e", fontweight="bold")
    cb = fig.colorbar(imB, ax=ax, fraction=0.046, pad=0.02)
    cb.set_label("% CHEWS closer", fontsize=18)
    cb.ax.tick_params(labelsize=16)

    for ax, tag in zip(axes, "AB"):
        ax.set_xticks(range(mean.shape[1]))
        ax.set_xticklabels(list(mean.columns))
        ax.set_yticks(range(mean.shape[0]))
        ax.set_yticklabels(order)
        ax.set_ylabel("independent\n$\\it{cis}$-eQTL signals")
        ax.grid(False)
        panel_tag(ax, tag, dx=-0.155, dy=1.30 if tag == "A" else 1.42)
    axes[0].set_xticklabels([])
    axes[1].set_xlabel("half-A GCTA $\\hat{h}^2$ stratum")
    fig.subplots_adjust(hspace=0.62)
    save(fig, "P6_better_or_different")


if __name__ == "__main__":
    import sys

    round4 = "--round" in sys.argv and sys.argv[sys.argv.index("--round") + 1] == "4"
    print("building poster figures ->", OUT)
    fig_polygenicity()
    fig_landscape()
    fig_sample_size()
    fig_sufficiency()
    fig_gtex()
    fig_gtex_architecture()
    fig_reliability()
    # --- revised-poster figures ---
    fig_workflow()
    fig_polygenicity_single()
    fig_agreement()
    fig_disjoint_strata()
    fig_reliability_sets()
    if round4:
        # round-4 figures are additive: new file names, nothing above changes
        # design §2.3: A3 (win rate over relative MAE) was prototyped and failed
        # the legibility rule - labels collide across cells at >=16 pt - so A2
        # ships.  The A3 prototype is kept for the record.
        fig_sample_size_winrate(two_line=False, name="P3_sample_size_r4")
        fig_sample_size_winrate(two_line=True, name="P3_sample_size_r4_A3_rejected")
        if (TBL / "sim_strata_win.csv").exists():
            fig_better_or_different()
        else:
            print("  (skipping P6_better_or_different: run revision_analyses.py f g first)")
    print("done")
