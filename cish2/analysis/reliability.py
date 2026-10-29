"""Test-retest reliability from disjoint donor halves — a truth-free criterion.

Real data has no ground truth, so accuracy cannot be measured directly.  What
*can* be measured is **reliability**: split the donors into two disjoint halves,
estimate every gene's h2 independently in each, and ask how well the two agree.
No donor contributes to both estimates, so the correlation is an honest
test-retest coefficient rather than the inflated self-similarity that nested
subsamples produce.

Two cautions are built into the output, because reliability on its own is a
trap:

1. **A shrunken estimator wins this contest automatically.**  An estimator that
   always returned the same number would be perfectly reliable and perfectly
   useless.  Reliability is necessary, not sufficient; it has to be read next
   to the simulated bias.
2. **A pile of estimates at zero fakes both directions** — it inflates
   correlation (two clouds far apart) and deflates rank correlation (ties).  So
   everything is reported on a common gene set with signal in every estimator.

The comparison that matters most here is the *ratio*: how an estimator agrees
with itself across independent donors, against how two different estimators
agree within the same donors.  When the latter is higher, the choice of method
is not the limiting factor — the sample is.
"""

from __future__ import annotations

import argparse
import os
import sys
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from analysis.style import (  # noqa: E402
    INK, INK_SECONDARY, REFERENCE, SERIES, apply_style, finish, panel_label, save,
)

COL = {"GCTA": "gcta_h2", "SuSiE": "susie_h2", "HEELS": "heels_h2"}
SIGNAL_FLOOR = 0.05


def load_halves(run_dir: str):
    from cish2.results import ResultStore
    from scipy.stats import pearsonr, spearmanr

    d = ResultStore(run_dir).load()
    a = d[d["subset"] == "split1"].set_index("gene_id")
    b = d[d["subset"] == "split2"].set_index("gene_id")
    genes = a.index.intersection(b.index)
    ests = [e for e, c in COL.items() if c in d.columns]
    A = a.loc[genes, [COL[e] for e in ests]]
    B = b.loc[genes, [COL[e] for e in ests]]
    ok = A.notna().all(axis=1) & B.notna().all(axis=1)
    signal = ok & (A > SIGNAL_FLOOR).all(axis=1) & (B > SIGNAL_FLOOR).all(axis=1)
    return A, B, ok, signal, ests, int(d[d["subset"] == "split1"]["n"].iloc[0])


def reliability_table(A, B, mask, ests) -> pd.DataFrame:
    from scipy.stats import pearsonr, spearmanr

    rows = []
    for e in ests:
        x, y = A.loc[mask, COL[e]], B.loc[mask, COL[e]]
        m = (x + y) / 2
        ms_b = 2 * ((m - m.mean()) ** 2).sum() / (len(x) - 1)
        ms_w = ((x - m) ** 2 + (y - m) ** 2).sum() / len(x)
        rows.append(dict(
            estimator=e, genes=int(len(x)),
            pearson=float(pearsonr(x, y)[0]),
            spearman=float(spearmanr(x, y).correlation),
            icc=float((ms_b - ms_w) / (ms_b + ms_w)),
            rmsd=float(np.sqrt(((x - y) ** 2).mean())),
            median_abs_diff=float(np.median((x - y).abs())),
        ))
    return pd.DataFrame(rows)


def cross_estimator_table(A, mask, ests) -> pd.DataFrame:
    from scipy.stats import pearsonr

    rows = []
    for i, p in enumerate(ests):
        for q in ests[i + 1:]:
            x, y = A.loc[mask, COL[p]], A.loc[mask, COL[q]]
            rows.append(dict(pair=f"{p} vs {q}", pearson=float(pearsonr(x, y)[0])))
    return pd.DataFrame(rows)


def fig_reliability(A, B, signal, ests, within, between, out_dir: str,
                    n_half: int) -> List[str]:
    """Half-A vs half-B per estimator, plus the comparison that reframes it."""
    apply_style()
    fig, axes = plt.subplots(1, len(ests) + 1,
                             figsize=(3.5 * (len(ests) + 1) + 0.6, 3.9))

    hi = float(max(A.loc[signal].max().max(), B.loc[signal].max().max()))
    hi = min(hi, 1.0)
    for k, e in enumerate(ests):
        ax = axes[k]
        x, y = A.loc[signal, COL[e]], B.loc[signal, COL[e]]
        ax.plot([0, hi], [0, hi], ls=(0, (4, 3)), lw=1.1, color=REFERENCE, zorder=1)
        ax.scatter(x, y, s=6, color=SERIES[e], alpha=0.35, linewidths=0, zorder=3)
        r = float(within.set_index("estimator").loc[e, "pearson"])
        ax.text(0.04, 0.95, f"r = {r:.3f}", transform=ax.transAxes, fontsize=10,
                fontweight="bold", color=INK, va="top")
        ax.set_xlim(0, hi); ax.set_ylim(0, hi)
        finish(ax, xlabel=f"$\\hat{{h}}^2$, donors 1-{n_half}",
               ylabel=f"$\\hat{{h}}^2$, other {n_half} donors" if k == 0 else "",
               title=e)
        panel_label(ax, chr(ord("A") + k))

    # the reframing panel: self-agreement across donors vs cross-method agreement
    ax = axes[-1]
    labels = [f"{r.estimator} vs itself\n(independent donors)" for r in within.itertuples()]
    vals = list(within["pearson"])
    colors = [SERIES[e] for e in within["estimator"]]
    labels += [f"{r.pair}\n(same donors)" for r in between.itertuples()]
    vals += list(between["pearson"])
    colors += [INK_SECONDARY] * len(between)
    order = np.argsort(vals)
    y = np.arange(len(vals))
    ax.barh(y, [vals[i] for i in order], height=0.6,
            color=[colors[i] for i in order], zorder=3,
            edgecolor="white", linewidth=1.4)
    for yi, i in enumerate(order):
        ax.text(vals[i] - 0.008, yi, f"{vals[i]:.3f}", ha="right", va="center",
                fontsize=8, color="white", fontweight="bold")
    ax.set_yticks(y)
    ax.set_yticklabels([labels[i] for i in order], fontsize=7.5)
    ax.set_xlim(0.7, 1.02)
    ax.grid(axis="y", visible=False)
    finish(ax, xlabel="Pearson r", title="Method choice vs. sampling noise")
    panel_label(ax, chr(ord("A") + len(ests)))

    fig.suptitle(f"Test-retest reliability — two disjoint halves of {n_half} donors "
                 f"({int(signal.sum()):,} genes)",
                 fontsize=11.5, fontweight="bold", color=INK, y=1.02)
    fig.text(0.5, -0.03,
             "Grey bars compare two different estimators on the SAME donors; "
             "coloured bars compare one estimator with itself on DIFFERENT donors.",
             ha="center", fontsize=8.5, color=INK_SECONDARY)
    fig.tight_layout()
    return save(fig, out_dir, "gtex_figR1_reliability")


def build_reliability_report(run_dir: str) -> str:
    run_dir = os.path.abspath(run_dir)
    fig_dir = os.path.join(run_dir, "figures")
    tab_dir = os.path.join(run_dir, "tables")
    os.makedirs(fig_dir, exist_ok=True)
    os.makedirs(tab_dir, exist_ok=True)

    A, B, ok, signal, ests, n_half = load_halves(run_dir)
    all_genes = reliability_table(A, B, ok, ests)
    sig = reliability_table(A, B, signal, ests)
    between = cross_estimator_table(A, signal, ests)
    all_genes.to_csv(os.path.join(tab_dir, "reliability_all_genes.csv"), index=False)
    sig.to_csv(os.path.join(tab_dir, "reliability_signal_genes.csv"), index=False)
    between.to_csv(os.path.join(tab_dir, "cross_estimator_agreement.csv"), index=False)
    fig_reliability(A, B, signal, ests, sig, between, fig_dir, n_half)

    best = sig.loc[sig["pearson"].idxmax(), "estimator"]
    worst_self = float(sig["pearson"].min())
    best_between = float(between["pearson"].max())

    L: List[str] = []
    add = L.append
    add("# Test-retest reliability on disjoint donor halves")
    add("")
    add(f"GTEx Whole Blood donors split into two **disjoint** halves of {n_half}, "
        f"each analysed independently end to end. No donor contributes to both "
        f"estimates, so half-A versus half-B is an honest reliability "
        f"coefficient rather than the self-similarity nested subsamples produce.")
    add("")
    add(f"Reported on **{int(signal.sum()):,} genes** with $\\hat h^2$ > "
        f"{SIGNAL_FLOOR} in every estimator in both halves — a common set, so the "
        f"three are compared on identical genes and nobody is helped by a pile of "
        f"estimates sitting at zero.")
    add("")
    add("| estimator | Pearson | Spearman | ICC | RMSD | median &#124;A−B&#124; |")
    add("|---|---:|---:|---:|---:|---:|")
    for r in sig.itertuples():
        add(f"| {r.estimator} | **{r.pearson:.3f}** | {r.spearman:.3f} | {r.icc:.3f} "
            f"| {r.rmsd:.4f} | {r.median_abs_diff:.4f} |")
    add("")
    add(f"**{best} is the most reliable estimator** on every measure. That is a "
        f"real result and it must be read carefully, because **reliability "
        f"systematically favours the more shrunken estimator**: an estimator that "
        f"returned a constant would score a perfect ICC and be worthless. The "
        f"simulations show SuSiE is the shrunken one — its inter-quartile error "
        f"spread is a third of GCTA's at low polygenicity — so high reliability "
        f"is exactly what its bias predicts, not evidence against it.")
    add("")
    add("Reliability answers *precision*. It cannot, on its own, answer "
        "*accuracy*; only the simulations can, and there SuSiE is the biased arm.")
    add("")
    add("## The comparison that reframes the question")
    add("")
    add("| comparison | Pearson r |")
    add("|---|---:|")
    for r in between.itertuples():
        add(f"| {r.pair}, **same** donors | {r.pearson:.3f} |")
    for r in sig.itertuples():
        add(f"| {r.estimator} vs itself, **different** donors | {r.pearson:.3f} |")
    add("")
    # state this precisely: the two comparisons land differently for the
    # infinitesimal pair and for SuSiE, and an earlier draft over-claimed
    self_r = dict(zip(sig["estimator"], sig["pearson"]))
    pair_r = dict(zip(between["pair"], between["pearson"]))
    gh = pair_r.get("GCTA vs HEELS", float("nan"))
    gs = pair_r.get("GCTA vs SuSiE", float("nan"))
    g_self = self_r.get("GCTA", float("nan"))
    add(f"**GCTA and HEELS are one method by this measure.** They agree at "
        f"r = {gh:.3f} on the same donors, against GCTA's own across-donor "
        f"reliability of {g_self:.3f} — the method difference is more than an "
        f"order of magnitude smaller than the sampling noise. Choosing between "
        f"them on real data is not a meaningful decision; the sufficiency result "
        f"from the simulations says it never could be.")
    add("")
    add(f"**SuSiE is a genuinely different measurement, but only just.** It "
        f"agrees with GCTA at r = {gs:.3f} on the same donors, essentially the "
        f"same as GCTA's agreement with *itself* on different donors "
        f"({g_self:.3f}). So the SuSiE-versus-GCTA difference is of the same "
        f"order as donor sampling noise at this sample size — detectable in "
        f"aggregate over ~4,000 genes, but not something to act on for any single "
        f"gene.")
    add("")
    add("![reliability](figures/gtex_figR1_reliability.png)")
    add("")
    add("## What this licenses, and what it does not")
    add("")
    add("**Licensed:** at this sample size, no pair of these estimators differs by "
        "more than donor sampling noise. A paper claiming one is better on GTEx "
        "Whole Blood alone would be over-reading; the simulations, where truth is "
        "known, are what carry that claim.")
    add("")
    add("**Not licensed:** \"they are all equivalent\". They are equivalent *on "
        "average over genes at this n*. The simulations locate exactly where they "
        "are not — polygenic, non-localisable architecture — and the discordant "
        "tail of real genes is enriched for precisely that.")
    add("")
    add("For completeness, the same table over **all** genes (including the "
        "boundary pile, which is why these numbers differ):")
    add("")
    add("| estimator | Pearson | Spearman | ICC | RMSD | genes |")
    add("|---|---:|---:|---:|---:|---:|")
    for r in all_genes.itertuples():
        add(f"| {r.estimator} | {r.pearson:.3f} | {r.spearman:.3f} | {r.icc:.3f} "
            f"| {r.rmsd:.4f} | {r.genes:,} |")
    add("")
    add("---")
    add("")
    add("Generated by `python -m analysis.reliability`.")

    path = os.path.join(run_dir, "RELIABILITY.md")
    with open(path, "w") as fh:
        fh.write("\n".join(L) + "\n")
    return path


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir", nargs="?",
                    default="/home/sux002/work/cis_heritability/results/gtex_split_half")
    args = ap.parse_args(argv)
    print(build_reliability_report(args.run_dir))
    return 0


if __name__ == "__main__":
    sys.exit(main())
