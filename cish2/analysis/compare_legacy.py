"""Check a new campaign against the legacy result tables.

The previous pipeline wrote one flat TSV per script run.  This module parses
those tables, puts them on the same footing as a `cish2` campaign, and reports
whether the rebuilt engine lands on the same numbers.

The comparison is between *cell means over independent replicate loci*, not
row-by-row: the two pipelines draw different random regions, so agreement is
expected within Monte-Carlo error, not exactly.

    python -m analysis.compare_legacy \\
        --new  ../results/baseline_n300 \\
        --legacy ../simulations/result_eur_num_people_100_300_regular_LD.txt

Note the legacy ``hsq`` column is GCTA's **V(G)**, not V(G)/Vp, so it is
compared against ``gcta_vg``.
"""

from __future__ import annotations

import argparse
import os
import sys
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from analysis.load import load_run  # noqa: E402
from analysis.style import (  # noqa: E402
    INK,
    INK_SECONDARY,
    REFERENCE,
    SERIES,
    apply_style,
    finish,
    panel_label,
    save,
)


# The legacy tables predate HEELS, so this comparison is fixed to the two
# estimators they recorded; a third arm has nothing to be compared against.
LEGACY_ESTIMATORS = ("GCTA", "SuSiE")


def read_legacy(path: str) -> pd.DataFrame:
    """Parse a legacy result table into (ncausal, h2_true, gcta_vg, susie_h2)."""
    df = pd.read_csv(path, sep="\t")
    if "parameter" in df.columns:  # "<ncausal>_<h2>"
        parts = df["parameter"].astype(str).str.split("_", n=1, expand=True)
        df["ncausal"] = parts[0].astype(int)
        df["h2_true"] = parts[1].astype(float)
    for col in ("ncausal", "h2"):
        if col in df.columns and col not in ("h2_true",):
            df[col] = pd.to_numeric(df[col], errors="coerce")
    if "h2" in df.columns and "h2_true" not in df.columns:
        df["h2_true"] = df["h2"].astype(float)

    out = pd.DataFrame({
        "ncausal": df["ncausal"].astype(int),
        "h2_true": df["h2_true"].astype(float),
        "gcta_vg": pd.to_numeric(df["hsq"], errors="coerce"),
        "gcta_vg_se": pd.to_numeric(df["hsq_se"], errors="coerce"),
        "susie_h2": pd.to_numeric(df["su_hsq"], errors="coerce"),
    })
    for col in ("effect_dist", "maf_bin"):
        if col in df.columns:
            out[col] = df[col].astype(str)
    # the legacy sentinel for a failed REML fit
    out.loc[out["gcta_vg_se"] >= 10000, ["gcta_vg", "gcta_vg_se"]] = np.nan
    return out


def cell_means(df: pd.DataFrame, gcta_col: str, susie_col: str,
               extra: Optional[List[str]] = None) -> pd.DataFrame:
    keys = ["ncausal", "h2_true"] + [c for c in (extra or []) if c in df.columns]
    grouped = df.groupby(keys, observed=True)
    out = grouped.agg(
        n_cells=(gcta_col, "size"),
        GCTA=(gcta_col, "mean"),
        GCTA_sd=(gcta_col, "std"),
        SuSiE=(susie_col, "mean"),
        SuSiE_sd=(susie_col, "std"),
    ).reset_index()
    for est in LEGACY_ESTIMATORS:
        out[f"{est}_mcse"] = out[f"{est}_sd"] / np.sqrt(out["n_cells"])
    return out


def compare(new_run: str, legacy_path: str, out_dir: Optional[str] = None,
            extra: Optional[List[str]] = None) -> Dict[str, object]:
    new = load_run(new_run, drop_failed=False)
    old = read_legacy(legacy_path)

    a = cell_means(new, "gcta_vg", "susie_h2", extra).add_suffix("_new")
    b = cell_means(old, "gcta_vg", "susie_h2", extra).add_suffix("_old")
    keys = ["ncausal", "h2_true"] + [c for c in (extra or []) if c in old.columns]
    a = a.rename(columns={f"{k}_new": k for k in keys})
    b = b.rename(columns={f"{k}_old": k for k in keys})
    merged = a.merge(b, on=keys, how="inner")

    for est in LEGACY_ESTIMATORS:
        merged[f"{est}_diff"] = merged[f"{est}_new"] - merged[f"{est}_old"]
        # combined Monte-Carlo SE of the difference of two independent means
        merged[f"{est}_diff_mcse"] = np.sqrt(
            merged[f"{est}_mcse_new"] ** 2 + merged[f"{est}_mcse_old"] ** 2
        )
        merged[f"{est}_z"] = merged[f"{est}_diff"] / merged[f"{est}_diff_mcse"]

    summary = {
        "cells_compared": int(len(merged)),
        "legacy_rows": int(len(old)),
        "new_rows": int(len(new)),
    }
    for est in LEGACY_ESTIMATORS:
        summary[f"{est}_mean_abs_diff"] = float(merged[f"{est}_diff"].abs().mean())
        summary[f"{est}_max_abs_diff"] = float(merged[f"{est}_diff"].abs().max())
        summary[f"{est}_corr"] = float(
            np.corrcoef(merged[f"{est}_new"], merged[f"{est}_old"])[0, 1]
        )
        summary[f"{est}_frac_within_2mcse"] = float((merged[f"{est}_z"].abs() <= 2).mean())

    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
        merged.to_csv(os.path.join(out_dir, "legacy_comparison.csv"), index=False)
        _plot(merged, out_dir, os.path.basename(legacy_path))
    return {"summary": summary, "table": merged}


def _plot(merged: pd.DataFrame, out_dir: str, legacy_name: str) -> List[str]:
    apply_style()
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.2))

    ax = axes[0]
    vals = pd.concat([merged[f"{e}_{s}"] for e in LEGACY_ESTIMATORS for s in ("new", "old")])
    lims = [float(vals.min()) - 0.03, float(vals.max()) + 0.03]
    ax.plot(lims, lims, ls=(0, (4, 3)), lw=1.1, color=REFERENCE, zorder=1)
    for est in LEGACY_ESTIMATORS:
        ax.errorbar(merged[f"{est}_old"], merged[f"{est}_new"],
                    xerr=1.96 * merged[f"{est}_mcse_old"],
                    yerr=1.96 * merged[f"{est}_mcse_new"],
                    fmt="o", ms=5, color=SERIES[est], elinewidth=0.9, capsize=0,
                    markeredgecolor="white", markeredgewidth=0.6, alpha=0.9, zorder=3)
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    finish(ax, xlabel="legacy pipeline, cell mean $\\hat{h}^2$",
           ylabel="cish2, cell mean $\\hat{h}^2$",
           title="Cell means agree across pipelines")
    panel_label(ax, "A")
    from matplotlib.lines import Line2D
    ax.legend(handles=[Line2D([], [], color=REFERENCE, ls=(0, (4, 3)), lw=1.1, label="identity")]
              + [Line2D([], [], color=SERIES[e], marker="o", ls="none", ms=5, label=e)
                 for e in LEGACY_ESTIMATORS], loc="upper left")

    ax = axes[1]
    ax.axhline(0, color=REFERENCE, lw=1.1, ls=(0, (4, 3)), zorder=1)
    for band, alpha in ((2, 0.10),):
        ax.axhspan(-band, band, color=INK_SECONDARY, alpha=alpha, lw=0, zorder=0)
    x = np.arange(len(merged))
    order = merged.sort_values(["ncausal", "h2_true"]).index
    for est in LEGACY_ESTIMATORS:
        ax.scatter(x, merged.loc[order, f"{est}_z"], s=26, color=SERIES[est],
                   edgecolor="white", linewidth=0.6, zorder=3)
    ax.set_xticks([])
    finish(ax, xlabel="design cell (ordered by causal SNPs, then $h^2$)",
           ylabel="difference / Monte-Carlo SE",
           title="Differences vs. sampling noise")
    ax.text(0.99, 0.03, "shaded: within 2 MC SE", transform=ax.transAxes, ha="right",
            va="bottom", fontsize=8, color=INK_SECONDARY)
    panel_label(ax, "B")

    fig.suptitle(f"Reproduction check against `{legacy_name}`", fontsize=11.5,
                 fontweight="bold", color=INK, y=1.01)
    fig.tight_layout()
    return save(fig, out_dir, "fig0_legacy_reproduction")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--new", required=True, help="cish2 run directory")
    ap.add_argument("--legacy", required=True, help="legacy result .txt")
    ap.add_argument("--out", default=None, help="output dir (default <new>/figures)")
    ap.add_argument("--extra", action="append", default=[],
                    help="extra grouping column present in both tables")
    args = ap.parse_args(argv)

    out_dir = args.out or os.path.join(args.new, "figures")
    result = compare(args.new, args.legacy, out_dir, args.extra)
    summary = result["summary"]
    width = max(len(k) for k in summary)
    for key, value in summary.items():
        shown = f"{value:.4f}" if isinstance(value, float) else f"{value:,}"
        print(f"{key:<{width}}  {shown}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
