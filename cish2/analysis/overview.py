"""Cross-campaign overview, reported one *condition* at a time.

    python -m analysis.overview ../results

A condition is one level of one swept axis (see :mod:`analysis.conditions`),
so nothing is averaged over a sweep: `n = 100 (below baseline)` and
`n = 1000 (above baseline)` are separate rows and separate panels, never a
single "sample_size" entry whose direction the reader has to guess.

Writes ``OVERVIEW.md``, ``overview.csv``, ``conditions.csv`` and the
comparison figures under ``<results_root>/figures/``.
"""

from __future__ import annotations

import argparse
import glob
import os
import sys
from typing import Dict, List

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from analysis.conditions import baseline_condition, collect  # noqa: E402
from analysis.figures import fig_condition_effects, fig_condition_grid  # noqa: E402
from analysis.load import load_run, varying_axes  # noqa: E402
from analysis.metrics import summarize_by, win_rate  # noqa: E402

SKIP = {"smoke", ".cache", "scratch", "figures"}


def is_simulation_run(path: str) -> bool:
    """True for a simulation campaign, False for a real-data run.

    Real-data runs (GTEx) share the shard layout but have no true heritability,
    so every simulation metric is undefined for them; they get their own report
    via :mod:`analysis.gtex_report`.
    """
    import gzip

    shards = sorted(glob.glob(os.path.join(path, "shards", "*.tsv.gz")))
    if not shards:
        return False
    with gzip.open(shards[0], "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
    return "h2_true" in header


def discover(results_root: str, simulation_only: bool = True) -> List[str]:
    """Campaign directories that hold at least one result shard."""
    runs = []
    for name in sorted(os.listdir(results_root)):
        path = os.path.join(results_root, name)
        if not os.path.isdir(path) or name in SKIP or name.startswith("."):
            continue
        shards = os.path.join(path, "shards")
        if not (os.path.isdir(shards) and os.listdir(shards)):
            continue
        if simulation_only and not is_simulation_run(path):
            continue
        runs.append(path)
    return runs


def campaign_table(runs: List[str]) -> pd.DataFrame:
    """One row per campaign directory — the provenance view."""
    rows = []
    for run in runs:
        wide = load_run(run)
        rows.append(dict(
            campaign=os.path.basename(run),
            regions=int(wide["region_id"].nunique()),
            cells=int(len(wide)),
            median_p=int(wide["p_snps"].median()),
            axes=", ".join(varying_axes(wide)),
            ld_for_susie=str(wide["ld_for_susie"].iloc[0]),
        ))
    return pd.DataFrame(rows)


def condition_table(conditions) -> pd.DataFrame:
    """One row per swept level — the scientific view."""
    rows = []
    for cond in conditions:
        overall = summarize_by(cond.long, []).set_index("estimator")
        sparse = summarize_by(cond.long[cond.long["ncausal"] == cond.long["ncausal"].min()], [])
        poly = summarize_by(cond.long[cond.long["ncausal"] == cond.long["ncausal"].max()], [])
        sparse = sparse.set_index("estimator")
        poly = poly.set_index("estimator")
        wins = win_rate(cond.wide, [])
        rec = dict(
            group=cond.group,
            condition=cond.label,
            direction=cond.direction,
            campaign=cond.campaign,
            cells=int(len(cond.wide)),
            regions=int(cond.wide["region_id"].nunique()),
        )
        for est in ("GCTA", "SuSiE"):
            if est in overall.index:
                rec[f"{est}_bias"] = float(overall.loc[est, "bias"])
                rec[f"{est}_bias_mcse"] = float(overall.loc[est, "bias_mcse"])
                rec[f"{est}_rmse"] = float(overall.loc[est, "rmse"])
                rec[f"{est}_mae"] = float(overall.loc[est, "mae"])
                rec[f"{est}_frac_impossible"] = float(overall.loc[est, "frac_impossible"])
                rec[f"{est}_max_abs_err"] = float(overall.loc[est, "max_abs_err"])
                rec[f"{est}_cov95"] = float(overall.loc[est, "coverage95"])
                rec[f"{est}_bias_sparse"] = float(sparse.loc[est, "bias"])
                rec[f"{est}_bias_polygenic"] = float(poly.loc[est, "bias"])
        rec["susie_win_rate"] = float(wins["susie_win_rate"].iloc[0])
        rec["susie_win_mcse"] = float(wins["susie_win_mcse"].iloc[0])
        rows.append(rec)
    return pd.DataFrame(rows)


def build_overview(results_root: str) -> str:
    results_root = os.path.abspath(results_root)
    runs = discover(results_root)
    if not runs:
        raise FileNotFoundError(f"no campaigns with results under {results_root}")

    conditions = collect(runs)
    baseline = baseline_condition(conditions)
    if baseline is None:
        raise RuntimeError("no baseline condition found — expected an unselected-LD, "
                           "n=300, normal-effects campaign")

    camps = campaign_table(runs)
    conds = condition_table(conditions)
    camps.to_csv(os.path.join(results_root, "overview.csv"), index=False)
    conds.to_csv(os.path.join(results_root, "conditions.csv"), index=False)

    fig_dir = os.path.join(results_root, "figures")
    fig_condition_grid(conditions, fig_dir)
    fig_condition_effects(conditions, baseline, fig_dir)

    path = os.path.join(results_root, "OVERVIEW.md")
    with open(path, "w") as fh:
        fh.write(_render(camps, conds, baseline, results_root))
    return path


def _render(camps: pd.DataFrame, conds: pd.DataFrame, baseline, results_root: str) -> str:
    n_causal_lo = "sparse (1 causal SNP)"
    n_causal_hi = "polygenic (20 causal)"
    lines = [
        "# cis-heritability benchmark — all conditions",
        "",
        f"{len(camps)} campaigns, **{len(conds)} experimental conditions**, "
        f"{camps['cells'].sum():,} grid cells over {camps['regions'].sum():,} "
        "simulated cis-regions.",
        "",
        "Every swept level is its own row: a sample-size sweep contributes one row "
        "per *n* with its direction relative to the baseline, never a single "
        "averaged \"sample size\" entry.",
        "",
        "## Conditions",
        "",
        f"| group | condition | direction | cells | GCTA bias | SuSiE bias "
        f"| SuSiE bias, {n_causal_lo} | SuSiE bias, {n_causal_hi} | SuSiE wins |",
        "|---|---|---|---:|---:|---:|---:|---:|---:|",
    ]
    for _, r in conds.iterrows():
        mark = " **←**" if r["condition"] == baseline.label else ""
        lines.append(
            f"| {r['group']} | {r['condition']}{mark} | {r['direction']} | {r['cells']:,} "
            f"| {r['GCTA_bias']:+.4f} | {r['SuSiE_bias']:+.4f} "
            f"| {r['SuSiE_bias_sparse']:+.4f} | {r['SuSiE_bias_polygenic']:+.4f} "
            f"| {r['susie_win_rate']:.0%} |"
        )

    lines += [
        "",
        "![every swept level](figures/fig8_condition_grid.png)",
        "",
        "![effect of each change](figures/fig9_condition_effects.png)",
        "",
        "## Campaigns (provenance)",
        "",
        "| campaign | regions | cells | median p | axes swept | LD to SuSiE |",
        "|---|---:|---:|---:|---|---|",
    ]
    for _, r in camps.iterrows():
        lines.append(
            f"| [{r['campaign']}]({r['campaign']}/REPORT.md) | {r['regions']:,} "
            f"| {r['cells']:,} | {r['median_p']:,} | {r['axes']} | `{r['ld_for_susie']}` |"
        )
    lines += [
        "",
        "Each campaign has its own `REPORT.md`, `tables/*.csv` and `figures/*.pdf`. "
        "The narrative walk-through is in [FINDINGS.md](FINDINGS.md).",
        "",
        "---",
        "",
        "Generated by `python -m analysis.overview`.",
    ]
    return "\n".join(lines) + "\n"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("results_root", nargs="?",
                    default="/home/sux002/work/cis_heritability/results")
    args = ap.parse_args(argv)
    print(build_overview(args.results_root))
    return 0


if __name__ == "__main__":
    sys.exit(main())
