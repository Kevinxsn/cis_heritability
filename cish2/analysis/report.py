"""Build the full analysis output for a campaign: tables, figures, report.

``build_report(run_dir)`` is the single entry point.  It writes

    <run_dir>/tables/*.csv     machine-readable metrics
    <run_dir>/figures/*.png    400 dpi review copies
    <run_dir>/figures/*.pdf    vector masters for the manuscript
    <run_dir>/REPORT.md        a readable summary with the headline numbers

and adapts to the campaign: only the design axes that actually vary get their
own figure.
"""

from __future__ import annotations

import os
from typing import Dict, List, Optional, Sequence

import pandas as pd

from . import distributions as D
from . import figures as F
from .load import attrition, estimators_in, load_run, per_cell_timings, to_long, varying_axes
from .metrics import paired_test, summarize_by, win_rate

CORE_AXES = ["ncausal", "h2_true"]


def build_report(run_dir: str, compare: Optional[Sequence[str]] = None) -> str:
    run_dir = os.path.abspath(run_dir)
    fig_dir = os.path.join(run_dir, "figures")
    tab_dir = os.path.join(run_dir, "tables")
    os.makedirs(fig_dir, exist_ok=True)
    os.makedirs(tab_dir, exist_ok=True)

    wide = load_run(run_dir)
    long = to_long(wide)
    axes = varying_axes(wide)
    extra_axes = [a for a in axes if a not in CORE_AXES]

    # ------------------------------------------------------------- tables --
    tables: Dict[str, pd.DataFrame] = {}
    tables["metrics_overall"] = summarize_by(long, [])
    tables["metrics_by_grid"] = summarize_by(long, CORE_AXES)
    for axis in extra_axes:
        tables[f"metrics_by_{axis}"] = summarize_by(long, [axis] + CORE_AXES)
    tables["win_rate_by_grid"] = win_rate(wide, CORE_AXES)
    tables["paired_tests_by_ncausal"] = paired_test(wide, ["ncausal"])
    tables["attrition"] = attrition(run_dir)
    tables["timing_by_region"] = _timing_table(wide)

    for name, table in tables.items():
        table.to_csv(os.path.join(tab_dir, f"{name}.csv"), index=False)

    # ------------------------------------------------------------ figures --
    written: List[str] = []
    written += F.fig_accuracy_grid(long, fig_dir)
    written += F.fig_bias_by_polygenicity(long, fig_dir)
    written += F.fig_rmse_landscape(long, wide, fig_dir)
    written += F.fig_calibration(long, fig_dir)
    written += F.fig_head_to_head(wide, fig_dir)
    for axis in extra_axes:
        written += F.fig_axis_effect(long, axis, fig_dir, wide_df=wide)
    # distribution views: means hide the spread, and spread is most of the story
    written += D.fig_estimate_boxplots(long, fig_dir)
    written += D.fig_error_boxplots(long, fig_dir)
    written += D.fig_relative_error_boxplots(long, fig_dir)
    for axis in extra_axes:
        written += D.fig_axis_boxplots(long, axis, fig_dir)
    written += F.fig_compute_profile(wide, fig_dir)

    if compare:
        # compare at the level of swept *conditions*, not whole campaigns, so a
        # sweep never collapses into one direction-less panel
        from .conditions import collect

        written += F.fig_condition_grid(collect([run_dir, *compare]), fig_dir)

    # ------------------------------------------------------------- report --
    path = os.path.join(run_dir, "REPORT.md")
    with open(path, "w") as fh:
        fh.write(_render(run_dir, wide, long, tables, extra_axes, written))
    return path


def _timing_table(wide: pd.DataFrame) -> pd.DataFrame:
    """Per-region mean seconds per grid cell, stage by stage."""
    timings = per_cell_timings(wide)
    if timings.empty:
        return pd.DataFrame()
    cols = list(timings.columns)
    frame = pd.concat([wide[["region_id", "p_snps"]], timings], axis=1)
    per_region = (frame.groupby(["region_id", "p_snps"], observed=True)[cols]
                  .mean().reset_index())
    per_region["t_total_cell_s"] = per_region[cols].sum(axis=1)
    return per_region.sort_values("p_snps")


def _render(run_dir, wide, long, tables, extra_axes, figures) -> str:
    overall = tables["metrics_overall"].set_index("estimator")
    grid = tables["metrics_by_grid"]
    name = os.path.basename(run_dir)

    lines = [
        f"# cis-heritability benchmark - `{name}`",
        "",
        f"- **cells analysed:** {len(wide):,} "
        f"across {wide['region_id'].nunique():,} independent cis-regions",
        f"- **region size:** median {wide['p_snps'].median():,.0f} SNPs "
        f"(range {wide['p_snps'].min():,.0f}-{wide['p_snps'].max():,.0f})",
        f"- **design axes varying:** {', '.join(varying_axes(wide)) or 'none'}",
        f"- **LD handed to SuSiE:** `{wide['ld_for_susie'].iloc[0]}`",
        "",
        "## Headline",
        "",
        "| estimator | mean bias | RMSE | 95% coverage | corr. with truth |",
        "|---|---:|---:|---:|---:|",
    ]
    for est in overall.index:
        row = overall.loc[est]
        lines.append(
            f"| {est} | {row['bias']:+.4f} ± {1.96 * row['bias_mcse']:.4f} "
            f"| {row['rmse']:.4f} | {row['coverage95']:.1%} | {row['corr_truth']:.3f} |"
        )

    ests = [e for e in ["GCTA", "SuSiE", "HEELS"] if e in set(estimators_in(wide))]
    lines += ["", "### Bias by causal-SNP count", "",
              "| causal SNPs | " + " | ".join(f"{e} bias" for e in ests)
              + " | " + " | ".join(f"{e} wins" for e in ests) + " |",
              "|---|" + "---:|" * (2 * len(ests))]
    by_nc = summarize_by(long, ["ncausal"])
    wins = win_rate(wide, ["ncausal"]).set_index("ncausal")
    for nc in sorted(by_nc["ncausal"].unique()):
        bias = " | ".join(
            f"{by_nc[(by_nc['estimator'] == e) & (by_nc['ncausal'] == nc)].iloc[0]['bias']:+.4f}"
            for e in ests)
        share = " | ".join(f"{wins.loc[nc, f'win_{e}']:.0%}" for e in ests)
        lines.append(f"| {nc} | {bias} | {share} |")

    att = tables["attrition"]
    lines += ["", "### Convergence and attrition", "",
              "| estimator | cells | failed | % | main reasons |", "|---|---:|---:|---:|---|"]
    for _, r in att.iterrows():
        lines.append(f"| {r['estimator']} | {r['cells']:,} | {r['failed']:,} "
                     f"| {r['failed_pct']:.2f}% | {r['reasons'] or '-'} |")

    timing = tables["timing_by_region"]
    if not timing.empty:
        total = timing["t_total_cell_s"].mean()
        lines += ["", "### Cost", "",
                  f"- mean **{total:.2f} s per grid cell** at median p = "
                  f"{wide['p_snps'].median():,.0f}",
                  f"- total simulated compute: **{len(wide) * total / 3600:.2f} core-hours**"]

    pngs = sorted({os.path.basename(p) for p in figures if p.endswith(".png")})
    lines += ["", "## Figures", ""]
    for png in pngs:
        lines.append(f"### {png.replace('.png', '').replace('_', ' ')}")
        lines.append(f"![{png}](figures/{png})")
        lines.append("")

    lines += ["## Tables", ""]
    for tname in sorted(tables):
        lines.append(f"- `tables/{tname}.csv`")
    lines += ["", "---", "",
              "Generated by `cish2` - see `cish2/README.md` for the pipeline description."]
    return "\n".join(lines) + "\n"
