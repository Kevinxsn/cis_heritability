"""Round-4 §7 text and its number pins (PI feedback, docs/pi_feedback_round4_design.md).

Which §7 the poster gets is decided by part_g's pre-registered rule 1, read
from tables/calibration_summary.json - not chosen by hand:

* calibration FAILED (the 2026-10-28 run: mean |sim - real| gap 0.0139 against
  a 0.010 limit) -> design §3.3 row 1: the simulation does not go on the poster.
  §7 keeps its round-3 figure; the text gains one sentence from the
  permutation null (model-free, so rule 1 does not apply to it), and §10
  gains the limitation the design prescribes.
* calibration passed -> not written: the run failed, and writing prose for an
  outcome that did not happen would be text nobody has checked.

The prose is typed by hand, as everywhere on the poster; check_sec7 re-derives
every number in it from the tables.
"""

from __future__ import annotations

import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def calibration_passed(tables: Path = HERE / "tables") -> bool:
    return bool(json.load(open(tables / "calibration_summary.json"))["calibration_passes"])


def sec7_blocks(Block, body_pt):
    if calibration_passed():
        raise SystemExit("calibration passed: write the §3.4 B1 text in round4_text.py first")
    return [
        Block("section", "7  A conditioning artefact"),
        Block("figure", path="P6_strata_compare.png", caption=(
            "Mean CHEWS − GCTA within joint strata; 50–9,572 genes per cell, <30 blank.")),
        Block("body",
              "Strata from the same donors manufacture much of the pattern (mean "
              "|A − B| = 0.023); (B) uses disjoint donors. **No Simpson's paradox is "
              "claimed.** With genotypes permuted (true h² = 0), GCTA and CHEWS **both "
              "average 0.013**, so a higher GCTA floor alone does not explain CHEWS's "
              "lower values below 0.05.",
              size=body_pt),
    ]


LIMITATION = ("Semi-synthetic GTEx simulation missed the real strata beyond its "
              "pre-set tolerance; §7 is disagreement, not error.")


def check_sec7(check, tables: Path):
    import pandas as pd

    print("\n§7  round 4: permutation null and calibration (perm_null.csv, calibration_summary.json)")
    null = pd.read_csv(tables / "perm_null.csv")
    allg = null[null["scope"].str.startswith("all genes")]
    for est in ("GCTA", "CHEWS"):
        for _, r in allg[allg.estimator == est].iterrows():
            check(f"permuted: {est} mean h2, {r['scope'][11:]}", r["mean"], 0.013, 0.0005)
    coup = pd.read_csv(tables / "gap_coupled_strata.csv", index_col=0)
    disj = pd.read_csv(tables / "gap_disjoint_strata.csv", index_col=0)
    check("mean |coupled - disjoint| gap", float((coup - disj).abs().stack().mean()), 0.023, 0.0005)
    s = json.load(open(tables / "calibration_summary.json"))
    check("calibration rule 1 passes (1 = yes)", float(s["calibration_passes"]), 0, 0)
    check("calibration: mean |sim - real| gap", s["mad_gap"], 0.0139, 0.0005)
