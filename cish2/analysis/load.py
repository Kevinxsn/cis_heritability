"""Load a campaign into tidy frames.

Two shapes are useful and both are produced here:

``wide``
    one row per grid cell, GCTA and SuSiE side by side -- the natural shape
    for paired comparisons (agreement scatter, Bland-Altman, win rate).
``long``
    one row per (cell, estimator) -- the natural shape for faceted plots and
    groupby metrics, where "estimator" is just another factor.
"""

from __future__ import annotations

import os
from typing import List

import numpy as np
import pandas as pd

import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from cish2.results import ResultStore  # noqa: E402

ESTIMATORS = {
    "GCTA": dict(h2="gcta_h2", se="gcta_h2_se", ok="gcta_converged"),
    "SuSiE": dict(h2="susie_h2", se="susie_h2_se", ok="susie_converged"),
    "HEELS": dict(h2="heels_h2", se="heels_h2_se", ok="heels_converged"),
}


def estimators_in(df: pd.DataFrame) -> List[str]:
    """Estimators actually present in a table.

    Campaigns run before HEELS was added have no `heels_*` columns; every
    consumer works off this list so old and new results both load.
    """
    return [name for name, cols in ESTIMATORS.items() if cols["h2"] in df.columns]

DESIGN_COLUMNS = ["n", "ncausal", "h2_true", "effect_dist", "maf_bin", "ld_arch"]

# These timings are recorded once per region and copied onto every row of that
# region; the per-cell figures divide them by the region's cell count so all
# ``t_*_s`` columns are comparable on one axis.
REGION_LEVEL_TIMINGS = ["t_region_draw_s", "t_ld_s"]


def per_cell_timings(df: pd.DataFrame) -> pd.DataFrame:
    """Return the ``t_*_s`` columns rescaled to a consistent per-cell basis."""
    cols = [c for c in df.columns if c.startswith("t_") and c.endswith("_s")]
    if not cols or "region_id" not in df.columns:
        return df[cols].copy() if cols else pd.DataFrame(index=df.index)
    out = df[cols].copy()
    cells_per_region = df.groupby("region_id", observed=True)["region_id"].transform("size")
    for col in REGION_LEVEL_TIMINGS:
        if col in out.columns:
            out[col] = out[col] / cells_per_region
    return out


def load_run(run_dir: str, drop_failed: bool = True) -> pd.DataFrame:
    """Load one campaign's cell table and attach derived error columns.

    :param drop_failed: drop cells where an estimator produced no usable
        number.  Failures are kept in ``attrition()`` either way, so nothing
        is silently discarded.
    """
    df = ResultStore(run_dir).load()
    df["run_dir"] = run_dir
    df["campaign"] = df.get("campaign", os.path.basename(run_dir))

    names = estimators_in(df)
    for name in names:
        cols = ESTIMATORS[name]
        h2 = df[cols["h2"]]
        df[f"err_{name}"] = h2 - df["h2_true"]
        df[f"abs_err_{name}"] = (h2 - df["h2_true"]).abs()
        # a z-score of the estimate against truth, for calibration/coverage
        se = df[cols["se"]].replace(0, np.nan)
        df[f"z_{name}"] = (h2 - df["h2_true"]) / se
        lo, hi = h2 - 1.96 * df[cols["se"]], h2 + 1.96 * df[cols["se"]]
        df[f"covers_{name}"] = (lo <= df["h2_true"]) & (df["h2_true"] <= hi)

    # the winner is whichever estimator lands closest to the truth on this
    # exact simulated phenotype - a paired comparison across all of them
    errs = df[[f"abs_err_{n}" for n in names]].to_numpy(dtype=float)
    with np.errstate(invalid="ignore"):
        best = np.full(len(df), "none", dtype=object)
        finite_any = np.isfinite(errs).any(axis=1)
        if finite_any.any():
            idx = np.nanargmin(np.where(np.isfinite(errs), errs, np.inf), axis=1)
            best[finite_any] = [names[i] for i in idx[finite_any]]
    df["winner"] = best

    if drop_failed:
        usable = np.ones(len(df), dtype=bool)
        for name in names:
            usable &= df[ESTIMATORS[name]["h2"]].notna().to_numpy()
        df = df.loc[usable].reset_index(drop=True)
    return df


def attrition(run_dir: str) -> pd.DataFrame:
    """Per-estimator failure accounting for the whole campaign."""
    df = ResultStore(run_dir).load()
    rows = []
    total = len(df)
    for name in estimators_in(df):
        cols = ESTIMATORS[name]
        missing = int(df[cols["h2"]].isna().sum())
        status_col = f"{name.lower()}_status"
        reasons = (
            df.loc[df[cols["h2"]].isna(), status_col].value_counts().to_dict()
            if status_col in df.columns else {}
        )
        rows.append(
            dict(estimator=name, cells=total, failed=missing,
                 failed_pct=100.0 * missing / total if total else np.nan,
                 reasons="; ".join(f"{k} ({v})" for k, v in list(reasons.items())[:3]))
        )
    if "attempts" in df.columns:
        rows.append(dict(estimator="retries", cells=total,
                         failed=int((df["attempts"] > 1).sum()),
                         failed_pct=100.0 * (df["attempts"] > 1).mean(),
                         reasons="cells needing >1 simulation attempt"))
    return pd.DataFrame(rows)


def to_long(df: pd.DataFrame) -> pd.DataFrame:
    """Reshape the wide cell table to one row per (cell, estimator)."""
    keep = [c for c in
            ["campaign", "region_id", "chrom", "p_snps", "p_susie", "ld_score",
             "slot", "h2_sample"] + DESIGN_COLUMNS if c in df.columns]
    frames: List[pd.DataFrame] = []
    names = estimators_in(df)
    for name in names:
        cols = ESTIMATORS[name]
        part = df[keep].copy()
        part["estimator"] = name
        part["h2_hat"] = df[cols["h2"]].to_numpy()
        part["h2_se"] = df[cols["se"]].to_numpy()
        part["error"] = df[f"err_{name}"].to_numpy()
        part["abs_error"] = df[f"abs_err_{name}"].to_numpy()
        part["covers"] = df[f"covers_{name}"].to_numpy()
        frames.append(part)
    out = pd.concat(frames, ignore_index=True)
    out["estimator"] = pd.Categorical(out["estimator"], categories=names, ordered=True)
    return out


def varying_axes(df: pd.DataFrame) -> List[str]:
    """Design columns that actually vary in this campaign."""
    return [c for c in DESIGN_COLUMNS if c in df.columns and df[c].nunique() > 1]
