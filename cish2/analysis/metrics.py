"""Estimator-quality metrics with Monte-Carlo uncertainty.

Every metric is a mean or a proportion over replicate loci, so each carries the
Monte-Carlo standard error that the replicate count buys: ``sd/sqrt(N)`` for a
mean, ``sqrt(p(1-p)/N)`` for a proportion.  Reporting the MCSE alongside the
point estimate is what makes a bias difference between two methods claimable
rather than eyeballed.
"""

from __future__ import annotations

from typing import Dict, Optional, Sequence

import numpy as np
import pandas as pd

METRIC_COLUMNS = [
    "n_cells", "mean_h2_hat", "bias", "bias_mcse", "rel_bias",
    "rmse", "mae", "median_abs_err", "sd_h2_hat", "coverage95", "coverage95_mcse",
    "corr_truth", "power05", "power05_mcse", "frac_impossible", "max_abs_err",
]

# A heritability outside [0, 1] is not a heritability.  `--reml-no-constrain`
# lets REML return negatives, and SuSiE on a rank-deficient panel LD very
# occasionally converges (cleanly, reporting success) to an inflated estimate.
# Those few points dominate RMSE, so they are counted explicitly and MAE and
# the median absolute error are reported beside it.
H2_RANGE = (0.0, 1.0)


def mcse_mean(x: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    return float(np.std(x, ddof=1) / np.sqrt(x.size)) if x.size > 1 else np.nan


def mcse_prop(p: float, n: int) -> float:
    if not np.isfinite(p) or n < 1:
        return np.nan
    return float(np.sqrt(max(p * (1 - p), 0.0) / n))


def summarize(group: pd.DataFrame) -> pd.Series:
    """All metrics for one long-format group (one estimator, one design cell)."""
    h2_hat = group["h2_hat"].to_numpy(dtype=float)
    truth = group["h2_true"].to_numpy(dtype=float)
    err = h2_hat - truth
    finite = np.isfinite(err)
    err, h2_hat, truth = err[finite], h2_hat[finite], truth[finite]
    n = err.size

    bias = float(np.mean(err)) if n else np.nan
    cov = group["covers"].to_numpy(dtype=float)
    cov = cov[np.isfinite(cov)]
    coverage = float(np.mean(cov)) if cov.size else np.nan

    # "power": the estimate is significantly above zero at alpha = 0.05
    se = group["h2_se"].to_numpy(dtype=float)[finite]
    with np.errstate(divide="ignore", invalid="ignore"):
        signif = np.where(se > 0, h2_hat / se, np.nan) > 1.96
    signif = signif[np.isfinite(signif.astype(float))]
    power = float(np.mean(signif)) if signif.size else np.nan

    mean_truth = float(np.mean(truth)) if n else np.nan
    return pd.Series(
        {
            "n_cells": n,
            "mean_h2_hat": float(np.mean(h2_hat)) if n else np.nan,
            "bias": bias,
            "bias_mcse": mcse_mean(err),
            "rel_bias": bias / mean_truth if mean_truth else np.nan,
            "rmse": float(np.sqrt(np.mean(err**2))) if n else np.nan,
            "mae": float(np.mean(np.abs(err))) if n else np.nan,
            "median_abs_err": float(np.median(np.abs(err))) if n else np.nan,
            "max_abs_err": float(np.max(np.abs(err))) if n else np.nan,
            "frac_impossible": (
                float(np.mean((h2_hat < H2_RANGE[0]) | (h2_hat > H2_RANGE[1]))) if n else np.nan
            ),
            "sd_h2_hat": float(np.std(h2_hat, ddof=1)) if n > 1 else np.nan,
            "coverage95": coverage,
            "coverage95_mcse": mcse_prop(coverage, cov.size),
            "corr_truth": float(np.corrcoef(h2_hat, truth)[0, 1]) if n > 2 and np.std(truth) > 0 else np.nan,
            "power05": power,
            "power05_mcse": mcse_prop(power, signif.size),
        }
    )


def summarize_by(long_df: pd.DataFrame, by: Sequence[str]) -> pd.DataFrame:
    """Metrics for every (estimator x design) group."""
    keys = ["estimator"] + [c for c in by if c in long_df.columns]
    rows = []
    for key, part in long_df.groupby(keys, observed=True, dropna=False):
        rec = dict(zip(keys, key if isinstance(key, tuple) else (key,)))
        rec.update(summarize(part).to_dict())
        rows.append(rec)
    return pd.DataFrame(rows)


def win_rate(wide_df: pd.DataFrame, by: Sequence[str],
             estimators: Optional[Sequence[str]] = None) -> pd.DataFrame:
    """Share of cells where each estimator lands closest to the truth.

    Paired by construction: every estimator saw the same simulated phenotype,
    so the shares are over the same denominator and sum to 1.  ``susie_win_rate``
    is kept as an alias so the two-estimator narrative still reads.
    """
    if estimators is None:
        estimators = [e for e in ("GCTA", "SuSiE", "HEELS")
                      if f"abs_err_{e}" in wide_df.columns]
    keys = [c for c in by if c in wide_df.columns]
    grouped = wide_df.groupby(keys, observed=True, dropna=False) if keys else [((), wide_df)]
    rows = []
    for key, part in grouped:
        n = len(part)
        rec: Dict[str, object] = dict(zip(keys, key if isinstance(key, tuple) else (key,)))
        rec["n_cells"] = n
        for est in estimators:
            share = float((part["winner"] == est).mean()) if n else np.nan
            rec[f"win_{est}"] = share
            rec[f"win_{est}_mcse"] = mcse_prop(share, n)
        rec["susie_win_rate"] = rec.get("win_SuSiE", np.nan)
        rec["susie_win_mcse"] = rec.get("win_SuSiE_mcse", np.nan)
        rows.append(rec)
    return pd.DataFrame(rows)


def pairwise_win_rate(wide_df: pd.DataFrame, by: Sequence[str],
                      a: str, b: str) -> pd.DataFrame:
    """Head-to-head share where ``a`` beats ``b``, ignoring any third method."""
    keys = [c for c in by if c in wide_df.columns]
    grouped = wide_df.groupby(keys, observed=True, dropna=False) if keys else [((), wide_df)]
    rows = []
    for key, part in grouped:
        ea = part[f"abs_err_{a}"].to_numpy(dtype=float)
        eb = part[f"abs_err_{b}"].to_numpy(dtype=float)
        ok = np.isfinite(ea) & np.isfinite(eb)
        n = int(ok.sum())
        share = float((ea[ok] < eb[ok]).mean()) if n else np.nan
        rec: Dict[str, object] = dict(zip(keys, key if isinstance(key, tuple) else (key,)))
        rec.update(n_cells=n, win_rate=share, win_mcse=mcse_prop(share, n))
        rows.append(rec)
    return pd.DataFrame(rows)


def paired_test(wide_df: pd.DataFrame, by: Sequence[str]) -> pd.DataFrame:
    """Paired comparison of absolute error, per design cell.

    Both estimators see the *same* simulated phenotype, so the comparison is
    paired and a Wilcoxon signed-rank test is the right non-parametric call.
    """
    from scipy import stats

    keys = [c for c in by if c in wide_df.columns]
    rows = []
    groups = wide_df.groupby(keys, observed=True, dropna=False) if keys else [((), wide_df)]
    for key, part in groups:
        a = part["abs_err_GCTA"].to_numpy(dtype=float)
        b = part["abs_err_SuSiE"].to_numpy(dtype=float)
        ok = np.isfinite(a) & np.isfinite(b)
        a, b = a[ok], b[ok]
        rec: Dict[str, object] = dict(zip(keys, key if isinstance(key, tuple) else (key,)))
        if a.size < 8 or np.allclose(a, b):
            rec.update(n_pairs=int(a.size), median_diff=np.nan, pvalue=np.nan)
        else:
            stat, pval = stats.wilcoxon(a, b)
            rec.update(n_pairs=int(a.size),
                       median_diff=float(np.median(b - a)),  # SuSiE minus GCTA
                       pvalue=float(pval))
        rows.append(rec)
    return pd.DataFrame(rows)
