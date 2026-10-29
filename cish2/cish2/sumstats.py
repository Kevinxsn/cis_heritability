"""Marginal (per-SNP) association statistics, vectorised.

The legacy pipeline fitted ``statsmodels.OLS`` in a Python loop over all p
SNPs, which is p model objects and p QR factorisations per grid cell.  With
column-standardised genotypes and a standardised phenotype the whole set of
simple regressions collapses to a single matrix-vector product, so this is one
BLAS call for the entire region.
"""

from __future__ import annotations

from typing import Tuple

import numpy as np
from scipy import stats


def marginal_ols(Z: np.ndarray, y: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Per-SNP simple-regression estimates of ``y ~ 1 + z_j``.

    Assumes ``Z`` is column-centred/scaled and ``y`` centred/scaled, which is
    how :mod:`cish2.simulate` produces them.  Standard errors and p-values are
    identical to a per-SNP ``statsmodels.OLS`` fit with an intercept.

    :return: ``(beta, se, pval)``, each length p.
    """
    n, p = Z.shape
    if y.shape[0] != n:
        raise ValueError("Z and y disagree on sample size")

    # var(z_j) == 1 after standardisation, so beta_j = z_j'y / n
    szz = np.einsum("ij,ij->j", Z, Z)  # == n up to float error
    beta = (Z.T @ y) / szz

    syy = float(y @ y)
    rss = syy - beta**2 * szz
    dof = n - 2
    sigma2 = np.maximum(rss / dof, 0.0)
    se = np.sqrt(sigma2 / szz)

    with np.errstate(divide="ignore", invalid="ignore"):
        tstat = np.where(se > 0, beta / se, 0.0)
    pval = 2.0 * stats.t.sf(np.abs(tstat), dof)
    return beta, se, pval
