"""LD matrices and their Cholesky factor.

Two matrices per region, both held in memory and reused by every grid cell:

``R_panel``
    the raw reference-panel correlation matrix -- what PLINK ``--r`` used to
    write out as a multi-GB text file, here a float64 array.
``R_sim``
    the ridge-regularised, renormalised version whose Cholesky factor is used
    to sample genotypes.

Computing both costs one GEMM; the legacy pipeline paid for the PLINK version
once per *grid cell* and shipped it through a text file.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.linalg as linalg


@dataclass
class RegionLD:
    R_panel: np.ndarray  # p x p, raw panel correlation
    R_sim: np.ndarray  # p x p, regularised (used for simulation)
    L: np.ndarray  # p x p, lower Cholesky factor of R_sim
    ridge: float

    @property
    def p(self) -> int:
        return self.R_sim.shape[0]


def build_ld(Z_ref: np.ndarray, ridge: float = 0.1, jitter_max: float = 1e-2) -> RegionLD:
    """Build the panel and simulation LD matrices plus the Cholesky factor.

    :param Z_ref: ``n_ref x p`` column-standardised reference genotypes.
    :param ridge: ridge added to the diagonal before renormalising, which is
        what makes the matrix comfortably positive definite.
    :param jitter_max: if the factorisation still fails, escalate a jitter up
        to this value rather than discarding the region.
    """
    n = Z_ref.shape[0]
    R_panel = (Z_ref.T @ Z_ref) / n
    np.fill_diagonal(R_panel, 1.0)

    R_sim = R_panel.copy()
    if ridge > 0:
        R_sim.flat[:: R_sim.shape[0] + 1] += ridge
        R_sim /= 1.0 + ridge

    jitter = 0.0
    while True:
        try:
            target = R_sim if jitter == 0.0 else R_sim + np.eye(R_sim.shape[0]) * jitter
            L = linalg.cholesky(target, lower=True)
            break
        except linalg.LinAlgError:
            jitter = 1e-8 if jitter == 0.0 else jitter * 10
            if jitter > jitter_max:
                raise
    return RegionLD(R_panel=R_panel, R_sim=R_sim, L=L, ridge=ridge)


def compute_s2g(L: np.ndarray, beta: np.ndarray) -> float:
    """Genetic variance ``beta' R beta`` via the Cholesky factor."""
    Ltb = L.T @ beta
    return float(Ltb @ Ltb)


def write_ld_binary(R: np.ndarray, path: str) -> None:
    """Dump a symmetric matrix as raw little-endian float64 for R's readBin.

    A p=8k region is 512 MB of binary here versus 2.8 GB of text from PLINK,
    and R reads it with one ``readBin`` instead of parsing 34M lines.
    """
    np.ascontiguousarray(R, dtype="<f8").tofile(path)
