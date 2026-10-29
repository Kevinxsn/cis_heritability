"""Genotype, effect-size and phenotype simulation.

Follows the twas_sim construction: sample genotypes from an MVN whose
covariance is the region LD, draw sparse causal effects, rescale them so the
realised genetic variance equals the target h2 exactly, then add Gaussian
noise to hit that h2 in the phenotype.
"""

from __future__ import annotations

from typing import Optional, Tuple

import numpy as np

from .ld import compute_s2g


def sim_geno(L: np.ndarray, n: int, rng: np.random.Generator) -> np.ndarray:
    """Sample ``n`` standardised genotypes from ``N(0, LL')``.

    :param L: ``p x p`` lower Cholesky factor of the region LD.
    :return: ``n x p`` column-standardised genotype matrix.
    """
    p = L.shape[0]
    Z = (L @ rng.standard_normal(size=(p, n))).T
    Z -= Z.mean(axis=0)
    sd = Z.std(axis=0)
    sd[sd == 0] = 1.0
    Z /= sd
    return Z


def _draw_effects(k: int, h2: float, rng: np.random.Generator,
                  effect_dist: str = "normal") -> np.ndarray:
    """Raw effect sizes for ``k`` causal SNPs, before any rescaling.

    Shared by :func:`sim_beta` and :func:`sim_beta_empirical`; the RNG calls are
    exactly those ``sim_beta`` always made, so existing campaigns reproduce.
    """
    if effect_dist == "normal":
        return rng.normal(loc=0.0, scale=np.sqrt(h2 / k), size=k)
    if effect_dist == "fixed":
        return rng.choice([-1.0, 1.0], size=k) * np.sqrt(h2 / k)
    if effect_dist == "pareto":
        return rng.choice([-1.0, 1.0], size=k) * rng.pareto(a=2.0, size=k)
    raise ValueError(f"unknown effect_dist: {effect_dist}")


def sim_beta(
    L: np.ndarray,
    ncausal: int,
    h2: float,
    rng: np.random.Generator,
    maf: Optional[np.ndarray] = None,
    maf_range: Tuple[float, float] = (0.0, 0.5),
    effect_dist: str = "normal",
    rescale: bool = True,
    replace: bool = False,
) -> Tuple[np.ndarray, np.ndarray]:
    """Draw causal effects under one architecture.

    :param maf: folded MAF per SNP; causal SNPs are restricted to
        ``maf_range`` when given (falls back to all SNPs if the region has too
        few variants in the bin).
    :param effect_dist: ``normal`` (infinitesimal-per-causal), ``fixed``
        (equal magnitude, random sign) or ``pareto`` (heavy tailed, a=2).
    :return: ``(beta, causal_idx)``.
    """
    p = L.shape[0]
    beta = np.zeros(p)
    if h2 <= 0 or ncausal <= 0:
        return beta, np.array([], dtype=int)

    pool = np.arange(p)
    if maf is not None and maf_range != (0.0, 0.5):
        lo, hi = maf_range
        candidates = np.flatnonzero((maf >= lo) & (maf <= hi))
        if candidates.size >= ncausal:
            pool = candidates
    causal = rng.choice(pool, size=min(ncausal, pool.size), replace=replace)

    effects = _draw_effects(causal.size, h2, rng, effect_dist)
    np.add.at(beta, causal, effects)

    if rescale:
        s2g = compute_s2g(L, beta)
        if s2g > 0:
            beta *= np.sqrt(h2 / s2g)
    return beta, np.unique(causal)


def sim_trait(g: np.ndarray, h2: float, rng: np.random.Generator) -> Tuple[np.ndarray, float]:
    """Add environmental noise to a genetic value to realise heritability ``h2``.

    :return: ``(y, h2_sample)`` - the standardised phenotype (mean 0, sd 1) and
        the heritability actually realised in this finite sample, which is the
        quantity an unbiased estimator should recover.
    """
    n = g.size
    if h2 > 0:
        s2g = np.var(g, ddof=1)
        s2e = s2g * ((1.0 / h2) - 1.0)
        y = g + rng.normal(0.0, np.sqrt(max(s2e, 0.0)), n)
    else:
        y = rng.normal(0.0, 1.0, n)
    var_y = np.var(y, ddof=1)
    h2_sample = float(np.var(g, ddof=1) / var_y) if var_y > 0 else np.nan
    y = y - y.mean()
    sd = y.std()
    return y / (sd if sd > 0 else 1.0), h2_sample


def sim_beta_empirical(
    Z: np.ndarray,
    ncausal: int,
    h2: float,
    rng: np.random.Generator,
    effect_dist: str = "normal",
) -> Tuple[np.ndarray, np.ndarray]:
    """Causal effects on *real* standardised genotypes (semi-synthetic GTEx).

    Same draw as :func:`sim_beta` - causal SNPs without replacement, effects from
    :func:`_draw_effects` - but there is no LD Cholesky factor here, so the
    rescaling is empirical: ``var(Z @ beta)`` equals ``h2`` exactly in the
    sample ``Z`` came from.

    :param Z: ``n x p`` column-standardised genotypes (all donors).
    :return: ``(beta, causal_idx)``; all-zero beta when ``h2 <= 0``.
    """
    p = Z.shape[1]
    beta = np.zeros(p)
    if h2 <= 0 or ncausal <= 0 or p == 0:
        return beta, np.array([], dtype=int)
    causal = rng.choice(p, size=min(ncausal, p), replace=False)
    beta[causal] = _draw_effects(causal.size, h2, rng, effect_dist)
    s2g = float(np.var(Z @ beta))
    if s2g > 0:
        beta *= np.sqrt(h2 / s2g)
    return beta, np.sort(causal)
