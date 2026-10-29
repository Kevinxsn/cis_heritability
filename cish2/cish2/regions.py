"""Random cis-region sampling with LD-architecture rejection.

A "region" is a +/- ``window_bp`` window around a randomly chosen variant.
The three LD architectures of the study (``regular`` / ``high`` / ``low``) are
implemented as an accept-reject rule on a cheap LD score computed from a SNP
subsample, so a rejected region costs a .bed slice rather than a full p x p
correlation matrix.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict
from typing import Any, Dict, Optional

import numpy as np

from .config import CampaignConfig
from .genotypes import ChromIndexCache, folded_maf, read_bed_slice, read_bim_rows, standardize


@dataclass
class Region:
    """A drawn cis-region and its standardised reference genotypes."""

    region_id: str
    chrom: int
    center_bp: int
    start_bp: int
    end_bp: int
    i0: int
    i1: int
    n_snps: int
    n_ref: int
    ld_score: float
    draws: int
    Z_ref: np.ndarray  # n_ref x p, column standardised
    maf: np.ndarray  # p, folded
    a0: np.ndarray  # p, first bim allele
    a1: np.ndarray  # p, second bim allele
    snp: np.ndarray  # p, rsIDs

    def meta(self) -> Dict[str, Any]:
        d = asdict(self)
        for key in ("Z_ref", "maf", "a0", "a1", "snp"):
            d.pop(key)
        return d


def ld_score(Z: np.ndarray, r2_threshold: float, max_snps: int, rng: np.random.Generator) -> float:
    """Fraction of SNP pairs with r^2 >= threshold, on a SNP subsample.

    Subsampling keeps this O(max_snps^2) instead of O(p^2); the fraction is an
    unbiased estimate of the region-wide statistic under random subsampling.
    """
    n, p = Z.shape
    m = min(max_snps, p)
    if m < 2:
        return 0.0
    sub = Z if p <= m else Z[:, rng.choice(p, size=m, replace=False)]
    R = (sub.T @ sub) / n
    iu = np.triu_indices(R.shape[0], k=1)
    r2 = R[iu] ** 2
    return float(np.mean(r2 >= r2_threshold)) if r2.size else 0.0


def _accepts(arch: str, score: float, cfg: CampaignConfig) -> bool:
    if arch == "regular":
        return True
    if arch == "high":
        return score >= cfg.region.frac_min
    if arch == "low":
        return score <= cfg.region.frac_max
    raise ValueError(arch)


def draw_region(
    cfg: CampaignConfig,
    rng: np.random.Generator,
    index_cache: ChromIndexCache,
    slot: int,
) -> Optional[Region]:
    """Draw one region matching the configured LD architecture.

    :return: an accepted :class:`Region`, or ``None`` if ``region.max_draws``
        candidates were all rejected.
    """
    data, reg = cfg.data, cfg.region
    for draw in range(1, reg.max_draws + 1):
        chrom = int(rng.choice(data.chromosomes))
        idx = index_cache.get(chrom)
        v = int(rng.integers(0, idx.n_variants))
        center = int(idx.pos[v])
        start, end = max(0, center - data.window_bp), center + data.window_bp

        i0, i1 = np.searchsorted(idx.pos, [start, end + 1])
        i0, i1 = int(i0), int(i1)
        if i1 - i0 < data.min_snps:
            continue
        if i1 - i0 > data.max_snps:
            continue

        G = read_bed_slice(data.geno_prefix, chrom, i0, i1, idx.n_samples)
        Z, freq, keep = standardize(G)
        if Z.shape[1] < data.min_snps:
            continue

        score = ld_score(Z, reg.r2_threshold, reg.score_max_snps, rng)
        if not _accepts(reg.ld_arch, score, cfg):
            continue

        bim = read_bim_rows(data.geno_prefix, chrom, i0, i1)
        bim = bim.loc[keep].reset_index(drop=True)
        return Region(
            region_id=f"{reg.ld_arch}_chr{chrom}_{center}_s{slot:05d}",
            chrom=chrom,
            center_bp=center,
            start_bp=start,
            end_bp=end,
            i0=i0,
            i1=i1,
            n_snps=int(Z.shape[1]),
            n_ref=int(Z.shape[0]),
            ld_score=score,
            draws=draw,
            Z_ref=Z,
            maf=folded_maf(freq),
            a0=bim["a0"].to_numpy(),
            a1=bim["a1"].to_numpy(),
            snp=bim["snp"].to_numpy(),
        )
    return None
