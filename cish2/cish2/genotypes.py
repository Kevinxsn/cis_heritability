"""Direct PLINK1 .bed/.bim readers.

The legacy pipeline shelled out to PLINK to cut every cis-region into its own
three-file set, then read it back with ``pandas_plink``.  That is one
subprocess plus three files per region.  Here we slice the chromosome-level
.bed in place: a variant-major .bed stores ``ceil(n/4)`` bytes per variant, so
a contiguous SNP range is a contiguous byte range and can be read with a
single seek.  No temp files, no subprocess, no PLINK dependency.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

_BED_MAGIC = b"\x6c\x1b\x01"

# PLINK1 2-bit codes -> dosage of the SECOND bim allele (column 6, "a1" in
# pandas_plink), so values match ``pandas_plink.read_plink`` exactly.  Code 1
# is the missing marker.
_CODE_TO_DOSAGE = np.array([0.0, np.nan, 1.0, 2.0], dtype=np.float32)


@dataclass(frozen=True)
class ChromIndex:
    """Positions of every variant on one chromosome, plus file geometry."""

    chrom: int
    prefix: str
    pos: np.ndarray  # int64, ascending
    n_samples: int

    @property
    def n_variants(self) -> int:
        return int(self.pos.size)


def _bed_path(prefix: str, chrom: int) -> str:
    return f"{prefix}{chrom}.bed"


def count_samples(prefix: str, chrom: int) -> int:
    with open(f"{prefix}{chrom}.fam") as fh:
        return sum(1 for _ in fh)


def build_chrom_index(prefix: str, chrom: int, cache_dir: Optional[str] = None) -> ChromIndex:
    """Parse (and cache) the variant positions of one chromosome.

    The .bim is the only large text file we ever parse; caching it as .npy
    turns a ~30 MB text scan into a memory-mapped read.
    """
    n_samples = count_samples(prefix, chrom)
    cache = None
    if cache_dir:
        os.makedirs(cache_dir, exist_ok=True)
        tag = os.path.basename(prefix.rstrip(".")) or "geno"
        cache = os.path.join(cache_dir, f"{tag}.chr{chrom}.pos.npy")
        if os.path.exists(cache):
            pos = np.load(cache, mmap_mode="r")
            return ChromIndex(chrom=chrom, prefix=prefix, pos=np.asarray(pos), n_samples=n_samples)

    pos = pd.read_csv(
        f"{prefix}{chrom}.bim",
        sep=r"\s+",
        header=None,
        usecols=[3],
        dtype={3: np.int64},
    )[3].to_numpy()

    if cache:
        tmp = f"{cache}.tmp{os.getpid()}.npy"
        np.save(tmp, pos)
        os.replace(tmp, cache)
    return ChromIndex(chrom=chrom, prefix=prefix, pos=pos, n_samples=n_samples)


def read_bim_rows(prefix: str, chrom: int, i0: int, i1: int) -> pd.DataFrame:
    """Read .bim rows ``[i0, i1)`` without materialising the whole file."""
    rows: List[List[str]] = []
    with open(f"{prefix}{chrom}.bim") as fh:
        for idx, line in enumerate(fh):
            if idx < i0:
                continue
            if idx >= i1:
                break
            rows.append(line.split())
    df = pd.DataFrame(rows, columns=["chrom", "snp", "cm", "pos", "a0", "a1"])
    df["pos"] = df["pos"].astype(np.int64)
    df["chrom"] = df["chrom"].astype(str)
    return df.reset_index(drop=True)


def read_bed_slice(
    prefix: str,
    chrom: int,
    i0: int,
    i1: int,
    n_samples: int,
) -> np.ndarray:
    """Variant range ``[i0, i1)`` from a per-chromosome fileset ``<prefix><chrom>``."""
    return read_bed_range(_bed_path(prefix, chrom), i0, i1, n_samples)


def read_bed_range(bed_path: str, i0: int, i1: int, n_samples: int) -> np.ndarray:
    """Return an ``n_samples x (i1-i0)`` float32 dosage matrix (NaN = missing).

    Dosages count the second .bim allele, matching ``pandas_plink.read_plink``.
    Works on any variant-major PLINK1 .bed, per-chromosome or whole-genome.
    """
    if i1 <= i0:
        raise ValueError("empty variant range")
    bytes_per_variant = (n_samples + 3) // 4
    n_var = i1 - i0
    offset = 3 + i0 * bytes_per_variant
    with open(bed_path, "rb") as fh:
        magic = fh.read(3)
        if magic != _BED_MAGIC:
            raise ValueError(f"{bed_path} is not a variant-major PLINK1 .bed")
        fh.seek(offset)
        raw = np.frombuffer(fh.read(bytes_per_variant * n_var), dtype=np.uint8)
    if raw.size != bytes_per_variant * n_var:
        raise IOError("short read from .bed - file truncated or index out of range")

    raw = raw.reshape(n_var, bytes_per_variant)
    # unpack 4 two-bit codes per byte, least-significant pair first
    codes = np.empty((n_var, bytes_per_variant * 4), dtype=np.uint8)
    for k in range(4):
        codes[:, k::4] = (raw >> (2 * k)) & 0b11
    codes = codes[:, :n_samples]
    return _CODE_TO_DOSAGE[codes].T.copy()  # samples x variants


def standardize(G: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Mean-impute, centre and scale genotypes.

    :return: ``(Z, maf, keep_mask)`` where ``Z`` is column-standardised over
        the variants that carry any variation and ``maf`` is the alternate
        allele frequency of the retained variants.
    """
    G = np.asarray(G, dtype=np.float64)
    col_mean = np.nanmean(G, axis=0)
    missing = np.isnan(G)
    if missing.any():
        G[missing] = np.take(col_mean, np.where(missing)[1])
    maf_raw = col_mean / 2.0
    sd = G.std(axis=0)
    keep = sd > 0
    Z = (G[:, keep] - col_mean[keep]) / sd[keep]
    return Z, maf_raw[keep], keep


def folded_maf(freq: np.ndarray) -> np.ndarray:
    """Fold an allele frequency onto [0, 0.5]."""
    return np.minimum(freq, 1.0 - freq)


class ChromIndexCache:
    """Per-process cache of chromosome indices (workers touch a few each)."""

    def __init__(self, prefix: str, cache_dir: Optional[str] = None):
        self.prefix = prefix
        self.cache_dir = cache_dir
        self._store: Dict[int, ChromIndex] = {}

    def get(self, chrom: int) -> ChromIndex:
        if chrom not in self._store:
            self._store[chrom] = build_chrom_index(self.prefix, chrom, self.cache_dir)
        return self._store[chrom]
