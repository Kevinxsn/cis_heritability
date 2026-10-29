"""Real-data layer: GTEx expression, covariates and cis-windows.

The simulation half of this project can measure *bias*, because it knows the
truth.  Real data cannot.  What it can test is whether the **patterns** the
simulations predict actually show up:

* the GCTA-minus-SuSiE gap should shrink as sample size grows — testable here
  by **subsampling individuals** from one tissue (nested subsamples, so the
  comparison stays paired across the sample-size axis);
* the gap should be near zero for genes whose signal is a single concentrated
  association, and large for genes whose signal is diffuse;
* fine-mapping should go out of range rarely but wildly, REML often but mildly.

Two departures from the previous GTEx script, both deliberate:

1. **Covariates are actually applied.**  The legacy script accepted a
   ``--covariate_file`` argument and never used it, so PEER factors, genotype
   PCs, sex, platform and PCR batch were left in the phenotype.  Here both the
   expression *and* the genotypes are residualised on the covariates, which is
   what ``susie_rss`` asks for when the sumstats were covariate-adjusted.
2. **LD is in-sample.**  Genotypes and expression come from the same
   individuals, so the LD handed to the fine-mapper is the exact LD of the data
   that produced the sumstats — no reference-panel mismatch.
"""

from __future__ import annotations

import gzip
import os
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from .genotypes import folded_maf, read_bed_range


@dataclass
class GeneWindow:
    """One gene's cis-window, resolved to a contiguous .bed variant range."""

    gene_id: str
    chrom: int
    tss: int
    start_bp: int
    end_bp: int
    i0: int
    i1: int
    expr_row: int

    @property
    def n_snps(self) -> int:
        return self.i1 - self.i0


@dataclass
class GtexDataset:
    """Everything shared across genes, built once and inherited by workers."""

    bed_path: str
    n_fam: int
    sample_ids: List[str]
    geno_rows: np.ndarray            # index of each sample into .fam order
    expr: np.ndarray                 # genes x samples, already sample-aligned
    genes: pd.DataFrame              # gene_id, chrom, tss, expr_row
    covar: np.ndarray                # samples x k, no intercept column
    covar_names: List[str]
    bim_chrom: np.ndarray
    bim_pos: np.ndarray
    chrom_bounds: Dict[int, Tuple[int, int]] = field(default_factory=dict)

    @property
    def n_samples(self) -> int:
        return len(self.sample_ids)

    def window(self, row: pd.Series, window_bp: int) -> Optional[GeneWindow]:
        """Resolve one gene's TSS +/- window to a .bed variant range."""
        chrom = int(row["chrom"])
        bounds = self.chrom_bounds.get(chrom)
        if bounds is None:
            return None
        lo, hi = bounds
        pos = self.bim_pos[lo:hi]
        tss = int(row["tss"])
        i0 = lo + int(np.searchsorted(pos, max(0, tss - window_bp)))
        i1 = lo + int(np.searchsorted(pos, tss + window_bp + 1))
        if i1 <= i0:
            return None
        return GeneWindow(gene_id=str(row["gene_id"]), chrom=chrom, tss=tss,
                          start_bp=max(0, tss - window_bp), end_bp=tss + window_bp,
                          i0=i0, i1=i1, expr_row=int(row["expr_row"]))

    def genotypes(self, win: GeneWindow) -> np.ndarray:
        """Dosages for one window, restricted and ordered to our samples."""
        G = read_bed_range(self.bed_path, win.i0, win.i1, self.n_fam)
        return G[self.geno_rows, :]


# ---------------------------------------------------------------------------
def _read_expression_header(path: str) -> List[str]:
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        return fh.readline().rstrip("\n").split("\t")


def load_dataset(
    geno_prefix: str,
    expression_path: str,
    covariate_path: Optional[str] = None,
    chromosomes: Sequence[int] = tuple(range(1, 23)),
) -> GtexDataset:
    """Load genotypes' index, expression and covariates onto a common sample set."""
    fam = pd.read_csv(f"{geno_prefix}.fam", sep=r"\s+", header=None,
                      usecols=[1], names=["iid"], dtype={1: str})
    fam_ids = fam["iid"].tolist()

    header = _read_expression_header(expression_path)
    expr_ids = header[4:]

    covar_df = None
    if covariate_path:
        covar_df = pd.read_csv(covariate_path, sep="\t")
        id_col = "ID" if "ID" in covar_df.columns else covar_df.columns[1]
        covar_df = covar_df.set_index(id_col)
        covar_df = covar_df.select_dtypes(include=[np.number])
        covar_ids = set(covar_df.index)
    else:
        covar_ids = set(expr_ids)

    fam_set = set(fam_ids)
    keep = [s for s in expr_ids if s in fam_set and s in covar_ids]
    if not keep:
        raise ValueError("no samples shared between genotypes, expression and covariates")

    fam_pos = {s: i for i, s in enumerate(fam_ids)}
    geno_rows = np.array([fam_pos[s] for s in keep], dtype=np.int64)
    expr_pos = {s: i for i, s in enumerate(expr_ids)}
    expr_cols = np.array([expr_pos[s] for s in keep], dtype=np.int64)

    frame = pd.read_csv(expression_path, sep="\t")
    frame.columns = ["chrom", "start", "end", "gene_id"] + expr_ids
    frame["chrom"] = (frame["chrom"].astype(str).str.replace("chr", "", regex=False))
    frame = frame[frame["chrom"].isin({str(c) for c in chromosomes})].reset_index(drop=True)

    expr = frame[expr_ids].to_numpy(dtype=np.float64)[:, expr_cols]
    genes = pd.DataFrame({
        "gene_id": frame["gene_id"].astype(str),
        "chrom": frame["chrom"].astype(int),
        "tss": frame["start"].astype(np.int64),
        "expr_row": np.arange(len(frame), dtype=np.int64),
    })

    covar = (covar_df.loc[keep].to_numpy(dtype=np.float64)
             if covar_df is not None else np.zeros((len(keep), 0)))
    covar_names = list(covar_df.columns) if covar_df is not None else []

    bim = pd.read_csv(f"{geno_prefix}.bim", sep=r"\s+", header=None,
                      usecols=[0, 3], names=["chrom", "pos"])
    bim_chrom = bim["chrom"].to_numpy(dtype=np.int64)
    bim_pos = bim["pos"].to_numpy(dtype=np.int64)
    bounds: Dict[int, Tuple[int, int]] = {}
    for chrom in np.unique(bim_chrom):
        idx = np.flatnonzero(bim_chrom == chrom)
        bounds[int(chrom)] = (int(idx[0]), int(idx[-1]) + 1)

    return GtexDataset(
        bed_path=f"{geno_prefix}.bed",
        n_fam=len(fam_ids),
        sample_ids=keep,
        geno_rows=geno_rows,
        expr=expr,
        genes=genes,
        covar=covar,
        covar_names=covar_names,
        bim_chrom=bim_chrom,
        bim_pos=bim_pos,
        chrom_bounds=bounds,
    )


# ---------------------------------------------------------------------------
def gtex_n_inferred_factors(n: int) -> int:
    """GTEx v8's own rule for how many inferred expression factors to use.

    The covariate count has to scale with sample size or a subsampled analysis
    is over-parameterised: at n = 100, using all 60 inferred factors plus 5
    genotype PCs plus intercept consumes 69 of 100 degrees of freedom and REML
    stops converging.  These are the thresholds the GTEx v8 eQTL pipeline uses.
    """
    if n < 150:
        return 15
    if n < 250:
        return 30
    if n < 350:
        return 45
    return 60


def select_covariates(names: Sequence[str], n: int) -> List[int]:
    """Column indices of the covariates to use at sample size ``n``.

    Always: genotype PCs and the technical covariates (sex, platform, PCR).
    Scaled: the inferred expression factors, per :func:`gtex_n_inferred_factors`.
    """
    keep_inferred = gtex_n_inferred_factors(n)
    idx = []
    for j, name in enumerate(names):
        low = str(name).lower()
        if low.startswith("inferredcov"):
            try:
                rank = int(low.replace("inferredcov", ""))
            except ValueError:
                rank = 10**6
            if rank <= keep_inferred:
                idx.append(j)
        else:
            idx.append(j)
    return idx


def full_rank_columns(block: np.ndarray, tol: float = 1e-8) -> np.ndarray:
    """Indices of a maximal linearly independent set of covariate columns.

    GTEx's covariate matrix is rank-deficient as shipped (at least one constant
    column), and GCTA refuses to run when ``X' V^-1 X`` is singular.  A
    pivoted QR picks a full-rank subset without changing the column space.
    """
    block = np.asarray(block, dtype=np.float64)
    if block.size == 0:
        return np.zeros(0, dtype=int)
    keep = np.flatnonzero(block.std(axis=0) > tol)
    if keep.size == 0:
        return keep
    sub = block[:, keep]
    import scipy.linalg as sla

    _, r, piv = sla.qr(sub, mode="economic", pivoting=True)
    diag = np.abs(np.diag(r))
    rank = int(np.sum(diag > tol * max(diag[0], 1.0)))
    return np.sort(keep[piv[:rank]])


class Residualizer:
    """Projects covariates (plus an intercept) out of phenotypes and genotypes.

    Built once per sample subset: the projection is ``I - QQ'`` for the thin QR
    factor ``Q`` of ``[1, C]``, so applying it to a p-column genotype block is
    one small GEMM rather than p separate regressions.
    """

    def __init__(self, covar: np.ndarray, n: int, names: Optional[Sequence[str]] = None):
        design = np.ones((n, 1))
        if covar is not None and covar.size:
            block = np.asarray(covar, dtype=np.float64)
            if names is not None:
                block = block[:, select_covariates(names, n)]
            # drop covariates that are constant in this subset
            keep = block.std(axis=0) > 0
            if keep.any():
                design = np.column_stack([design, block[:, keep]])
        self.Q, _ = np.linalg.qr(design)
        self.n_covar = design.shape[1]

    def __call__(self, X: np.ndarray) -> np.ndarray:
        X = np.asarray(X, dtype=np.float64)
        return X - self.Q @ (self.Q.T @ X)

    @property
    def dof_used(self) -> int:
        return self.n_covar


def impute_columns(X: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Fill missing calls with the column mean.

    This has to happen *before* covariate residualisation: the projection
    ``X - Q(Q'X)`` turns a column holding even one NaN into an all-NaN column,
    which would silently drop that SNP.  GTEx missingness is sparse (order one
    SNP per window) but not zero.
    """
    X = np.asarray(X, dtype=np.float64)
    missing = np.isnan(X)
    if not missing.any():
        return X, np.nanmean(X, axis=0) / 2.0
    col_mean = np.nanmean(X, axis=0)
    col_mean = np.where(np.isfinite(col_mean), col_mean, 0.0)
    X = X.copy()
    X[missing] = np.take(col_mean, np.where(missing)[1])
    return X, col_mean / 2.0


def standardize_columns(X: np.ndarray, min_sd: float = 1e-8
                        ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Centre and scale; drop columns with no variation.  Input must be complete."""
    X = np.asarray(X, dtype=np.float64)
    col_mean = X.mean(axis=0)
    sd = X.std(axis=0)
    keep = sd > min_sd
    Z = (X[:, keep] - col_mean[keep]) / sd[keep]
    return Z, folded_maf(col_mean[keep] / 2.0), keep


def nested_subsamples(n_total: int, levels: Sequence[int], seed: int
                      ) -> Dict[int, np.ndarray]:
    """Nested random subsets: the n=100 set is inside the n=200 set, and so on.

    Nesting keeps the sample-size axis paired — a change between levels is the
    added individuals, not a different draw of people.
    """
    rng = np.random.default_rng(seed)
    order = rng.permutation(n_total)
    out: Dict[int, np.ndarray] = {}
    for level in sorted(set(int(x) for x in levels)):
        size = min(level, n_total)
        out[size] = np.sort(order[:size])
    return out


def disjoint_subsamples(n_total: int, n_splits: int, seed: int
                         ) -> Dict[str, np.ndarray]:
    """Partition the donors into ``n_splits`` **disjoint** groups.

    Nested subsamples share individuals, so an estimate from n = 100 and one
    from n = 574 are correlated by construction and any "consistency" between
    them is inflated.  Disjoint halves give a genuine test-retest: no donor
    contributes to both estimates, so the correlation between them is an honest
    reliability coefficient.
    """
    rng = np.random.default_rng(seed)
    order = rng.permutation(n_total)
    return {f"split{i+1}": np.sort(chunk)
            for i, chunk in enumerate(np.array_split(order, n_splits))}


def clump_independent_signals(Z: np.ndarray, pvals: np.ndarray,
                              p_threshold: float = 1e-4,
                              r2_threshold: float = 0.1) -> int:
    """Count independent marginal signals by greedy LD clumping.

    A SuSiE-independent proxy for how many distinct associations a gene has —
    needed because stratifying SuSiE's own error by SuSiE's own credible-set
    count would be circular.
    """
    n, p = Z.shape
    order = np.argsort(pvals)
    order = order[pvals[order] < p_threshold]
    if order.size == 0:
        return 0
    available = np.ones(p, dtype=bool)
    clumps = 0
    for idx in order:
        if not available[idx]:
            continue
        clumps += 1
        r = (Z.T @ Z[:, idx]) / n
        available &= (r**2 < r2_threshold)
        available[idx] = False
    return clumps
