"""GCTA/REML benchmark estimator.

Two changes from the legacy harness, both pure cost removal:

1. The GRM is written **once per (region, sample size)** rather than once per
   grid cell.  Every cell of a region reuses the same simulated genotype
   matrix, so its GRM is bit-identical across the 30+ cells of the grid.
2. The GRM is written in GCTA's **binary** format (float32 lower triangle)
   instead of gzipped text, which is ~15x less CPU and ~4x less disk for the
   same numbers.

All phenotypes of a region go into one multi-column .phen file selected with
``--mpheno``, so REML is the only per-cell cost left.
"""

from __future__ import annotations

import os
import subprocess
from dataclasses import dataclass
from typing import Dict, List, Optional

import numpy as np


@dataclass
class GctaResult:
    """One REML fit.  ``h2`` is V(G)/Vp; ``vg`` is the raw genetic variance."""

    h2: float = np.nan
    h2_se: float = np.nan
    vg: float = np.nan
    vg_se: float = np.nan
    ve: float = np.nan
    vp: float = np.nan
    pval: float = np.nan
    logl: float = np.nan
    converged: bool = False
    message: str = ""

    @classmethod
    def failed(cls, message: str) -> "GctaResult":
        return cls(message=message)


def write_grm_binary(grm: np.ndarray, n_snps: int, prefix: str, ids: Optional[List[str]] = None) -> None:
    """Write ``prefix.grm.bin`` / ``.grm.N.bin`` / ``.grm.id``.

    GCTA stores the lower triangle (including the diagonal) row by row as
    little-endian float32.
    """
    n = grm.shape[0]
    tril = grm[np.tril_indices(n)].astype("<f4")
    tril.tofile(f"{prefix}.grm.bin")
    np.full(tril.size, float(n_snps), dtype="<f4").tofile(f"{prefix}.grm.N.bin")
    ids = ids or [f"I{i}" for i in range(n)]
    with open(f"{prefix}.grm.id", "w") as fh:
        fh.write("".join(f"{i}\t{i}\n" for i in ids))


def write_qcovar(covar: np.ndarray, prefix: str, ids: Optional[List[str]] = None) -> None:
    """Write GCTA's quantitative-covariate file (headerless FID IID v1 v2 ...)."""
    n = covar.shape[0]
    ids = ids or [f"I{i}" for i in range(n)]
    with open(f"{prefix}.qcovar", "w") as fh:
        for i in range(n):
            vals = "\t".join(f"{v:.8g}" for v in covar[i])
            fh.write(f"{ids[i]}\t{ids[i]}\t{vals}\n")


def write_pheno_matrix(phenos: np.ndarray, prefix: str, ids: Optional[List[str]] = None) -> None:
    """Write an ``n x K`` phenotype matrix as one GCTA .phen file."""
    n, k = phenos.shape
    ids = ids or [f"I{i}" for i in range(n)]
    with open(f"{prefix}.phen", "w") as fh:
        for i in range(n):
            vals = "\t".join(f"{v:.8g}" for v in phenos[i])
            fh.write(f"{ids[i]}\t{ids[i]}\t{vals}\n")


def parse_hsq(path: str) -> GctaResult:
    """Parse a GCTA ``.hsq`` file by row label (robust to layout changes)."""
    values: Dict[str, List[float]] = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2 or parts[0] in ("Source",):
                continue
            try:
                nums = [float(x) for x in parts[1:] if x != ""]
            except ValueError:
                continue
            if nums:
                values[parts[0]] = nums

    def get(key: str, idx: int = 0) -> float:
        v = values.get(key)
        return float(v[idx]) if v is not None and len(v) > idx else np.nan

    return GctaResult(
        h2=get("V(G)/Vp"),
        h2_se=get("V(G)/Vp", 1),
        vg=get("V(G)"),
        vg_se=get("V(G)", 1),
        ve=get("V(e)"),
        vp=get("Vp"),
        pval=get("Pval"),
        logl=get("logL"),
        converged=True,
    )


class GctaRunner:
    """Runs REML for every cell of a region against a single shared GRM."""

    def __init__(self, gcta_exe: str = "gcta64", reml_no_constrain: bool = True, threads: int = 1):
        self.gcta_exe = gcta_exe
        self.reml_no_constrain = reml_no_constrain
        self.threads = threads

    def prepare(self, Z: np.ndarray, prefix: str) -> None:
        """Build and write the GRM for one simulated genotype matrix."""
        n, p = Z.shape
        grm = (Z @ Z.T) / p
        write_grm_binary(grm, p, prefix)

    def run_batch(self, phenos: np.ndarray, prefix: str,
                  qcovar: Optional[np.ndarray] = None) -> List[GctaResult]:
        """Fit REML for each column of ``phenos`` against the prepared GRM.

        :param qcovar: quantitative covariates to fit as fixed effects.  Passing
            them here rather than pre-residualising the phenotype is what lets
            REML account for the degrees of freedom they consume; the covariate
            block must be full rank or GCTA refuses to invert ``X' V^-1 X``.
        """
        write_pheno_matrix(phenos, prefix)
        if qcovar is not None and qcovar.size:
            write_qcovar(qcovar, prefix)
        results: List[GctaResult] = []
        env = dict(os.environ, OMP_NUM_THREADS=str(self.threads))
        for k in range(phenos.shape[1]):
            out = f"{prefix}.c{k}"
            cmd = [
                self.gcta_exe,
                "--grm", prefix,
                "--pheno", f"{prefix}.phen",
                "--mpheno", str(k + 1),
                "--reml",
                "--thread-num", str(self.threads),
                "--out", out,
            ]
            if qcovar is not None and qcovar.size:
                cmd += ["--qcovar", f"{prefix}.qcovar"]
            if self.reml_no_constrain:
                cmd.append("--reml-no-constrain")
            proc = subprocess.run(cmd, capture_output=True, text=True)
            hsq = f"{out}.hsq"
            if os.path.exists(hsq):
                results.append(parse_hsq(hsq))
            else:
                tail = (proc.stdout or proc.stderr or "").strip().splitlines()
                results.append(GctaResult.failed(tail[-1] if tail else "gcta produced no .hsq"))
            for suffix in (".hsq", ".log"):
                try:
                    os.remove(out + suffix)
                except FileNotFoundError:
                    pass
        for suffix in (".phen", ".qcovar"):
            try:
                os.remove(prefix + suffix)
            except FileNotFoundError:
                pass
        return results

    def cleanup(self, prefix: str) -> None:
        for suffix in (".grm.bin", ".grm.N.bin", ".grm.id", ".phen", ".qcovar"):
            try:
                os.remove(prefix + suffix)
            except FileNotFoundError:
                pass
