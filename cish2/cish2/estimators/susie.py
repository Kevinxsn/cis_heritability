"""Python side of the batched SuSiE fine-mapping bridge.

One ``Rscript`` per (region, sample size) instead of one per grid cell, and
the LD matrix crosses the boundary as raw float64 instead of a PLINK
long-format text file.  On a p = 8,000 region that replaces a 2.8 GB text
write + parse *per cell* with a single 512 MB binary write per region.
"""

from __future__ import annotations

import os
import subprocess
from dataclasses import dataclass, field
from typing import Any, List, Sequence

import numpy as np
import pandas as pd

R_SCRIPT = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "R", "susie_batch.R")

@dataclass
class SusieResult:
    h2: float = np.nan
    h2_se: float = np.nan
    h2_quad: float = np.nan
    sum_mu2: float = np.nan
    sum_sd2: float = np.nan
    n_cs: float = np.nan
    cs_size: float = np.nan
    max_pip: float = np.nan
    sum_pip: float = np.nan
    niter: float = np.nan
    sigma2: float = np.nan
    converged: bool = False
    status: str = "not_run"

    @classmethod
    def failed(cls, status: str) -> "SusieResult":
        return cls(status=status)


@dataclass
class SusieBatch:
    """Accumulates the per-cell sumstats of one region before one R call."""

    cell_ids: List[str] = field(default_factory=list)
    n_values: List[int] = field(default_factory=list)
    bhat: List[np.ndarray] = field(default_factory=list)
    shat: List[np.ndarray] = field(default_factory=list)

    def add(self, cell_id: str, n: int, beta: np.ndarray, se: np.ndarray) -> None:
        self.cell_ids.append(cell_id)
        self.n_values.append(int(n))
        self.bhat.append(np.asarray(beta, dtype=np.float64))
        self.shat.append(np.asarray(se, dtype=np.float64))

    def __len__(self) -> int:
        return len(self.cell_ids)


class SusieRunner:
    """Drives ``R/susie_batch.R`` over a whole region's grid at once."""

    def __init__(
        self,
        rscript: str = "Rscript",
        susie_L: int = 10,
        max_iter: int = 100,
        coverage: float = 0.95,
        threads: int = 1,
        timeout_s: int = 7200,
        script: str = R_SCRIPT,
    ):
        self.rscript = rscript
        self.susie_L = susie_L
        self.max_iter = max_iter
        self.coverage = coverage
        self.threads = threads
        self.timeout_s = timeout_s
        self.script = script

    def write_ld(self, R: np.ndarray, workdir: str) -> None:
        """Persist the region LD once; reused by every batch in the workdir."""
        os.makedirs(workdir, exist_ok=True)
        np.ascontiguousarray(R, dtype="<f8").tofile(os.path.join(workdir, "ld.f64"))

    def run(self, batch: SusieBatch, workdir: str, p: int) -> List[SusieResult]:
        if len(batch) == 0:
            return []
        os.makedirs(workdir, exist_ok=True)

        # R reads these back with matrix(v, nrow=p, ncol=K), i.e. column-major,
        # and ndarray.tofile always emits C order -- so stack cells as *rows*
        # here and each cell lands in its own R column.
        _write_cells(batch.bhat, os.path.join(workdir, "bhat.f64"))
        _write_cells(batch.shat, os.path.join(workdir, "shat.f64"))
        pd.DataFrame({"cell_id": batch.cell_ids, "n": batch.n_values}).to_csv(
            os.path.join(workdir, "cells.tsv"), sep="\t", index=False
        )
        params = {
            "p": p,
            "susie_L": self.susie_L,
            "susie_max_iter": self.max_iter,
            "susie_coverage": self.coverage,
        }
        pd.DataFrame({"key": list(params), "value": [str(v) for v in params.values()]}).to_csv(
            os.path.join(workdir, "params.tsv"), sep="\t", index=False
        )

        env = dict(
            os.environ,
            OMP_NUM_THREADS=str(self.threads),
            OPENBLAS_NUM_THREADS=str(self.threads),
            MKL_NUM_THREADS=str(self.threads),
        )
        try:
            proc = subprocess.run(
                [self.rscript, "--vanilla", self.script, "--dir", workdir],
                capture_output=True,
                text=True,
                env=env,
                timeout=self.timeout_s,
            )
        except subprocess.TimeoutExpired:
            return [SusieResult.failed("timeout") for _ in batch.cell_ids]

        out_path = os.path.join(workdir, "susie.tsv")
        if not os.path.exists(out_path):
            tail = (proc.stderr or proc.stdout or "").strip().splitlines()
            msg = f"rscript_failed: {tail[-1] if tail else 'no output'}"
            return [SusieResult.failed(msg[:200]) for _ in batch.cell_ids]

        table = pd.read_csv(out_path, sep="\t")
        by_id = {row["cell_id"]: row for _, row in table.iterrows()}
        results = []
        for cid in batch.cell_ids:
            row = by_id.get(cid)
            if row is None:
                results.append(SusieResult.failed("missing_from_r_output"))
                continue
            results.append(
                SusieResult(
                    h2=_f(row.get("susie_h2")),
                    h2_se=_f(row.get("susie_h2_se")),
                    h2_quad=_f(row.get("susie_h2_quad")),
                    sum_mu2=_f(row.get("susie_sum_mu2")),
                    sum_sd2=_f(row.get("susie_sum_sd2")),
                    n_cs=_f(row.get("susie_n_cs")),
                    cs_size=_f(row.get("susie_cs_size")),
                    max_pip=_f(row.get("susie_max_pip")),
                    sum_pip=_f(row.get("susie_sum_pip")),
                    niter=_f(row.get("susie_niter")),
                    sigma2=_f(row.get("susie_sigma2")),
                    converged=bool(row.get("susie_converged", False)),
                    status=str(row.get("susie_status", "unknown")),
                )
            )
        for name in ("bhat.f64", "shat.f64", "cells.tsv", "params.tsv", "susie.tsv"):
            _unlink(os.path.join(workdir, name))
        return results

    def cleanup(self, workdir: str) -> None:
        _unlink(os.path.join(workdir, "ld.f64"))


def _write_cells(vectors: Sequence[np.ndarray], path: str) -> None:
    """Write ``K`` length-p vectors as the column-major image of a p x K matrix."""
    np.ascontiguousarray(np.vstack(vectors), dtype="<f8").tofile(path)


def _f(value: Any) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return float("nan")
    return out


def _unlink(path: str) -> None:
    try:
        os.remove(path)
    except FileNotFoundError:
        pass
