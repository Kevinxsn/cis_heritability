"""Per-region worker: draw a region, run the whole architecture grid on it.

The region is the unit of parallelism *and* the unit of caching.  Everything
that depends only on the region -- the .bed slice, the LD matrix, its Cholesky
factor, the binary LD handed to R -- is computed once and reused by every cell
of the grid.  Everything that depends only on (region, n) -- the simulated
genotype matrix and its GRM -- is computed once per sample size.  Only the
phenotype, the marginal sumstats, REML and SuSiE are genuinely per-cell.
"""

from __future__ import annotations

import os
import shutil
import time
import traceback
from dataclasses import dataclass
from typing import Any, Dict, List, Optional

import numpy as np

from .config import CampaignConfig, maf_range
from .estimators.gcta import GctaResult, GctaRunner
from .estimators.heels import HeelsResult, HeelsSolver
from .estimators.susie import SusieBatch, SusieResult, SusieRunner
from .genotypes import ChromIndexCache
from .ld import build_ld
from .regions import Region, draw_region
from .simulate import sim_beta, sim_geno, sim_trait
from .sumstats import marginal_ols

# strand-ambiguous allele pairs dropped by the legacy run_susie.R
_AMBIGUOUS = {("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")}
_ACGT = set("ACGT")


@dataclass
class RegionOutcome:
    slot: int
    rows: List[Dict[str, Any]]
    status: str
    seconds: float
    message: str = ""
    region_meta: Optional[Dict[str, Any]] = None


def susie_snp_mask(a0: np.ndarray, a1: np.ndarray, legacy_filter: bool) -> np.ndarray:
    """SNPs handed to the fine-mapper.

    Reproduces the legacy allele QC: drop strand-ambiguous pairs and any
    non-ACGT allele.  Because our sumstats and reference panel are the same
    .bim there is nothing to disambiguate, so this only ever removes usable
    SNPs -- keep it on solely to match historical runs.
    """
    if not legacy_filter:
        return np.ones(a0.size, dtype=bool)
    up0 = np.char.upper(a0.astype(str))
    up1 = np.char.upper(a1.astype(str))
    keep = np.array(
        [(x in _ACGT and y in _ACGT and (x, y) not in _AMBIGUOUS) for x, y in zip(up0, up1)],
        dtype=bool,
    )
    return keep


class RegionWorker:
    """Stateful worker; one per process, reused across region slots."""

    def __init__(self, cfg: CampaignConfig, scratch_root: str, cache_dir: str):
        self.cfg = cfg
        self.scratch_root = scratch_root
        self.index_cache = ChromIndexCache(cfg.data.geno_prefix, cache_dir)
        self.gcta = GctaRunner(
            gcta_exe=cfg.compute.gcta,
            reml_no_constrain=cfg.method.gcta_reml_no_constrain,
            threads=cfg.compute.threads_per_worker,
        )
        self.susie = SusieRunner(
            rscript=cfg.compute.rscript,
            susie_L=cfg.method.susie_L,
            max_iter=cfg.method.susie_max_iter,
            coverage=cfg.method.susie_coverage,
            threads=cfg.compute.threads_per_worker,
            timeout_s=cfg.compute.task_timeout_s,
        )

    # ------------------------------------------------------------------ run --
    def run_slot(self, slot: int) -> RegionOutcome:
        t0 = time.time()
        workdir = os.path.join(self.scratch_root, f"slot{slot:05d}")
        try:
            outcome = self._run_slot(slot, workdir)
        except Exception as exc:  # a bad region must not kill the campaign
            outcome = RegionOutcome(
                slot=slot,
                rows=[],
                status="error",
                seconds=time.time() - t0,
                message=f"{type(exc).__name__}: {exc}\n{traceback.format_exc(limit=4)}",
            )
        finally:
            if not self.cfg.compute.keep_scratch:
                shutil.rmtree(workdir, ignore_errors=True)
        outcome.seconds = time.time() - t0
        return outcome

    def _run_slot(self, slot: int, workdir: str) -> RegionOutcome:
        cfg = self.cfg
        os.makedirs(workdir, exist_ok=True)
        rng = np.random.default_rng(np.random.SeedSequence([cfg.seed, slot]))

        t = time.time()
        region = draw_region(cfg, rng, self.index_cache, slot)
        t_draw = time.time() - t
        if region is None:
            return RegionOutcome(slot, [], "no_region", 0.0, "no region matched the LD architecture")

        t = time.time()
        ld = build_ld(region.Z_ref, ridge=cfg.method.ld_regularization)
        t_ld = time.time() - t

        keep = susie_snp_mask(region.a0, region.a1, cfg.method.legacy_ambiguous_filter)
        R_susie = ld.R_panel if cfg.method.ld_for_susie == "panel" else ld.R_sim
        if not keep.all():
            R_susie = np.ascontiguousarray(R_susie[np.ix_(keep, keep)])
        p_susie = int(keep.sum())

        # From here only the Cholesky factor is needed for simulation, and
        # R_susie for the two summary-statistic estimators.  Release the other
        # p x p matrices *before* the eigendecomposition allocates its own, or
        # peak memory per worker is 5 x p^2 instead of 3 x p^2.
        ld.R_panel = None
        ld.R_sim = None

        if cfg.method.run_susie:
            t = time.time()
            self.susie.write_ld(R_susie, workdir)
            t_ld += time.time() - t

        del R_susie

        rows: List[Dict[str, Any]] = []
        base = dict(
            campaign=cfg.name,
            fingerprint=cfg.fingerprint,
            slot=slot,
            region_id=region.region_id,
            chrom=region.chrom,
            center_bp=region.center_bp,
            start_bp=region.start_bp,
            end_bp=region.end_bp,
            p_snps=region.n_snps,
            p_susie=p_susie,
            n_ref=region.n_ref,
            ld_arch=cfg.region.ld_arch,
            ld_score=region.ld_score,
            region_draws=region.draws,
            ld_for_susie=cfg.method.ld_for_susie,
            t_region_draw_s=t_draw,
            t_ld_s=t_ld,
        )

        cells = cfg.grid.cells()
        for n in cfg.grid.n_individuals:
            rows.extend(self._run_sample_size(region, ld, keep, n,
                                              [c for c in cells if c["n"] == n],
                                              workdir, rng, base))

        self.susie.cleanup(workdir)
        return RegionOutcome(slot, rows, "ok", 0.0, region_meta=region.meta())

    # -------------------------------------------------------- one sample size --
    def _run_sample_size(
        self,
        region: Region,
        ld,
        keep: np.ndarray,
        n: int,
        cells: List[Dict[str, Any]],
        workdir: str,
        rng: np.random.Generator,
        base: Dict[str, Any],
    ) -> List[Dict[str, Any]]:
        cfg = self.cfg
        if not cells:
            return []

        t = time.time()
        Z = sim_geno(ld.L, n, rng)
        t_geno = time.time() - t

        # HEELS needs the LD of the sample its summary statistics came from, so
        # its basis is built from the simulated genotypes, once per (region, n).
        # The SVD of an n x p matrix is far cheaper than eigendecomposing p x p.
        heels: Optional[HeelsSolver] = None
        t_heels_setup = 0.0
        if cfg.method.run_heels:
            t = time.time()
            heels = HeelsSolver.from_genotypes(Z[:, keep])
            t_heels_setup = time.time() - t

        grm_prefix = os.path.join(workdir, f"n{n}")
        t_gcta_prep = 0.0
        if cfg.method.run_gcta:
            t = time.time()
            self.gcta.prepare(Z, grm_prefix)
            t_gcta_prep = time.time() - t

        # per-cell state, refreshed on retry
        pending = list(range(len(cells)))
        draws: Dict[int, Dict[str, Any]] = {}
        gcta_out: Dict[int, GctaResult] = {}
        susie_out: Dict[int, SusieResult] = {}
        heels_out: Dict[int, HeelsResult] = {}
        attempts: Dict[int, int] = {i: 0 for i in pending}
        t_sim = t_ols = t_gcta = t_susie = t_heels = 0.0

        for attempt in range(1, max(cfg.method.max_retries, 1) + 1):
            if not pending:
                break
            t = time.time()
            for i in pending:
                draws[i] = self._simulate_cell(region, ld, Z, cells[i], rng)
                attempts[i] = attempt
            t_sim += time.time() - t

            if cfg.method.run_gcta:
                t = time.time()
                phenos = np.column_stack([draws[i]["y"] for i in pending])
                for i, res in zip(pending, self.gcta.run_batch(phenos, grm_prefix)):
                    gcta_out[i] = res
                t_gcta += time.time() - t

            # one set of marginal statistics feeds both summary-statistic
            # estimators, so they see identical input
            t = time.time()
            batch = SusieBatch()
            sumstats: Dict[int, np.ndarray] = {}
            for i in pending:
                beta, se, _ = marginal_ols(Z, draws[i]["y"])
                draws[i]["max_chisq"] = float(np.nanmax((beta / np.where(se > 0, se, np.nan)) ** 2))
                sumstats[i] = beta[keep]
                if cfg.method.run_susie:
                    batch.add(f"cell{i}", n, beta[keep], se[keep])
            t_ols += time.time() - t

            if heels is not None:
                t = time.time()
                for i in pending:
                    heels_out[i] = heels.solve(
                        sumstats[i], n,
                        tol=cfg.method.heels_tol,
                        max_iter=cfg.method.heels_max_iter,
                        init_h2=cfg.method.heels_init_h2,
                    )
                t_heels += time.time() - t

            if not cfg.method.run_susie:
                pending = []
                break

            t = time.time()
            results = self.susie.run(batch, workdir, p=int(keep.sum()))
            t_susie += time.time() - t
            for i, res in zip(pending, results):
                susie_out[i] = res

            pending = [i for i in pending if susie_out[i].status != "ok"]

        rows = []
        for i, cell in enumerate(cells):
            rows.append(
                self._row(
                    base,
                    cell,
                    draws.get(i, {}),
                    gcta_out.get(i, GctaResult.failed("not_run")),
                    susie_out.get(i, SusieResult.failed("not_run")),
                    heels_out.get(i, HeelsResult.failed("not_run")),
                    attempts.get(i, 0),
                    dict(
                        t_geno_s=t_geno / len(cells),
                        t_gcta_prep_s=t_gcta_prep / len(cells),
                        t_sim_s=t_sim / len(cells),
                        t_ols_s=t_ols / len(cells),
                        t_gcta_s=t_gcta / len(cells),
                        t_susie_s=t_susie / len(cells),
                        t_heels_s=(t_heels + t_heels_setup) / len(cells),
                    ),
                )
            )

        if cfg.method.run_gcta:
            self.gcta.cleanup(grm_prefix)
        return rows

    # ------------------------------------------------------------- one cell --
    def _simulate_cell(self, region, ld, Z: np.ndarray, cell: Dict[str, Any], rng) -> Dict[str, Any]:
        cfg = self.cfg
        beta, causal = sim_beta(
            ld.L,
            cell["ncausal"],
            cell["h2_true"],
            rng,
            maf=region.maf,
            maf_range=maf_range(cell["maf_bin"]),
            effect_dist=cell["effect_dist"],
            rescale=cfg.method.rescale_beta,
            replace=cfg.method.causal_replace,
        )
        g = Z @ beta
        y, h2_sample = sim_trait(g, cell["h2_true"], rng)
        return dict(y=y, h2_sample=h2_sample, n_causal_drawn=int(causal.size))

    @staticmethod
    def _row(base, cell, draw, gcta: GctaResult, susie: SusieResult,
             heels: HeelsResult, attempts: int, timings) -> Dict[str, Any]:
        row = dict(base)
        row.update(
            rep=int(cell.get("rep", 1)),
            n=cell["n"],
            ncausal=cell["ncausal"],
            h2_true=cell["h2_true"],
            effect_dist=cell["effect_dist"],
            maf_bin=cell["maf_bin"],
            n_causal_drawn=draw.get("n_causal_drawn", np.nan),
            h2_sample=draw.get("h2_sample", np.nan),
            max_chisq=draw.get("max_chisq", np.nan),
            attempts=attempts,
            gcta_h2=gcta.h2,
            gcta_h2_se=gcta.h2_se,
            gcta_vg=gcta.vg,
            gcta_vg_se=gcta.vg_se,
            gcta_ve=gcta.ve,
            gcta_vp=gcta.vp,
            gcta_pval=gcta.pval,
            gcta_logl=gcta.logl,
            gcta_converged=gcta.converged,
            gcta_status="ok" if gcta.converged else (gcta.message or "failed")[:120],
            susie_h2=susie.h2,
            susie_h2_se=susie.h2_se,
            susie_h2_quad=susie.h2_quad,
            susie_sum_mu2=susie.sum_mu2,
            susie_sum_sd2=susie.sum_sd2,
            susie_n_cs=susie.n_cs,
            susie_cs_size=susie.cs_size,
            susie_max_pip=susie.max_pip,
            susie_sum_pip=susie.sum_pip,
            susie_niter=susie.niter,
            susie_sigma2=susie.sigma2,
            susie_converged=susie.converged,
            susie_status=susie.status[:120],
            heels_h2=heels.h2,
            heels_h2_se=heels.h2_se,
            heels_s2g=heels.s2g,
            heels_s2e=heels.s2e,
            heels_niter=heels.niter,
            heels_loglik=heels.loglik,
            heels_converged=heels.converged,
            heels_status=heels.status[:120],
        )
        row.update(timings)
        return row
