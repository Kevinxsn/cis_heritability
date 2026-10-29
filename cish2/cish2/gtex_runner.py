"""Per-chunk worker for the GTEx real-data pass.

The unit of work is a **chunk of genes** rather than one gene: at a HapMap3
cis-window (~400 SNPs) the actual fitting takes a fraction of a second, so R
startup would dominate if it were paid per gene.  Within a chunk each gene is
fitted at every sample-size level, and the whole chunk goes to R in one call.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import time
import traceback
import zlib
from dataclasses import dataclass
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd

from .estimators.gcta import GctaResult, GctaRunner
from .estimators.heels import HeelsResult, HeelsSolver
from .estimators.susie import SusieResult, _f, _unlink
from .genotypes import folded_maf
from .gtex import (
    GtexDataset,
    Residualizer,
    clump_independent_signals,
    full_rank_columns,
    impute_columns,
    select_covariates,
    standardize_columns,
)
from .simulate import sim_beta_empirical
from .sumstats import marginal_ols

R_GENES_SCRIPT = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "R", "susie_genes.R"
)


@dataclass
class ChunkOutcome:
    chunk: int
    rows: List[Dict[str, Any]]
    status: str
    seconds: float
    message: str = ""


class GeneWorker:
    """Fits both estimators for a chunk of genes at every sample-size level."""

    def __init__(self, cfg, data: GtexDataset, scratch_root: str):
        self.cfg = cfg
        self.data = data
        self.scratch_root = scratch_root
        self.gcta = GctaRunner(
            gcta_exe=cfg["gcta"],
            reml_no_constrain=cfg["gcta_reml_no_constrain"],
            threads=1,
        )
        # round 4: "real" (default, the existing behaviour), "simulated" or
        # "permuted" - see gtex_campaign.DEFAULTS
        self.mode = cfg.get("phenotype_mode", "real")
        self.gcta_uc = (GctaRunner(gcta_exe=cfg["gcta"], reml_no_constrain=True, threads=1)
                        if cfg.get("gcta_also_unconstrained") else None)
        # subsets are keyed by *label*, not size: split-half mode has two
        # groups of identical size that must stay distinguishable
        self.subsets: Dict[str, np.ndarray] = cfg["subsets"]
        self.residualizers = {
            label: Residualizer(data.covar[idx], len(idx), data.covar_names)
            for label, idx in self.subsets.items()
        }
        # REML gets the covariates as fixed effects rather than a pre-residualised
        # phenotype, so it accounts for the degrees of freedom they consume.  The
        # block must be full rank or GCTA cannot invert X' V^-1 X.
        self.qcovar = {}
        for label, idx in self.subsets.items():
            block = data.covar[np.ix_(idx, select_covariates(data.covar_names, len(idx)))]
            self.qcovar[label] = block[:, full_rank_columns(block)]

    # ------------------------------------------------------------------ run --
    def run_chunk(self, chunk: int, gene_rows: List[int]) -> ChunkOutcome:
        t0 = time.time()
        workdir = os.path.join(self.scratch_root, f"chunk{chunk:05d}")
        try:
            rows = self._run_chunk(chunk, gene_rows, workdir)
            outcome = ChunkOutcome(chunk, rows, "ok", 0.0)
        except Exception as exc:
            outcome = ChunkOutcome(chunk, [], "error", 0.0,
                                   f"{type(exc).__name__}: {exc}\n{traceback.format_exc(limit=4)}")
        finally:
            shutil.rmtree(workdir, ignore_errors=True)
        outcome.seconds = time.time() - t0
        return outcome

    def _run_chunk(self, chunk: int, gene_rows: List[int], workdir: str) -> List[Dict[str, Any]]:
        cfg = self.cfg
        os.makedirs(workdir, exist_ok=True)

        pending: List[Dict[str, Any]] = []   # one entry per (gene, n) cell
        gcta_out: Dict[str, GctaResult] = {}
        gcta_uc_out: Dict[str, GctaResult] = {}

        for row_idx in gene_rows:
            grow = self.data.genes.iloc[row_idx]
            win = self.data.window(grow, cfg["window_bp"])
            if win is None or win.n_snps < cfg["min_snps"]:
                continue
            if win.n_snps > cfg["max_snps"]:
                continue

            raw_geno = self.data.genotypes(win)          # all samples x p
            raw_expr = self.data.expr[win.expr_row, :]   # all samples

            if self.mode == "real":
                for label, idx in self.subsets.items():
                    cell = self._prepare_cell(win, raw_geno, raw_expr, label, idx, workdir)
                    if cell is None:
                        continue
                    gcta_out[cell["cell_id"]] = self._run_gcta(cell, workdir)
                    if self.gcta_uc is not None:
                        gcta_uc_out[cell["cell_id"]] = self._run_gcta(cell, workdir, self.gcta_uc)
                    pending.append(cell)
                continue

            # round-4 modes: same subsets, same covariates, same estimators;
            # only the phenotype (simulated) or the genotype rows (permuted) change
            if self.mode == "simulated":
                draws = self._simulate_draws(win, raw_geno)
            else:
                draws = [dict(tag=None, y=raw_expr, extra=dict(phenotype_mode="permuted"))]
            for draw in draws:
                for j, (label, idx) in enumerate(self.subsets.items()):
                    geno_idx = None
                    extra = dict(draw["extra"])
                    if self.mode == "permuted":
                        geno_idx = self._permuted_rows(win.gene_id, j, idx)
                    cell = self._prepare_cell(win, raw_geno, draw["y"], label, idx, workdir,
                                              geno_idx=geno_idx, tag=draw["tag"])
                    if cell is None:
                        continue
                    if self.mode == "simulated":
                        extra.update(self._realised_h2(draw, label, idx))
                    cell["extra"] = extra
                    gcta_out[cell["cell_id"]] = self._run_gcta(cell, workdir)
                    if self.gcta_uc is not None:
                        gcta_uc_out[cell["cell_id"]] = self._run_gcta(cell, workdir, self.gcta_uc)
                    pending.append(cell)

        susie_out = self._run_susie(pending, workdir)

        rows = []
        for cell in pending:
            row = self._row(cell,
                            gcta_out.get(cell["cell_id"], GctaResult.failed("not_run")),
                            susie_out.get(cell["cell_id"], SusieResult.failed("not_run")))
            # new columns only when a round-4 option is on, so a real-mode run
            # with an existing config writes exactly the columns it always did
            row.update(cell.get("extra", {}))
            if self.gcta_uc is not None:
                uc = gcta_uc_out.get(cell["cell_id"], GctaResult.failed("not_run"))
                row.update(gcta_h2_uc=uc.h2, gcta_h2_uc_se=uc.h2_se, gcta_uc_pval=uc.pval,
                           gcta_uc_status="ok" if uc.converged else (uc.message or "failed")[:120])
            rows.append(row)
        return rows

    # ------------------------------------------------- round-4 phenotypes --
    def _gene_rng(self, gene_id: str, stream: int) -> np.random.Generator:
        """Seeded by (campaign seed, gene, stream): results never depend on
        chunking or worker order, the same principle as the main campaign."""
        return np.random.default_rng([int(self.cfg["seed"]),
                                      zlib.crc32(str(gene_id).encode()), int(stream)])

    def _simulate_draws(self, win, raw_geno) -> List[Dict[str, Any]]:
        """K semi-synthetic phenotypes on the real genotypes of all donors.

        One beta per draw, shared by every subset: truth is a property of the
        gene, as in the real-data question.  The design prior is drawn per
        (gene, draw) - ``null_mass`` at h2 = 0, the rest uniform over the
        non-null (h2, ncausal) grid.
        """
        design = self.cfg["sim_design"]
        G, _ = impute_columns(raw_geno)
        Z_full, _, _ = standardize_columns(G)
        nonnull_h2 = [float(h) for h in design["h2_grid"] if float(h) > 0]
        ks = [int(k) for k in design["ncausal_grid"]]
        out = []
        for d in range(int(design["draws_per_gene"])):
            rng = self._gene_rng(win.gene_id, d)
            if rng.random() < float(design["null_mass"]):
                h2, k = 0.0, 0
            else:
                h2 = nonnull_h2[int(rng.integers(len(nonnull_h2)))]
                k = ks[int(rng.integers(len(ks)))]
            beta, causal = sim_beta_empirical(Z_full, k, h2, rng,
                                              effect_dist=design["effect_dist"])
            g = Z_full @ beta
            y = g + rng.normal(0.0, np.sqrt(1.0 - h2), g.size)
            out.append(dict(
                tag=f"d{d}", y=y, g=g,
                extra=dict(phenotype_mode="simulated", draw=d, h2_target=h2,
                           ncausal_target=k, n_causal_snps=int(causal.size),
                           sim_beta_hash=f"{zlib.crc32(beta.tobytes()):08x}"),
            ))
        return out

    def _realised_h2(self, draw, label, idx) -> Dict[str, float]:
        """Per-subset realised h2: raw, and on the covariate-residualised scale
        that the summary-statistic estimators actually see."""
        g, y = draw["g"][idx], draw["y"][idx]
        resid = self.residualizers[label]
        vy = float(np.var(y))
        gr, yr = resid(g), resid(y)
        vyr = float(np.var(yr))
        return dict(h2_sample=float(np.var(g)) / vy if vy > 0 else np.nan,
                    h2_resid=float(np.var(gr)) / vyr if vyr > 0 else np.nan)

    def _permuted_rows(self, gene_id: str, subset_ord: int, idx: np.ndarray) -> np.ndarray:
        """Genotype rows permuted *within* the subset: expression and covariates
        stay paired, the genotype-expression link is broken, true cis-h2 = 0."""
        rng = self._gene_rng(gene_id, 1_000_000 + subset_ord)
        return idx[rng.permutation(len(idx))]

    # ------------------------------------------------------------ one cell --
    def _prepare_cell(self, win, raw_geno, raw_expr, label, idx, workdir,
                      geno_idx=None, tag=None) -> Optional[Dict[str, Any]]:
        """Residualise, standardise, compute sumstats and LD for one (gene, subset).

        ``geno_idx`` (default ``idx``) lets the permutation null take genotype
        rows from different donors than the expression; ``tag`` keeps several
        simulated draws of one gene distinct.
        """
        if geno_idx is None:
            geno_idx = idx
        cfg = self.cfg
        resid = self.residualizers[label]
        n = len(idx)

        y = resid(raw_expr[idx].astype(np.float64))
        sd_y = y.std()
        if not np.isfinite(sd_y) or sd_y <= 0:
            return None
        y = (y - y.mean()) / sd_y

        # impute first, then residualise: the projection would turn any column
        # holding a single NaN into an all-NaN column and drop that SNP.
        # Allele frequency is taken from the raw dosages, before residualising,
        # where it still means something.
        G_raw, freq = impute_columns(raw_geno[geno_idx, :])
        maf_all = folded_maf(freq)

        # REML works on the *unadjusted* phenotype and GRM, with the covariates
        # supplied as fixed effects; the summary-statistic estimators work on the
        # residualised versions.  Each method is then used the way it expects.
        y_raw = raw_expr[idx].astype(np.float64)
        sd_raw = y_raw.std()
        if not np.isfinite(sd_raw) or sd_raw <= 0:
            return None
        y_raw = (y_raw - y_raw.mean()) / sd_raw
        Z_gcta, _, _ = standardize_columns(G_raw)

        # residualise the genotypes too: susie_rss wants the LD of the same
        # design the sumstats came from, covariates included
        Z, _, keep = standardize_columns(resid(G_raw))
        maf = maf_all[keep]
        p = int(Z.shape[1])
        if p < cfg["min_snps"]:
            return None

        beta, se, pval = marginal_ols(Z, y)
        R = (Z.T @ Z) / Z.shape[0]
        np.fill_diagonal(R, 1.0)

        # HEELS needs in-sample LD, which is exactly what real data provides
        # here: the genotypes and the expression come from the same people, so
        # both summary-statistic estimators see the same, correctly specified LD
        heels_res = HeelsResult.failed("not_run")
        if cfg.get("run_heels", True):
          try:
            # n_eff, not n: residualising on the covariates spent that many
            # degrees of freedom, and the likelihood has to know
            heels_res = HeelsSolver.from_genotypes(Z).solve(
                beta, len(idx),
                n_eff=len(idx) - resid.dof_used,
                tol=cfg.get("heels_tol", 1e-5),
                max_iter=cfg.get("heels_max_iter", 200),
                init_h2=cfg.get("heels_init_h2", 0.25),
            )
          except Exception as exc:  # never let one gene take down the chunk
            heels_res = HeelsResult.failed(f"{type(exc).__name__}: {exc}"[:120])

        cell_id = f"{win.gene_id}|{label}" + (f"|{tag}" if tag else "")
        path = os.path.join(workdir, f"{abs(hash(cell_id)) % (10**12)}_{n}.f64")
        with open(path, "wb") as fh:
            np.ascontiguousarray(R, dtype="<f8").tofile(fh)
            np.ascontiguousarray(beta, dtype="<f8").tofile(fh)
            np.ascontiguousarray(se, dtype="<f8").tofile(fh)

        n_clumps = clump_independent_signals(Z, pval, cfg["clump_p"], cfg["clump_r2"])
        chi2 = np.where(se > 0, (beta / np.where(se > 0, se, np.nan)) ** 2, np.nan)

        return dict(
            cell_id=cell_id,
            gene_id=win.gene_id,
            chrom=win.chrom,
            tss=win.tss,
            n=n,
            subset=str(label),
            p_snps=p,
            file=os.path.basename(path),
            Z=Z,
            y=y,
            maf_mean=float(np.mean(maf)),
            max_chisq=float(np.nanmax(chi2)) if np.isfinite(chi2).any() else np.nan,
            min_pval=float(np.nanmin(pval)),
            n_clumps=int(n_clumps),
            n_covar=resid.dof_used,
            heels=heels_res,
            Z_gcta=Z_gcta,
            y_gcta=y_raw,
            qcovar=self.qcovar.get(label),
        )

    def _run_gcta(self, cell: Dict[str, Any], workdir: str,
                  runner: Optional[GctaRunner] = None) -> GctaResult:
        runner = runner or self.gcta
        prefix = os.path.join(workdir, f"g{abs(hash(cell['cell_id'])) % (10**12)}")
        try:
            runner.prepare(cell["Z_gcta"], prefix)
            results = runner.run_batch(cell["y_gcta"].reshape(-1, 1), prefix,
                                       qcovar=cell.get("qcovar"))
        finally:
            runner.cleanup(prefix)
        return results[0] if results else GctaResult.failed("no result")

    # ---------------------------------------------------------- susie batch --
    def _run_susie(self, cells: List[Dict[str, Any]], workdir: str) -> Dict[str, SusieResult]:
        if not cells:
            return {}
        cfg = self.cfg
        pd.DataFrame({
            "cell_id": [c["cell_id"] for c in cells],
            "n": [c["n"] for c in cells],
            "p": [c["p_snps"] for c in cells],
            "file": [c["file"] for c in cells],
        }).to_csv(os.path.join(workdir, "cells.tsv"), sep="\t", index=False)

        params = {"susie_L": cfg["susie_L"], "susie_max_iter": cfg["susie_max_iter"],
                  "susie_coverage": cfg["susie_coverage"]}
        pd.DataFrame({"key": list(params), "value": [str(v) for v in params.values()]}).to_csv(
            os.path.join(workdir, "params.tsv"), sep="\t", index=False)

        env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
                   MKL_NUM_THREADS="1")
        try:
            proc = subprocess.run([cfg["rscript"], "--vanilla", R_GENES_SCRIPT, "--dir", workdir],
                                  capture_output=True, text=True, env=env,
                                  timeout=cfg["task_timeout_s"])
        except subprocess.TimeoutExpired:
            return {c["cell_id"]: SusieResult.failed("timeout") for c in cells}

        out_path = os.path.join(workdir, "susie.tsv")
        if not os.path.exists(out_path):
            tail = (proc.stderr or proc.stdout or "").strip().splitlines()
            msg = f"rscript_failed: {tail[-1] if tail else 'no output'}"
            return {c["cell_id"]: SusieResult.failed(msg[:200]) for c in cells}

        table = pd.read_csv(out_path, sep="\t")
        results: Dict[str, SusieResult] = {}
        for _, row in table.iterrows():
            results[str(row["cell_id"])] = SusieResult(
                h2=_f(row.get("susie_h2")), h2_se=_f(row.get("susie_h2_se")),
                h2_quad=_f(row.get("susie_h2_quad")), sum_mu2=_f(row.get("susie_sum_mu2")), sum_sd2=_f(row.get("susie_sum_sd2")),
                n_cs=_f(row.get("susie_n_cs")), cs_size=_f(row.get("susie_cs_size")),
                max_pip=_f(row.get("susie_max_pip")), sum_pip=_f(row.get("susie_sum_pip")),
                niter=_f(row.get("susie_niter")), sigma2=_f(row.get("susie_sigma2")),
                converged=bool(row.get("susie_converged", False)),
                status=str(row.get("susie_status", "unknown")),
            )
        for name in ("cells.tsv", "params.tsv", "susie.tsv"):
            _unlink(os.path.join(workdir, name))
        return results

    # ------------------------------------------------------------- one row --
    def _row(self, cell: Dict[str, Any], gcta: GctaResult, susie: SusieResult) -> Dict[str, Any]:
        heels: HeelsResult = cell.get("heels", HeelsResult.failed("not_run"))
        return dict(
            dataset=self.cfg["name"],
            tissue=self.cfg["tissue"],
            gene_id=cell["gene_id"],
            chrom=cell["chrom"],
            tss=cell["tss"],
            n=cell["n"],
            subset=cell["subset"],
            p_snps=cell["p_snps"],
            n_covar=cell["n_covar"],
            gcta_covar_mode="qcovar",
            maf_mean=cell["maf_mean"],
            max_chisq=cell["max_chisq"],
            min_pval=cell["min_pval"],
            n_clumps=cell["n_clumps"],
            gcta_h2=gcta.h2, gcta_h2_se=gcta.h2_se, gcta_vg=gcta.vg,
            gcta_vg_se=gcta.vg_se, gcta_ve=gcta.ve, gcta_vp=gcta.vp,
            gcta_pval=gcta.pval, gcta_logl=gcta.logl,
            gcta_converged=gcta.converged,
            gcta_status="ok" if gcta.converged else (gcta.message or "failed")[:120],
            susie_h2=susie.h2, susie_h2_se=susie.h2_se,
            susie_h2_quad=susie.h2_quad, susie_sum_mu2=susie.sum_mu2, susie_sum_sd2=susie.sum_sd2,
            susie_n_cs=susie.n_cs, susie_cs_size=susie.cs_size,
            susie_max_pip=susie.max_pip, susie_sum_pip=susie.sum_pip,
            susie_niter=susie.niter, susie_sigma2=susie.sigma2,
            susie_converged=susie.converged, susie_status=susie.status[:120],
            heels_h2=heels.h2, heels_h2_se=heels.h2_se, heels_s2g=heels.s2g,
            heels_s2e=heels.s2e, heels_niter=heels.niter, heels_loglik=heels.loglik,
            heels_converged=heels.converged, heels_status=heels.status[:120],
        )
