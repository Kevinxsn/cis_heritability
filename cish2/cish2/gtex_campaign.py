"""Driver for the GTEx real-data pass.

Loads the dataset once in the parent process so the 100 MB expression matrix is
shared with forked workers by copy-on-write, then fans chunks of genes across
the pool.  Same shard/manifest layout as the simulation campaigns, so the whole
analysis stack reads it unchanged.
"""

from __future__ import annotations

import logging
import os
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Any, Dict, List, Optional

import numpy as np
import yaml

from .gtex import GtexDataset, disjoint_subsamples, load_dataset, nested_subsamples
from .gtex_runner import ChunkOutcome, GeneWorker
from .results import ResultStore
from .scratch import describe, resolve_scratch

_LOG = logging.getLogger("cish2")
_WORKER: Optional[GeneWorker] = None

DEFAULTS: Dict[str, Any] = {
    "name": "gtex_whole_blood",
    "tissue": "Whole_Blood",
    "seed": 20260804,
    "results_root": "/home/sux002/work/cis_heritability/results",
    "geno_prefix": "/expanse/lustre/projects/ddp412/sux002/cis_heritability_data/gtex/"
                   "GTEx_v8_genotype_EUR_HM3_exclude_dups.allchr.reorder",
    "expression": "/expanse/lustre/projects/ddp412/sux002/cis_heritability_data/gtex/"
                  "Whole_Blood.v8.normalized_expression.bed.gz",
    "covariates": "/expanse/lustre/projects/ddp412/sux002/cis_heritability_data/gtex/Covar_all_Whole_Blood.txt",
    "chromosomes": list(range(1, 23)),
    "window_bp": 500_000,
    "min_snps": 50,
    "max_snps": 25_000,
    "n_individuals": [100, 200, 300, 450, 0],   # 0 means "all available"
    # "nested" = subsets of increasing size sharing donors (sample-size axis);
    # "split_half" = disjoint groups, for an honest test-retest reliability
    "subset_mode": "nested",
    "n_splits": 2,
    "max_genes": 0,                              # 0 means all
    "genes_per_chunk": 50,
    "susie_L": 10,
    "susie_max_iter": 100,
    "susie_coverage": 0.95,
    "clump_p": 1e-4,
    "clump_r2": 0.1,
    "gcta": "gcta64",
    "rscript": "Rscript",
    "gcta_reml_no_constrain": True,
    "run_heels": True,
    "heels_tol": 1e-5,
    "heels_max_iter": 200,
    "heels_init_h2": 0.25,
    "workers": 48,
    "task_timeout_s": 7200,
    "scratch_root": None,
    "description": "",
    # --- round-4 additions (PI feedback); all default-off -------------------
    # "real" = GTEx expression (every existing campaign).  "simulated" replaces
    # the expression with a phenotype simulated on the real genotypes of the
    # same donors; "permuted" keeps real expression and covariates but permutes
    # genotype rows within each subset, so every gene's true cis-h2 is 0.
    "phenotype_mode": "real",
    "sim_design": {
        "h2_grid": [0, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.4],
        "ncausal_grid": [1, 2, 3, 5, 10, 20],
        "null_mass": 0.25,          # P(h2 = 0); the rest uniform over the non-null grid
        "draws_per_gene": 2,
        "effect_dist": "normal",
    },
    # a second, unconstrained REML fit -> gcta_h2_uc / gcta_h2_uc_se / gcta_uc_status
    "gcta_also_unconstrained": False,
}

PHENOTYPE_MODES = ("real", "simulated", "permuted")


def load_config(path: str, overrides: Optional[List[str]] = None) -> Dict[str, Any]:
    cfg = dict(DEFAULTS)
    with open(path) as fh:
        cfg.update(yaml.safe_load(fh) or {})
    for ov in overrides or []:
        key, _, value = ov.partition("=")
        if key not in cfg:
            raise ValueError(f"unknown config key: {key}")
        cfg[key] = yaml.safe_load(value)
    unknown = set(cfg) - set(DEFAULTS)
    if unknown:
        raise ValueError(f"unknown config keys: {sorted(unknown)}")
    if cfg["phenotype_mode"] not in PHENOTYPE_MODES:
        raise ValueError(f"phenotype_mode must be one of {PHENOTYPE_MODES}")
    design = dict(DEFAULTS["sim_design"])
    design.update(cfg.get("sim_design") or {})
    bad = set(design) - set(DEFAULTS["sim_design"])
    if bad:
        raise ValueError(f"unknown sim_design keys: {sorted(bad)}")
    cfg["sim_design"] = design
    return cfg


def _init_worker(cfg: Dict[str, Any], data: GtexDataset, scratch_root: str) -> None:
    global _WORKER
    for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        os.environ[var] = "1"
    pid_scratch = os.path.join(scratch_root, f"w{os.getpid()}")
    os.makedirs(pid_scratch, exist_ok=True)
    _WORKER = GeneWorker(cfg, data, pid_scratch)


def _run_chunk(args) -> ChunkOutcome:
    chunk, rows = args
    assert _WORKER is not None
    return _WORKER.run_chunk(chunk, rows)


def run_gtex(cfg: Dict[str, Any], resume: bool = True) -> str:
    from .campaign import setup_logging

    run_dir = os.path.join(cfg["results_root"], cfg["name"])
    store = ResultStore(run_dir)
    setup_logging(os.path.join(store.log_dir, "run.log"))
    with open(os.path.join(run_dir, "config.yaml"), "w") as fh:
        yaml.safe_dump(cfg, fh, sort_keys=False)

    _LOG.info("loading GTEx dataset ...")
    data = load_dataset(cfg["geno_prefix"], cfg["expression"], cfg["covariates"],
                        cfg["chromosomes"])
    _LOG.info("%d samples with genotype + expression + covariates; %d genes; %d covariates",
              data.n_samples, len(data.genes), data.covar.shape[1])

    cfg = dict(cfg)
    if cfg["subset_mode"] == "split_half":
        cfg["subsets"] = disjoint_subsamples(data.n_samples, int(cfg["n_splits"]),
                                             cfg["seed"])
        _LOG.info("disjoint splits: %s (sizes %s) - no donor is in two splits",
                  sorted(cfg["subsets"]),
                  [len(v) for v in cfg["subsets"].values()])
    else:
        levels = [data.n_samples if int(x) in (0,) else int(x) for x in cfg["n_individuals"]]
        levels = sorted({min(x, data.n_samples) for x in levels})
        cfg["subsets"] = {n: idx for n, idx in
                          nested_subsamples(data.n_samples, levels, cfg["seed"]).items()}
        _LOG.info("sample-size levels: %s (nested subsets)", sorted(cfg["subsets"]))

    gene_rows = list(range(len(data.genes)))
    if cfg["max_genes"]:
        rng = np.random.default_rng(cfg["seed"])
        gene_rows = sorted(rng.choice(len(gene_rows), size=min(cfg["max_genes"], len(gene_rows)),
                                      replace=False).tolist())

    per_chunk = int(cfg["genes_per_chunk"])
    chunks = [(i, gene_rows[k:k + per_chunk])
              for i, k in enumerate(range(0, len(gene_rows), per_chunk))]

    done = store.completed_slots() if resume else set()
    todo = [c for c in chunks if c[0] not in done]

    scratch_root = resolve_scratch(explicit=cfg["scratch_root"],
                                   fallback=os.path.join(run_dir, "scratch"),
                                   tag=f"cish2-{cfg['name']}")
    _LOG.info("%d genes in %d chunks | %d already done, %d to run | workers=%d | scratch=%s",
              len(gene_rows), len(chunks), len(done), len(todo), cfg["workers"],
              describe(scratch_root))
    if not todo:
        _LOG.info("nothing to do")
        return run_dir

    t0 = time.time()
    n_ok = n_rows = n_fail = n_empty = 0
    with ProcessPoolExecutor(max_workers=min(cfg["workers"], len(todo)),
                             initializer=_init_worker,
                             initargs=(cfg, data, scratch_root)) as pool:
        futures = {pool.submit(_run_chunk, c): c[0] for c in todo}
        for i, fut in enumerate(as_completed(futures), start=1):
            chunk = futures[fut]
            try:
                outcome = fut.result()
            except Exception as exc:  # pragma: no cover
                _LOG.error("chunk %d crashed: %s", chunk, exc)
                n_fail += 1
                continue
            if outcome.status == "ok" and outcome.rows:
                store.write_shard(outcome.chunk, outcome.rows)
                n_ok += 1
                n_rows += len(outcome.rows)
            elif outcome.status == "ok":
                # every gene in the chunk was filtered out (too few SNPs in the
                # window, or no expression variance) - not an error
                n_empty += 1
                _LOG.info("chunk %d: no eligible genes", chunk)
            else:
                n_fail += 1
                _LOG.warning("chunk %d: %s %s", chunk, outcome.status, outcome.message[:300])
            store.append_manifest(dict(slot=outcome.chunk, status=outcome.status,
                                       n_rows=len(outcome.rows),
                                       seconds=round(outcome.seconds, 2),
                                       message=outcome.message[:300]))
            elapsed = time.time() - t0
            eta = (len(todo) - i) * elapsed / i
            _LOG.info("[%d/%d] chunk %d %s in %.1fs | %d rows | %.1f min elapsed, ETA %.1f min",
                      i, len(todo), chunk, outcome.status, outcome.seconds, n_rows,
                      elapsed / 60, eta / 60)

    _LOG.info("done: %d chunks ok, %d empty, %d failed, %d rows in %.1f min",
              n_ok, n_empty, n_fail, n_rows, (time.time() - t0) / 60)
    try:
        _LOG.info("consolidated -> %s", store.consolidate())
    except FileNotFoundError:
        _LOG.warning("no shards to consolidate")
    return run_dir
