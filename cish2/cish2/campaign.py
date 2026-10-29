"""Campaign driver: fan region slots across a process pool, stream results.

Region slots are independent, so the pool is a plain
``ProcessPoolExecutor`` with each worker pinned to one BLAS thread -- packing
many single-threaded tasks beats a few multi-threaded ones for this workload,
which is memory-bandwidth bound in its dense p x p kernels.

Results are streamed to disk as each slot finishes rather than collected at
the end, so a campaign that dies mid-flight loses nothing and resumes from its
manifest.
"""

from __future__ import annotations

import logging
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import List, Optional

from .config import CampaignConfig
from .results import ResultStore
from .runner import RegionOutcome, RegionWorker
from .scratch import describe, resolve_scratch

_WORKER: Optional[RegionWorker] = None
_LOG = logging.getLogger("cish2")


def _pin_threads(n: int = 1) -> None:
    for var in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
    ):
        os.environ[var] = str(n)


def _init_worker(cfg: CampaignConfig, scratch_root: str, cache_dir: str) -> None:
    global _WORKER
    _pin_threads(cfg.compute.threads_per_worker)
    pid_scratch = os.path.join(scratch_root, f"w{os.getpid()}")
    os.makedirs(pid_scratch, exist_ok=True)
    _WORKER = RegionWorker(cfg, pid_scratch, cache_dir)


def _run_slot(slot: int) -> RegionOutcome:
    assert _WORKER is not None, "worker not initialised"
    return _WORKER.run_slot(slot)


def setup_logging(log_path: str, verbose: bool = True) -> None:
    fmt = logging.Formatter("%(asctime)s %(levelname)-7s %(message)s", "%Y-%m-%d %H:%M:%S")
    _LOG.setLevel(logging.DEBUG)
    _LOG.handlers.clear()
    fh = logging.FileHandler(log_path)
    fh.setFormatter(fmt)
    _LOG.addHandler(fh)
    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(fmt)
    sh.setLevel(logging.INFO if verbose else logging.WARNING)
    _LOG.addHandler(sh)


def run_campaign(cfg: CampaignConfig, resume: bool = True) -> str:
    """Execute (or resume) a campaign.  Returns the run directory."""
    store = ResultStore(cfg.run_dir)
    setup_logging(os.path.join(store.log_dir, "run.log"))
    cfg.save(os.path.join(cfg.run_dir, "config.yaml"))

    scratch_root = resolve_scratch(
        explicit=cfg.compute.scratch_root,
        fallback=os.path.join(cfg.run_dir, "scratch"),
        tag=f"cish2-{cfg.name}",
    )
    cache_dir = os.path.join(cfg.results_root, ".cache", "bim_index")
    os.makedirs(cache_dir, exist_ok=True)

    done = store.completed_slots() if resume else set()
    slots: List[int] = [s for s in range(cfg.n_regions) if s not in done]

    _LOG.info("campaign %s (fingerprint %s)", cfg.name, cfg.fingerprint)
    _LOG.info(
        "grid: %d cells/region x %d regions = %d cells | %d already done, %d to run",
        cfg.grid.size, cfg.n_regions, cfg.total_cells, len(done), len(slots),
    )
    _LOG.info("workers=%d threads/worker=%d scratch=%s",
              cfg.compute.workers, cfg.compute.threads_per_worker, describe(scratch_root))
    if not slots:
        _LOG.info("nothing to do")
        return cfg.run_dir

    # warm the .bim position cache once, serially, so N workers do not race
    from .genotypes import build_chrom_index

    for chrom in cfg.data.chromosomes:
        build_chrom_index(cfg.data.geno_prefix, chrom, cache_dir)
    _LOG.info("bim index ready for %d chromosomes", len(cfg.data.chromosomes))

    t0 = time.time()
    n_ok = n_rows = n_fail = 0
    workers = min(cfg.compute.workers, len(slots))

    with ProcessPoolExecutor(
        max_workers=workers,
        initializer=_init_worker,
        initargs=(cfg, scratch_root, cache_dir),
    ) as pool:
        futures = {pool.submit(_run_slot, s): s for s in slots}
        for done_count, fut in enumerate(as_completed(futures), start=1):
            slot = futures[fut]
            try:
                outcome = fut.result()
            except Exception as exc:  # pragma: no cover - pool-level failure
                _LOG.error("slot %d crashed the worker: %s", slot, exc)
                n_fail += 1
                continue

            if outcome.status == "ok" and outcome.rows:
                store.write_shard(outcome.slot, outcome.rows)
                if outcome.region_meta:
                    store.append_region(outcome.region_meta)
                n_ok += 1
                n_rows += len(outcome.rows)
            else:
                n_fail += 1
                _LOG.warning("slot %d: %s %s", slot, outcome.status, outcome.message[:300])

            store.append_manifest(
                dict(slot=outcome.slot, status=outcome.status, n_rows=len(outcome.rows),
                     seconds=round(outcome.seconds, 2), message=outcome.message[:300])
            )

            elapsed = time.time() - t0
            rate = done_count / elapsed
            eta = (len(slots) - done_count) / rate if rate > 0 else float("nan")
            _LOG.info(
                "[%d/%d] slot %d %s in %.1fs | %d rows | %.1f min elapsed, ETA %.1f min",
                done_count, len(slots), slot, outcome.status, outcome.seconds,
                n_rows, elapsed / 60, eta / 60,
            )

    _LOG.info("done: %d regions ok, %d failed, %d rows in %.1f min",
              n_ok, n_fail, n_rows, (time.time() - t0) / 60)

    try:
        out = store.consolidate()
        _LOG.info("consolidated -> %s", out)
    except FileNotFoundError:
        _LOG.warning("no shards to consolidate")
    return cfg.run_dir
