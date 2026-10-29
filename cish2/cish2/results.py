"""Result layout, sharded writes and resumability.

Every campaign owns one directory::

    results/<name>/
        config.yaml         exact configuration this run was launched with
        manifest.jsonl      one line per finished region slot (the resume log)
        regions.tsv         one row per drawn region (chrom, window, p, LD score)
        shards/*.tsv.gz     one gzipped TSV per region slot, written on completion
        cells.tsv.gz        consolidated long table, one row per grid cell
        logs/run.log        driver log
        tables/, figures/   analysis outputs

Sharding by slot is what makes the campaign restartable: a slot already in
``manifest.jsonl`` is skipped, so an interrupted run resumes at zero cost and
two runs of the same config never duplicate rows.
"""

from __future__ import annotations

import glob
import gzip
import json
import os
from typing import Any, Dict, List, Optional, Set

import pandas as pd

# stable column order: identity, design, truth, estimators, diagnostics
COLUMN_ORDER = [
    "campaign", "fingerprint", "slot", "region_id", "chrom", "center_bp",
    "start_bp", "end_bp", "p_snps", "p_susie", "n_ref", "ld_arch", "ld_score",
    "region_draws", "ld_for_susie",
    "n", "subset", "rep", "ncausal", "h2_true", "effect_dist", "maf_bin", "gcta_covar_mode",
    "n_causal_drawn", "h2_sample", "max_chisq", "attempts",
    "gcta_h2", "gcta_h2_se", "gcta_vg", "gcta_vg_se", "gcta_ve", "gcta_vp",
    "gcta_pval", "gcta_logl", "gcta_converged", "gcta_status",
    "susie_h2", "susie_h2_se", "susie_h2_quad", "susie_sum_mu2", "susie_sum_sd2", "susie_n_cs",
    "susie_cs_size", "susie_max_pip", "susie_sum_pip", "susie_niter",
    "susie_sigma2", "susie_converged", "susie_status",
    "heels_h2", "heels_h2_se", "heels_s2g", "heels_s2e", "heels_niter",
    "heels_loglik", "heels_converged", "heels_status",
    "t_region_draw_s", "t_ld_s", "t_geno_s", "t_gcta_prep_s",
    "t_sim_s", "t_ols_s", "t_gcta_s", "t_susie_s", "t_heels_s",
]


class ResultStore:
    """Writes and reads one campaign's results directory."""

    def __init__(self, run_dir: str):
        self.run_dir = run_dir
        self.shard_dir = os.path.join(run_dir, "shards")
        self.manifest_path = os.path.join(run_dir, "manifest.jsonl")
        self.regions_path = os.path.join(run_dir, "regions.tsv")
        self.log_dir = os.path.join(run_dir, "logs")
        for d in (self.run_dir, self.shard_dir, self.log_dir,
                  os.path.join(run_dir, "tables"), os.path.join(run_dir, "figures")):
            os.makedirs(d, exist_ok=True)

    # ------------------------------------------------------------- resume --
    def completed_slots(self) -> Set[int]:
        """Slots that need no rerun.

        A slot that *errored* is deliberately not counted: re-running the
        campaign should retry it, which is the whole point of the manifest
        surviving a crash.  Slots that finished with no rows (every gene or
        region filtered out) are complete and stay skipped.
        """
        if not os.path.exists(self.manifest_path):
            return set()
        done: Set[int] = set()
        failed: Set[int] = set()
        with open(self.manifest_path) as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                try:
                    record = json.loads(line)
                    slot = int(record["slot"])
                except (json.JSONDecodeError, KeyError, ValueError):
                    continue
                if str(record.get("status", "ok")) == "ok":
                    done.add(slot)
                else:
                    failed.add(slot)
        return done - failed

    # -------------------------------------------------------------- write --
    def write_shard(self, slot: int, rows: List[Dict[str, Any]]) -> Optional[str]:
        if not rows:
            return None
        df = pd.DataFrame(rows)
        ordered = [c for c in COLUMN_ORDER if c in df.columns]
        extra = [c for c in df.columns if c not in ordered]
        df = df[ordered + extra]
        path = os.path.join(self.shard_dir, f"slot{slot:05d}.tsv.gz")
        tmp = path + ".tmp"
        df.to_csv(tmp, sep="\t", index=False, compression="gzip")
        os.replace(tmp, path)
        return path

    def append_manifest(self, record: Dict[str, Any]) -> None:
        with open(self.manifest_path, "a") as fh:
            fh.write(json.dumps(record, default=str) + "\n")
            fh.flush()

    def append_region(self, meta: Dict[str, Any]) -> None:
        header = not os.path.exists(self.regions_path)
        pd.DataFrame([meta]).to_csv(
            self.regions_path, sep="\t", index=False, header=header, mode="a"
        )

    # --------------------------------------------------------- consolidate --
    def shard_paths(self) -> List[str]:
        return sorted(glob.glob(os.path.join(self.shard_dir, "slot*.tsv.gz")))

    def consolidate(self, out_name: str = "cells.tsv.gz") -> str:
        """Merge every shard into one long table."""
        paths = self.shard_paths()
        if not paths:
            raise FileNotFoundError(f"no shards under {self.shard_dir}")
        frames = [pd.read_csv(p, sep="\t") for p in paths]
        df = pd.concat(frames, ignore_index=True)
        ordered = [c for c in COLUMN_ORDER if c in df.columns]
        df = df[ordered + [c for c in df.columns if c not in ordered]]
        out = os.path.join(self.run_dir, out_name)
        df.to_csv(out, sep="\t", index=False, compression="gzip")
        return out

    def load(self, name: str = "cells.tsv.gz") -> pd.DataFrame:
        """Load the consolidated table, falling back to the shards."""
        path = os.path.join(self.run_dir, name)
        if os.path.exists(path):
            return pd.read_csv(path, sep="\t")
        paths = self.shard_paths()
        if not paths:
            raise FileNotFoundError(f"no results in {self.run_dir}")
        return pd.concat([pd.read_csv(p, sep="\t") for p in paths], ignore_index=True)
