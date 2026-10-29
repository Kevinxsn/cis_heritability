"""Campaign configuration: typed, YAML-backed, self-validating.

One ``CampaignConfig`` fully determines a run.  It is serialised into the
results directory so every output table can be traced back to its inputs.
"""

from __future__ import annotations

import copy
import dataclasses
import hashlib
import json
import os
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence, Tuple

import yaml

# --------------------------------------------------------------------------
# axis vocabularies
# --------------------------------------------------------------------------
LD_ARCHITECTURES = ("regular", "high", "low")
MAF_BINS = {"all": (0.0, 0.5), "common": (0.05, 0.5), "rare": (0.01, 0.05)}
EFFECT_DISTS = ("normal", "fixed", "pareto")
LD_FOR_SUSIE = ("panel", "sim")


@dataclass
class DataConfig:
    """Where the reference genotypes live and how a cis-region is cut out."""

    geno_prefix: str = "/expanse/lustre/projects/ddp412/sux002/cis_heritability_data/1000G.EUR."
    chromosomes: List[int] = field(default_factory=lambda: list(range(1, 23)))
    window_bp: int = 500_000
    min_snps: int = 50
    max_snps: int = 25_000


@dataclass
class RegionConfig:
    """Rejection-sampling rules that define the LD architecture of a region."""

    ld_arch: str = "regular"
    r2_threshold: float = 0.8
    # `high`: keep region if frac(pairs with r2>=thr) >= frac_min
    frac_min: float = 0.02
    # `low` : keep region if frac(pairs with r2>=thr) <= frac_max
    frac_max: float = 0.001
    score_max_snps: int = 500
    max_draws: int = 500


@dataclass
class GridConfig:
    """The fully-crossed architecture grid evaluated on every region."""

    n_individuals: List[int] = field(default_factory=lambda: [300])
    ncausal: List[int] = field(default_factory=lambda: [1, 3, 5, 10, 20])
    h2: List[float] = field(default_factory=lambda: [0.05, 0.1, 0.2, 0.4, 0.6, 0.8])
    effect_dist: List[str] = field(default_factory=lambda: ["normal"])
    maf_bin: List[str] = field(default_factory=lambda: ["all"])
    # Independent phenotype draws per (region, architecture).  The default of 1
    # is the historical behaviour; >1 is what makes a WITHIN-region sampling SD
    # estimable, as opposed to dispersion across regions of differing LD.
    replicates: int = 1

    def cells(self) -> List[Dict[str, Any]]:
        """Expand the grid.  Order is stable so shards are comparable."""
        out = []
        for n in self.n_individuals:
            for maf_bin in self.maf_bin:
                for dist in self.effect_dist:
                    for ncausal in self.ncausal:
                        for h2 in self.h2:
                            for rep in range(1, int(self.replicates) + 1):
                                out.append(
                                    dict(
                                        n=int(n),
                                        maf_bin=maf_bin,
                                        effect_dist=dist,
                                        ncausal=int(ncausal),
                                        h2_true=float(h2),
                                        rep=rep,
                                    )
                                )
        return out

    @property
    def size(self) -> int:
        return (
            len(self.n_individuals)
            * len(self.maf_bin)
            * len(self.effect_dist)
            * len(self.ncausal)
            * len(self.h2)
            * int(self.replicates)
        )


@dataclass
class MethodConfig:
    """Estimator / simulation semantics.

    ``ld_for_susie``
        ``panel`` - hand SuSiE the raw reference-panel correlation matrix.
        This is what the legacy pipeline did (PLINK ``--r inter-chr``) and is
        the setting to use when reproducing historical results.
        ``sim``   - hand SuSiE the *same* regularised LD the genotypes were
        sampled from.  This is the design principle stated in the proposal
        ("the LD used to simulate must be the LD given to the fine-mapper")
        and removes LD mismatch as a confounder.

    ``ld_for_susie`` also governs the LD given to HEELS.

    ``legacy_ambiguous_filter``
        The legacy ``run_susie.R`` dropped strand-ambiguous SNPs (A/T, C/G)
        before fine-mapping even though sumstats and reference come from the
        same .bim, so there is nothing to disambiguate.  Kept as a switch for
        bit-comparable reproduction; ``false`` is the scientifically correct
        setting.
    """

    ld_regularization: float = 0.1
    ld_for_susie: str = "panel"
    legacy_ambiguous_filter: bool = True
    rescale_beta: bool = True
    causal_replace: bool = False
    susie_L: int = 10
    susie_max_iter: int = 100
    susie_coverage: float = 0.95
    max_retries: int = 5
    run_gcta: bool = True
    run_susie: bool = True
    run_heels: bool = True
    gcta_reml_no_constrain: bool = True
    # HEELS shares `ld_for_susie`: both summary-statistic estimators are handed
    # the same LD matrix, so any difference between them is the method, not the
    # input.
    heels_tol: float = 1e-5
    heels_max_iter: int = 200
    heels_init_h2: float = 0.25


@dataclass
class ComputeConfig:
    workers: int = 32
    threads_per_worker: int = 1
    scratch_root: Optional[str] = None  # defaults to <results>/scratch
    plink: str = "plink"
    gcta: str = "gcta64"
    rscript: str = "Rscript"
    keep_scratch: bool = False
    flush_every: int = 1  # regions per shard flush
    task_timeout_s: int = 7200


@dataclass
class CampaignConfig:
    name: str = "campaign"
    seed: int = 20260802
    n_regions: int = 100
    results_root: str = "/home/sux002/work/cis_heritability/results"
    description: str = ""
    data: DataConfig = field(default_factory=DataConfig)
    region: RegionConfig = field(default_factory=RegionConfig)
    grid: GridConfig = field(default_factory=GridConfig)
    method: MethodConfig = field(default_factory=MethodConfig)
    compute: ComputeConfig = field(default_factory=ComputeConfig)

    # ---------------------------------------------------------------- io ---
    @classmethod
    def from_yaml(cls, path: str, overrides: Optional[Sequence[str]] = None) -> "CampaignConfig":
        with open(path) as fh:
            raw = yaml.safe_load(fh) or {}
        for ov in overrides or []:
            _apply_override(raw, ov)
        return cls.from_dict(raw)

    @classmethod
    def from_dict(cls, raw: Dict[str, Any]) -> "CampaignConfig":
        raw = copy.deepcopy(raw)
        sub = {
            "data": DataConfig,
            "region": RegionConfig,
            "grid": GridConfig,
            "method": MethodConfig,
            "compute": ComputeConfig,
        }
        kwargs: Dict[str, Any] = {}
        for key, klass in sub.items():
            kwargs[key] = klass(**(raw.pop(key, {}) or {}))
        top = {f.name for f in dataclasses.fields(cls)}
        unknown = set(raw) - top
        if unknown:
            raise ValueError(f"unknown config keys: {sorted(unknown)}")
        kwargs.update(raw)
        cfg = cls(**kwargs)
        cfg.validate()
        return cfg

    def to_dict(self) -> Dict[str, Any]:
        return dataclasses.asdict(self)

    def save(self, path: str) -> None:
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w") as fh:
            yaml.safe_dump(self.to_dict(), fh, sort_keys=False)

    # ---------------------------------------------------------- validate ---
    def validate(self) -> None:
        if self.region.ld_arch not in LD_ARCHITECTURES:
            raise ValueError(f"region.ld_arch must be one of {LD_ARCHITECTURES}")
        for b in self.grid.maf_bin:
            if b not in MAF_BINS:
                raise ValueError(f"grid.maf_bin entries must be in {sorted(MAF_BINS)}")
        for d in self.grid.effect_dist:
            if d not in EFFECT_DISTS:
                raise ValueError(f"grid.effect_dist entries must be in {EFFECT_DISTS}")
        if self.method.ld_for_susie not in LD_FOR_SUSIE:
            raise ValueError(f"method.ld_for_susie must be one of {LD_FOR_SUSIE}")
        if not (0.0 <= self.method.ld_regularization < 1.0):
            raise ValueError("method.ld_regularization must be in [0, 1)")
        if self.n_regions < 1:
            raise ValueError("n_regions must be >= 1")
        if self.grid.replicates < 1:
            raise ValueError("grid.replicates must be >= 1")
        if self.compute.workers < 1:
            raise ValueError("compute.workers must be >= 1")

    # ------------------------------------------------------------ derived --
    @property
    def run_dir(self) -> str:
        return os.path.join(self.results_root, self.name)

    @property
    def fingerprint(self) -> str:
        """Stable hash of the scientific content of the config (not compute)."""
        payload = self.to_dict()
        payload.pop("compute", None)
        payload.pop("results_root", None)
        payload.pop("description", None)
        blob = json.dumps(payload, sort_keys=True).encode()
        return hashlib.sha1(blob).hexdigest()[:12]

    @property
    def total_cells(self) -> int:
        return self.n_regions * self.grid.size


def _apply_override(raw: Dict[str, Any], override: str) -> None:
    """Apply a ``dotted.key=value`` CLI override onto a raw config dict."""
    if "=" not in override:
        raise ValueError(f"override must look like key=value, got {override!r}")
    key, value = override.split("=", 1)
    node = raw
    parts = key.split(".")
    for part in parts[:-1]:
        node = node.setdefault(part, {})
    try:
        node[parts[-1]] = yaml.safe_load(value)
    except yaml.YAMLError:
        node[parts[-1]] = value


def maf_range(bin_name: str) -> Tuple[float, float]:
    return MAF_BINS[bin_name]
