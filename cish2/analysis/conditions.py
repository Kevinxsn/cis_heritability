"""Split campaigns into labelled *conditions* — one per level of the axis varied.

A campaign like `sample_size` sweeps five values of n inside one directory.
Reporting it as a single row or panel averages the sweep away, so the reader
cannot tell whether "sample size" meant more or less of it.  A **condition** is
one level of one axis, carrying an explicit direction relative to the baseline
(n = 300, normal effect prior, unrestricted MAF, unselected LD).

Note on the effect-size axis: `sim_beta` rescales every draw so that
``beta' R beta`` equals the target h² exactly, so all three priors carry the
*same* total genetic effect.  What they change is how concentrated that effect
is across the causal SNPs, which is the direction labelled here.  Effect
*magnitude* is already two other axes: total magnitude is ``h2_true``, and
per-causal-SNP magnitude is sqrt(h²/k), i.e. ``ncausal``.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import pandas as pd

from .load import load_run, to_long, varying_axes

# the reference cell every direction is stated against
BASELINE = {"n": 300, "effect_dist": "normal", "maf_bin": "all", "ld_arch": "regular"}

GROUP_ORDER = ["baseline", "sample size", "effect concentration", "causal MAF",
               "LD architecture"]

# level -> (short label, direction relative to baseline), per splittable axis
_EFFECT_DIST = {
    "fixed": ("equal-magnitude effects", "most spread out"),
    "normal": ("normal effects", "baseline"),
    "pareto": ("heavy-tailed effects", "most concentrated"),
}
_EFFECT_DIST_ORDER = ["fixed", "normal", "pareto"]

_MAF_BIN = {
    "common": ("common causal SNPs", "MAF 0.05-0.5"),
    "rare": ("rare causal SNPs", "MAF 0.01-0.05"),
    "all": ("any causal SNP", "baseline"),
}
_MAF_BIN_ORDER = ["common", "rare", "all"]

_LD_ARCH = {
    "low": ("low-LD regions", "bottom decile"),
    "regular": ("unselected regions", "baseline"),
    "high": ("high-LD regions", "top decile"),
}

# which axis a campaign is split on, if any
SPLIT_AXES = ["n", "effect_dist", "maf_bin"]


@dataclass
class Condition:
    """One level of one design axis, with its data and its direction."""

    key: str
    label: str
    direction: str
    group: str
    campaign: str
    wide: pd.DataFrame = field(repr=False)
    long: pd.DataFrame = field(repr=False)

    @property
    def is_baseline(self) -> bool:
        return self.group == "baseline"

    @property
    def title(self) -> str:
        return f"{self.label}\n{self.direction}"


def _describe(axis: str, level) -> Tuple[str, str, str]:
    """Return ``(label, direction, group)`` for one level of one axis."""
    if axis == "n":
        base = BASELINE["n"]
        if level == base:
            direction = "baseline"
        else:
            direction = f"{'above' if level > base else 'below'} baseline (n = {base})"
        return f"n = {level:,}", direction, "sample size"
    if axis == "effect_dist":
        label, direction = _EFFECT_DIST.get(str(level), (str(level), ""))
        return label, direction, "effect concentration"
    if axis == "maf_bin":
        label, direction = _MAF_BIN.get(str(level), (str(level), ""))
        return label, direction, "causal MAF"
    if axis == "ld_arch":
        label, direction = _LD_ARCH.get(str(level), (str(level), ""))
        return label, direction, "LD architecture"
    return f"{axis} = {level}", "", axis


def _levels_in_order(axis: str, values) -> List:
    vals = list(values)
    if axis == "effect_dist":
        return [v for v in _EFFECT_DIST_ORDER if v in vals] + [v for v in vals if v not in _EFFECT_DIST_ORDER]
    if axis == "maf_bin":
        return [v for v in _MAF_BIN_ORDER if v in vals] + [v for v in vals if v not in _MAF_BIN_ORDER]
    return sorted(vals)


def conditions_for_run(run_dir: str) -> List[Condition]:
    """Split one campaign into conditions, one per level of the axis it varies."""
    campaign = os.path.basename(os.path.abspath(run_dir))
    wide = load_run(run_dir)
    extra = [a for a in varying_axes(wide) if a in SPLIT_AXES]

    if not extra:
        # a whole-campaign condition: its identity is the LD architecture it fixed
        arch = str(wide["ld_arch"].iloc[0]) if "ld_arch" in wide.columns else "regular"
        label, direction, group = _describe("ld_arch", arch)
        if arch == BASELINE["ld_arch"]:
            group = "baseline"
            label, direction = "baseline", f"n = {BASELINE['n']}, normal effects"
        return [Condition(key=campaign, label=label, direction=direction, group=group,
                          campaign=campaign, wide=wide, long=to_long(wide))]

    if len(extra) > 1:  # not produced by the shipped configs, but keep it honest
        extra = extra[:1]
    axis = extra[0]

    out: List[Condition] = []
    for level in _levels_in_order(axis, wide[axis].dropna().unique()):
        part = wide[wide[axis] == level].reset_index(drop=True)
        if part.empty:
            continue
        label, direction, group = _describe(axis, level)
        out.append(Condition(key=f"{campaign}:{axis}={level}", label=label,
                             direction=direction, group=group, campaign=campaign,
                             wide=part, long=to_long(part)))
    return out


def collect(run_dirs: List[str]) -> List[Condition]:
    """Build every condition across campaigns, ordered baseline-first by group."""
    found: List[Condition] = []
    for run in run_dirs:
        found.extend(conditions_for_run(run))

    def sort_key(c: Condition) -> Tuple[int, int, str]:
        group_rank = GROUP_ORDER.index(c.group) if c.group in GROUP_ORDER else len(GROUP_ORDER)
        # inside a group, keep the order conditions_for_run produced
        return (group_rank, found.index(c), c.key)

    return sorted(found, key=sort_key)


def baseline_condition(conditions: List[Condition]) -> Optional[Condition]:
    for c in conditions:
        if c.is_baseline:
            return c
    return None
