"""Analyses added in response to PI review of the IGVF 2026 poster.

Three things the review asked for that the existing tables could not answer:

  A. Agreement, not just correlation.  Pearson r is insensitive to offset and
     scale, so "HEELS and GCTA are one method" was not supported by r alone.
     Computes mean difference, limits of agreement, Deming-style calibration
     slope/intercept and RMSE, per sample size.

  B. Heatmap conditioning.  The published heatmap stratified genes by GCTA h2
     and plotted CHEWS - GCTA, so GCTA sat on both axes and regression to the
     mean can manufacture the pattern.  This re-runs the same stratification
     using a DISJOINT half of the donors to define the strata, which removes
     the shared-noise coupling entirely.

  C. Split-half selection.  The reliability set conditioned on h2 > 0.05 in
     all six estimates.  Reports the same coefficients on prespecified gene
     sets chosen independently of the split-half estimates.

Writes tables/ and prints a summary.  Read by make_figures.py.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent / "results"
OUT = HERE / "tables"
OUT.mkdir(exist_ok=True)


# ---------------------------------------------------------------------------
# A. Agreement statistics
# ---------------------------------------------------------------------------
def deming(x, y, lam=1.0):
    """Errors-in-variables slope. Both estimates are noisy, so ordinary least
    squares would attenuate the slope toward zero and understate agreement."""
    mx, my = x.mean(), y.mean()
    sxx = ((x - mx) ** 2).mean()
    syy = ((y - my) ** 2).mean()
    sxy = ((x - mx) * (y - my)).mean()
    d = syy - lam * sxx
    slope = (d + np.sqrt(d ** 2 + 4 * lam * sxy ** 2)) / (2 * sxy)
    return slope, my - slope * mx


def agreement(a, b):
    """Bland-Altman style summary of b (test) against a (reference)."""
    ok = np.isfinite(a) & np.isfinite(b)
    a, b = np.asarray(a)[ok], np.asarray(b)[ok]
    d = b - a
    slope, icpt = deming(a, b)
    return dict(
        n=len(a),
        mean_diff=d.mean(),
        sd_diff=d.std(ddof=1),
        loa_lo=d.mean() - 1.96 * d.std(ddof=1),
        loa_hi=d.mean() + 1.96 * d.std(ddof=1),
        rmse=float(np.sqrt((d ** 2).mean())),
        pearson=float(np.corrcoef(a, b)[0, 1]),
        slope=slope,
        intercept=icpt,
    )


def part_a():
    g = pd.read_csv(RESULTS / "gtex_whole_blood" / "cells.tsv.gz", sep="\t")
    rows = []
    for n, part in g.groupby("n"):
        for name, col in (("CHEWS", "susie_h2"), ("HEELS", "heels_h2")):
            r = agreement(part["gcta_h2"].to_numpy(), part[col].to_numpy())
            r.update(n=n, estimator=name, genes=r.pop("n"))
            rows.append(r)
    t = pd.DataFrame(rows)[["n", "estimator", "genes", "mean_diff", "sd_diff",
                            "loa_lo", "loa_hi", "rmse", "pearson", "slope", "intercept"]]
    t.to_csv(OUT / "gtex_agreement.csv", index=False)
    print("\nA. Agreement vs GCTA (GTEx Whole Blood)")
    print(t.round(4).to_string(index=False))
    return t


# ---------------------------------------------------------------------------
# B. Heatmap stratified on a DISJOINT donor half
# ---------------------------------------------------------------------------
H2_BINS = [0.0, 0.05, 0.1, 0.2, 0.4, 1.0]
H2_LABELS = ["0–.05", ".05–.1", ".1–.2", ".2–.4", ".4–1"]


def clump_bin(v):
    return "4+" if v >= 4 else str(int(v))


def part_b():
    s = pd.read_csv(RESULTS / "gtex_split_half" / "cells.tsv.gz", sep="\t")
    a = s[s["subset"] == "split1"].set_index("gene_id")
    b = s[s["subset"] == "split2"].set_index("gene_id")
    common = a.index.intersection(b.index)
    a, b = a.loc[common], b.loc[common]

    # strata come from half A; the gap is measured in half B.  No donor is
    # shared, so the stratification variable carries none of half B's noise.
    d = pd.DataFrame({
        "strat_h2": a["gcta_h2"].to_numpy(),
        "strat_clumps": a["n_clumps"].to_numpy(),
        "gap": (b["susie_h2"] - b["gcta_h2"]).to_numpy(),
    }, index=common).dropna()
    d["h2_bin"] = pd.cut(d["strat_h2"], H2_BINS, labels=H2_LABELS)
    d["clump_bin"] = d["strat_clumps"].map(clump_bin)

    order = ["0", "1", "2", "3", "4+"]
    mean = d.pivot_table(index="clump_bin", columns="h2_bin", values="gap",
                         aggfunc="mean", observed=False).reindex(order)
    cnt = d.pivot_table(index="clump_bin", columns="h2_bin", values="gap",
                        aggfunc="size", observed=False).reindex(order)
    sem = d.pivot_table(index="clump_bin", columns="h2_bin", values="gap",
                        aggfunc="sem", observed=False).reindex(order)
    mean = mean.where(cnt >= 30)
    sem = sem.where(cnt >= 30)

    mean.to_csv(OUT / "gap_disjoint_strata.csv")
    cnt.to_csv(OUT / "gap_disjoint_counts.csv")
    sem.to_csv(OUT / "gap_disjoint_sem.csv")

    print("\nB. CHEWS - GCTA in half B, stratified on half A (disjoint donors)")
    print(mean.round(4).to_string())
    print("\n   genes per cell")
    print(cnt.fillna(0).astype(int).to_string())

    # the same thing done the coupled way, for direct comparison
    dc = pd.DataFrame({
        "strat_h2": b["gcta_h2"].to_numpy(),
        "strat_clumps": b["n_clumps"].to_numpy(),
        "gap": (b["susie_h2"] - b["gcta_h2"]).to_numpy(),
    }, index=common).dropna()
    dc["h2_bin"] = pd.cut(dc["strat_h2"], H2_BINS, labels=H2_LABELS)
    dc["clump_bin"] = dc["strat_clumps"].map(clump_bin)
    coupled = dc.pivot_table(index="clump_bin", columns="h2_bin", values="gap",
                             aggfunc="mean", observed=False).reindex(order)
    coupled = coupled.where(cnt >= 30)
    coupled.to_csv(OUT / "gap_coupled_strata.csv")
    print("\n   same-half (coupled) version, for comparison")
    print(coupled.round(4).to_string())
    print(f"\n   mean |coupled - disjoint| = "
          f"{(coupled - mean).abs().stack().mean():.4f}")
    return mean, cnt, sem, coupled


# ---------------------------------------------------------------------------
# C. Split-half reliability without conditioning on the split-half estimates
# ---------------------------------------------------------------------------
def icc_a1(a, b):
    """Two-way ICC(A,1) - absolute agreement, so a constant offset is penalised."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    n = len(a)
    m = (a + b) / 2
    msr = 2 * ((m - m.mean()) ** 2).sum() / (n - 1)
    mse = ((a - b) ** 2).sum() / (2 * n)
    msc = n * ((a.mean() - (a.mean() + b.mean()) / 2) ** 2
               + (b.mean() - (a.mean() + b.mean()) / 2) ** 2)
    return (msr - mse) / (msr + mse + 2 * (msc - mse) / n)


def part_c():
    s = pd.read_csv(RESULTS / "gtex_split_half" / "cells.tsv.gz", sep="\t")
    a = s[s["subset"] == "split1"].set_index("gene_id")
    b = s[s["subset"] == "split2"].set_index("gene_id")
    common = a.index.intersection(b.index)
    a, b = a.loc[common], b.loc[common]

    full = pd.read_csv(RESULTS / "gtex_whole_blood" / "cells.tsv.gz", sep="\t")
    full = full[full["n"] == 574].set_index("gene_id")
    # a gene set fixed BEFORE the split: has a genome-wide cis-eQTL in the full
    # sample, which is a statement about the data, not about either half's h2.
    has_eqtl = set(full.index[full["n_clumps"] >= 1])

    sets = {
        "all genes (no filter)": common,
        "prespecified: >=1 cis-eQTL in the full 574": [g for g in common if g in has_eqtl],
        "post-hoc: h2>0.05 in all six estimates":
            common[((a[["gcta_h2", "susie_h2", "heels_h2"]] > 0.05).all(axis=1)
                    & (b[["gcta_h2", "susie_h2", "heels_h2"]] > 0.05).all(axis=1)).to_numpy()],
    }

    rows = []
    for label, idx in sets.items():
        idx = pd.Index(idx)
        for name, col in (("GCTA", "gcta_h2"), ("CHEWS", "susie_h2"), ("HEELS", "heels_h2")):
            x, y = a.loc[idx, col].to_numpy(), b.loc[idx, col].to_numpy()
            ok = np.isfinite(x) & np.isfinite(y)
            rows.append(dict(gene_set=label, estimator=name, genes=int(ok.sum()),
                             pearson=float(np.corrcoef(x[ok], y[ok])[0, 1]),
                             icc=icc_a1(x[ok], y[ok]),
                             rmsd=float(np.sqrt(((x[ok] - y[ok]) ** 2).mean()))))
    t = pd.DataFrame(rows)
    t.to_csv(OUT / "reliability_by_gene_set.csv", index=False)
    print("\nC. Split-half reliability under three gene-selection rules")
    print(t.round(4).to_string(index=False))

    # Is CHEWS's reliability edge over GCTA larger than its own uncertainty?
    # Paired bootstrap over genes on the prespecified set.
    idx = pd.Index(sets["prespecified: >=1 cis-eQTL in the full 574"])
    xg, yg = a.loc[idx, "gcta_h2"].to_numpy(), b.loc[idx, "gcta_h2"].to_numpy()
    xs, ys = a.loc[idx, "susie_h2"].to_numpy(), b.loc[idx, "susie_h2"].to_numpy()
    ok = np.isfinite(xg) & np.isfinite(yg) & np.isfinite(xs) & np.isfinite(ys)
    xg, yg, xs, ys = xg[ok], yg[ok], xs[ok], ys[ok]
    rng = np.random.default_rng(20260908)
    diffs = []
    for _ in range(2000):
        j = rng.integers(0, len(xg), len(xg))
        diffs.append(np.corrcoef(xs[j], ys[j])[0, 1] - np.corrcoef(xg[j], yg[j])[0, 1])
    diffs = np.array(diffs)
    lo, hi = np.percentile(diffs, [2.5, 97.5])
    print(f"\n   CHEWS - GCTA split-half r, prespecified set: "
          f"{np.corrcoef(xs, ys)[0,1] - np.corrcoef(xg, yg)[0,1]:+.4f} "
          f"(95% bootstrap CI {lo:+.4f} to {hi:+.4f}, {len(xg)} genes)")
    return t, (lo, hi)


# ---------------------------------------------------------------------------
# D. Gene inclusion flow
# ---------------------------------------------------------------------------
def part_d():
    full = pd.read_csv(RESULTS / "gtex_whole_blood" / "cells.tsv.gz", sep="\t")
    print("\nD. Gene inclusion flow (GTEx Whole Blood)")
    print(f"   genes with expression, genotypes and covariates : "
          f"{full['gene_id'].nunique():,}")
    for n, part in full.groupby("n"):
        print(f"   n = {n:>3}:  attempted {len(part):>6,}   "
              f"GCTA ok {part['gcta_h2'].notna().sum():>6,}   "
              f"CHEWS ok {part['susie_h2'].notna().sum():>6,}   "
              f"HEELS ok {part['heels_h2'].notna().sum():>6,}   "
              f"all three {int((part[['gcta_h2','susie_h2','heels_h2']].notna().all(axis=1)).sum()):>6,}")


# ===========================================================================
# Round 4 (PI feedback) - docs/pi_feedback_round4_design.md
# ===========================================================================
def _poster_helpers():
    """pairwise_win / region_bootstrap_ci live in make_figures; reuse them so
    the win rate here is by construction the one the figures draw."""
    import make_figures as mf
    return mf


def part_e():
    """§5 panel B: how OFTEN CHEWS is closer to truth, by (causal SNPs, n).

    Pairwise CHEWS-vs-GCTA (never the pipeline's three-way susie_win_rate,
    trap 15).  Pooled over the six true h2 values with equal weight per fit,
    which is how the round-3 MAE panel pooled.  CIs resample whole regions
    (2,000 draws): a region supplies every h2 and causal count, so a binomial
    SE would be too small.

    rel_mae = MAE_CHEWS / MAE_GCTA - 1 is the magnitude companion: scale-free
    because it divides by the competitor's error, not by h2 (per-fit "% of
    truth" error reweights toward h2 = 0.05 and can flip a cell's sign).
    """
    mf = _poster_helpers()
    cells = pd.read_csv(RESULTS / "sample_size" / "cells.tsv.gz", sep="\t")
    d = cells.dropna(subset=["susie_h2", "gcta_h2"]).reset_index(drop=True)
    d["_ec"] = (d["susie_h2"] - d["h2_true"]).abs()
    d["_eg"] = (d["gcta_h2"] - d["h2_true"]).abs()
    d["_win"] = (d["_ec"] < d["_eg"]).astype(float)
    d["_cell"] = d["ncausal"].astype(int) * 10000 + d["n"].astype(int)

    # pairwise_win names its count column "n", so the sample size goes in as n_ind
    win = mf.pairwise_win(d.rename(columns={"n": "n_ind"}), ["ncausal", "n_ind"])
    win = win.drop(columns=["n"]).rename(columns={"n_ind": "n"})
    keys, lo, hi = mf.region_bootstrap_ci(d, lambda b: b["_win"], "_cell")

    # ratio of means needs its own resampling loop - same seed and the same
    # region draws as region_bootstrap_ci, so the two CIs are paired
    rng = np.random.default_rng(20260908)
    regions = d["region_id"].unique()
    by_region = {r: g.index.to_numpy() for r, g in d.groupby("region_id")}
    ratios = {k: [] for k in keys}
    for _ in range(2000):
        pick = rng.choice(regions, size=len(regions), replace=True)
        b = d.loc[np.concatenate([by_region[r] for r in pick])]
        m = b.groupby("_cell")[["_ec", "_eg"]].mean()
        r = m["_ec"] / m["_eg"] - 1
        for k in keys:
            ratios[k].append(r.get(k, np.nan))

    rows = []
    for _, w in win.iterrows():
        k = int(w["ncausal"]) * 10000 + int(w["n"])
        part = d[d["_cell"] == k]
        mae_c, mae_g = part["_ec"].mean(), part["_eg"].mean()
        rows.append(dict(
            ncausal=int(w["ncausal"]), n=int(w["n"]),
            rate=float(w["rate"]), lo=lo[k], hi=hi[k],
            mae_chews=mae_c, mae_gcta=mae_g, rel_mae=mae_c / mae_g - 1,
            rel_mae_lo=float(np.nanpercentile(ratios[k], 2.5)),
            rel_mae_hi=float(np.nanpercentile(ratios[k], 97.5)),
            n_pairs=int(len(part)),
        ))
    t = pd.DataFrame(rows).sort_values(["ncausal", "n"])
    t.to_csv(OUT / "win_by_n_causal.csv", index=False)
    print("\nE. Pairwise win rate (CHEWS closer to truth) by causal SNPs x n")
    print((100 * t.pivot(index="ncausal", columns="n", values="rate")).round(0).to_string())
    print("\n   relative MAE change, MAE_C / MAE_G - 1 (%)")
    print((100 * t.pivot(index="ncausal", columns="n", values="rel_mae")).round(0).to_string())
    print("\n   cells whose 95% CI excludes 50%:",
          int(((t["lo"] > 0.5) | (t["hi"] < 0.5)).sum()), "of", len(t))
    print("   paired fits per cell:", t["n_pairs"].min(), "-", t["n_pairs"].max())
    return t


# ---------------------------------------------------------------------------
# helpers shared by F and G: half A defines the stratum, half B is measured
# ---------------------------------------------------------------------------
GCTA_FLOOR = 1e-4     # constrained REML pinned at its 0 bound reports ~1e-6
CHEWS_FLOOR = 0.005   # "CHEWS sets the gene to ~0" (design §3.1 table)
ORDER5 = ["0", "1", "2", "3", "4+"]


def _halves(s, key=("gene_id",)):
    """Pair split1 (half A) with split2 (half B) on ``key``."""
    key = list(key)
    a = s[s["subset"] == "split1"].set_index(key)
    b = s[s["subset"] == "split2"].set_index(key)
    common = a.index.intersection(b.index)
    return a.loc[common], b.loc[common]


def _strata(a):
    """Exactly part_b's strata: half-A constrained GCTA h2 x half-A clump count."""
    h2_bin = pd.cut(a["gcta_h2"], H2_BINS, labels=H2_LABELS)
    clump = a["n_clumps"].map(clump_bin)
    return h2_bin.astype(object), clump


def part_f():
    """B2. Permutation null on real expression, and a Storey-type pi0.

    Genotype rows are permuted within each half (cish2 phenotype_mode=permuted),
    so every gene's true cis-h2 is 0 while the expression and covariates are the
    real ones.  Reports where each estimator sits at h2 = 0 - overall, and in
    half B within the permuted half-A strata - and, per REAL stratum,

        pi0_hat = P_real(p in E) / P_perm(p in E)

    on GCTA's LRT p-value in half B, with the permutation null taken on the SAME
    genes (constrained-REML p-values have a point mass at the bound, so the null
    is empirical, not Uniform(0,1)).  Capped at 1.

    Deviation from the design (which said E = {p > 0.5}): GCTA reports p = 0.5
    exactly when REML sits at the 0 bound and never above it, so P(p > 0.5) = 0
    in both runs and that ratio is 0/0.  Primary E = {p = 0.5} ("at the bound",
    the most null-specific event there is); E = {p > 0.25} and {p > 0.1} are
    reported as sensitivity.
    """
    perm_path = RESULTS / "gtex_perm_split_half" / "cells.tsv.gz"
    if not perm_path.exists():
        print("\nF. skipped: results/gtex_perm_split_half has not been run")
        return None
    perm = pd.read_csv(perm_path, sep="\t")
    real = pd.read_csv(RESULTS / "gtex_split_half" / "cells.tsv.gz", sep="\t")

    rows = []
    for est, col, floor in (("GCTA", "gcta_h2", GCTA_FLOOR), ("CHEWS", "susie_h2", CHEWS_FLOOR),
                            ("HEELS", "heels_h2", CHEWS_FLOOR)):
        for half, part in perm.groupby("subset"):
            v = part[col].dropna()
            rows.append(dict(scope=f"all genes, {half}", estimator=est, genes=len(v),
                             mean=v.mean(), median=v.median(), p90=v.quantile(0.9),
                             pct_at_floor=100 * (v < floor).mean()))
    pa, pb = _halves(perm)
    h2b, cb = _strata(pa)
    d = pd.DataFrame({"h2_bin": h2b, "clump_bin": cb, "gcta": pb["gcta_h2"],
                      "chews": pb["susie_h2"]}).dropna()
    for (hb, clb), g in d.groupby(["h2_bin", "clump_bin"]):
        if hb == H2_LABELS[0] and clb in ("0", "1"):
            for est, col, floor in (("GCTA", "gcta", GCTA_FLOOR), ("CHEWS", "chews", CHEWS_FLOOR)):
                rows.append(dict(scope=f"half B | permuted half-A stratum {clb} signals, {hb}",
                                 estimator=est, genes=len(g), mean=g[col].mean(),
                                 median=g[col].median(), p90=g[col].quantile(0.9),
                                 pct_at_floor=100 * (g[col] < floor).mean()))
    null = pd.DataFrame(rows)
    null.to_csv(OUT / "perm_null.csv", index=False)
    print("\nF. Permutation null (true h2 = 0 on real expression)")
    print(null.round(4).to_string(index=False))

    # pi0 per real stratum: real strata from half A, p-values from half B
    ra, rb = _halves(real)
    h2r, cr = _strata(ra)
    pperm = pb["gcta_pval"]
    t = pd.DataFrame({"h2_bin": h2r, "clump_bin": cr, "p_real": rb["gcta_pval"]})
    t["p_perm"] = pperm.reindex(t.index)
    t = t.dropna()
    events = {"pi0": lambda p: p >= 0.5 - 1e-9,          # at the bound (primary)
              "pi0_lam25": lambda p: p > 0.25, "pi0_lam10": lambda p: p > 0.10}

    def ratio(g):
        r = dict(genes=len(g), at_bound_real=(g["p_real"] >= 0.5 - 1e-9).mean(),
                 at_bound_perm=(g["p_perm"] >= 0.5 - 1e-9).mean())
        for name, ev in events.items():
            den = ev(g["p_perm"]).mean()
            r[name] = float(np.minimum(1.0, ev(g["p_real"]).mean() / den)) if den > 0 else np.nan
        return r

    out = [dict(h2_bin=hb, clump_bin=clb, **ratio(g))
           for (hb, clb), g in t.groupby(["h2_bin", "clump_bin"])]
    out.append(dict(h2_bin="all", clump_bin="all", **ratio(t)))
    pi0 = pd.DataFrame(out)
    pi0.to_csv(OUT / "pi0_by_stratum.csv", index=False)
    print("\n   Storey-type pi0 per real stratum (E = REML at the bound, half B; "
          "real strata from half A)")
    piv = pi0[pi0.h2_bin != "all"].pivot(index="clump_bin", columns="h2_bin", values="pi0")
    print(piv.reindex(ORDER5)[H2_LABELS].round(3).to_string())
    last = pi0.iloc[-1]
    print(f"   all genes: pi0 = {last['pi0']:.3f}  (sensitivity: p>0.25 -> "
          f"{last['pi0_lam25']:.3f}, p>0.1 -> {last['pi0_lam10']:.3f})")
    return null, pi0


def part_g():
    """B1 + B3. Semi-synthetic GTEx strata, calibrated by importance reweighting.

    WRITTEN BEFORE ANY B1 OUTPUT WAS READ (2026-10-28).  The decision rules and
    every threshold below are fixed here so the interpretation cannot drift.

    Data.  results/gtex_sim_split_half: each GTEx gene gets K = 2 simulated
    phenotypes on its real genotypes (cish2 phenotype_mode=simulated), fitted on
    the same two donor halves as the real pass.  A (gene, draw) is one simulated
    gene.  Strata exactly as part_b: half-A constrained GCTA h2 (H2_BINS) x
    half-A clump count (clump_bin).  Everything is measured in half B.  Truth =
    h2_target (the population parameter); h2_resid is the robustness check.

    Win = |CHEWS_B - truth| < |GCTA_B - truth| (pairwise; constrained GCTA, the
    estimator the real panel uses).

    Reweighting.  Design prior: P(h2=0) = null_mass, the rest uniform over the
    7 x 6 non-null (h2, ncausal) grid.  A target prior on the same grid is
        pi0                       at h2 = 0
        (1 - pi0) * a(h2) * b(k)  otherwise,
        a(h2) ~ exp(-h2 / h2_scale)  (normalised over the non-null h2 grid)
        b(k)  ~ q ** rank(k)         (rank 0..5 over ncausal_grid; q < 1 = sparse)
    and each simulated gene gets w = pi_target / pi_design.  Grid:
    pi0 in 0, 0.1, ..., 0.9; h2_scale in {0.005, 0.01, 0.02, 0.05, 0.1, 0.2,
    0.5, inf (= uniform)}; q in {0.25, 0.5, 0.75, 1, 1.5, 2}.

    Score (lower is better), over the real figure's unmasked cells (real
    count >= 30; the masking rule of part_b):
        S = sum_c (gap_sim - gap_real)^2 / (sem_real^2 + sem_sim^2)
          + sum_c (n_real - N * p_sim)^2 / (N * p_sim)
    i.e. a chi-square on the half-B gap means plus a Pearson chi-square on the
    half-A stratum counts (N = real genes in unmasked cells, p_sim the weighted
    share of simulated genes in that cell among unmasked cells).  sem_sim is the
    weighted SE of the simulated gap.  Acceptable set:
        S <= S_min + 7.81 * max(1, S_min / (n_terms - 3))
    (chi-square 95% for 3 fitted knobs, inflated by the reduced chi-square when
    the best model under-fits - a Birge-ratio correction).

    Per stratum and prior: win rate, MAE_C - MAE_G, bias of each method, Kish
    ESS; win-rate 95% CI from a Poisson bootstrap over GENES (500 draws; both
    draws of a gene share its LD, so genes are the cluster).  Cells with ESS < 30
    under the best prior are masked.

    DECISION RULES (design §3.3), applied to the 0-0.05 column (0 and 1 signal
    cells are the poster sentence) and then to the other strata:

      1. Calibration fails if, at the best-fit prior, the sign of the simulated
         gap agrees with the real gap in < 80% of unmasked cells, OR the mean
         absolute deviation of simulated from real gap means exceeds 0.010, OR
         the total-variation distance between simulated and real stratum shares
         exceeds 0.15.  -> B1 does NOT go on the poster; limitations say "the
         simulation could not reproduce the real pattern"; §7 stays as in round 3.
      2. Win-rate CI > 50% for EVERY prior in the acceptable set -> "CHEWS's
         lower values there are closer to truth: GCTA's constrained floor sits
         above it."
      3. Win-rate CI < 50% for every acceptable prior -> "CHEWS's
         underestimation there is a real loss: it zeroes genes with small
         diffuse h2."
      4. Otherwise (the acceptable set straddles 50%) -> "depends on the null
         fraction, which the data do not pin down": show win rate against pi0
         and give part_f's pi0_hat.
    The hatch in P6_better_or_different marks cells under rule 4.
    """
    sim_path = RESULTS / "gtex_sim_split_half" / "cells.tsv.gz"
    if not sim_path.exists():
        print("\nG. skipped: results/gtex_sim_split_half has not been run")
        return None
    import json
    import yaml

    cfg = yaml.safe_load(open(RESULTS / "gtex_sim_split_half" / "config.yaml"))
    design = cfg["sim_design"]
    sim = pd.read_csv(sim_path, sep="\t")
    a, b = _halves(sim, key=("gene_id", "draw"))
    h2b, cb = _strata(a)
    d = pd.DataFrame({
        "h2_bin": h2b, "clump_bin": cb,
        "h2": b["h2_target"], "k": b["ncausal_target"], "h2_resid": b["h2_resid"],
        "gcta": b["gcta_h2"], "gcta_uc": b.get("gcta_h2_uc"), "chews": b["susie_h2"],
    })
    d = d.dropna(subset=["h2_bin", "gcta", "chews"]).reset_index()
    d["cell"] = d["clump_bin"] + "|" + d["h2_bin"]
    d["eC"] = (d["chews"] - d["h2"]).abs()
    d["eG"] = (d["gcta"] - d["h2"]).abs()
    d["win"] = (d["eC"] < d["eG"]).astype(float)
    d["win_resid"] = ((d["chews"] - d["h2_resid"]).abs()
                      < (d["gcta"] - d["h2_resid"]).abs()).astype(float)
    d["gap"] = d["chews"] - d["gcta"]
    print(f"\nG. simulated genes (gene x draw) with both halves: {len(d):,} "
          f"({d['gene_id'].nunique():,} genes)")

    # ---- real targets (part_b's tables) ------------------------------------
    real_gap = pd.read_csv(OUT / "gap_disjoint_strata.csv", index_col=0)
    real_cnt = pd.read_csv(OUT / "gap_disjoint_counts.csv", index_col=0)
    real_sem = pd.read_csv(OUT / "gap_disjoint_sem.csv", index_col=0)
    for t in (real_gap, real_cnt, real_sem):
        t.index = t.index.astype(str)
    cells = [f"{r}|{c}" for r in ORDER5 for c in H2_LABELS
             if real_cnt.loc[r, c] >= 30 and np.isfinite(real_gap.loc[r, c])]
    rg = np.array([real_gap.loc[c.split("|")[0], c.split("|")[1]] for c in cells])
    rs = np.array([real_sem.loc[c.split("|")[0], c.split("|")[1]] for c in cells])
    rn = np.array([real_cnt.loc[c.split("|")[0], c.split("|")[1]] for c in cells], float)
    N = rn.sum()

    # ---- priors --------------------------------------------------------------
    h2_nonnull = [float(h) for h in design["h2_grid"] if float(h) > 0]
    ks = [int(k) for k in design["ncausal_grid"]]
    null_mass = float(design["null_mass"])
    p_design = np.where(d["h2"] <= 0, null_mass,
                        (1 - null_mass) / (len(h2_nonnull) * len(ks)))
    h2_idx = d["h2"].map({h: i for i, h in enumerate(h2_nonnull)})
    k_idx = d["k"].map({k: i for i, k in enumerate(ks)})

    def target_weights(pi0, scale, q):
        ah = np.ones(len(h2_nonnull)) if np.isinf(scale) else np.exp(-np.array(h2_nonnull) / scale)
        ah /= ah.sum()
        bk = q ** np.arange(len(ks))
        bk /= bk.sum()
        pt = np.where(d["h2"] <= 0, pi0,
                      (1 - pi0) * ah[h2_idx.fillna(0).astype(int)] * bk[k_idx.fillna(0).astype(int)])
        return pt / p_design

    cell_codes = {c: i for i, c in enumerate(cells)}
    ci = d["cell"].map(cell_codes)
    inside = ci.notna().to_numpy()
    ci = ci.fillna(-1).astype(int).to_numpy()
    onehot = np.zeros((len(d), len(cells)))
    onehot[np.flatnonzero(inside), ci[inside]] = 1.0

    def cell_stats(w, col):
        sw = w @ onehot
        mean = ((w * d[col].to_numpy()) @ onehot) / np.where(sw > 0, sw, np.nan)
        return sw, mean

    grid = [(pi0, sc, q) for pi0 in np.round(np.arange(0, 1.0, 0.1), 1)
            for sc in (0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, np.inf)
            for q in (0.25, 0.5, 0.75, 1.0, 1.5, 2.0)]
    fits = []
    gap = d["gap"].to_numpy()
    for pi0, sc, q in grid:
        w = target_weights(pi0, sc, q)
        sw, gmean = cell_stats(w, "gap")
        # weighted SE of a weighted mean, with Kish ESS
        sw2 = (w ** 2) @ onehot
        ess = sw ** 2 / np.where(sw2 > 0, sw2, np.nan)
        dev = np.nan_to_num((gap - gmean[np.maximum(ci, 0)]) ** 2) * inside
        var = ((w * dev) @ onehot) / np.where(sw > 0, sw, np.nan)
        sem = np.sqrt(var / ess)
        psim = sw / sw.sum()
        chi_gap = np.nansum((gmean - rg) ** 2 / (rs ** 2 + sem ** 2))
        chi_cnt = np.nansum((rn - N * psim) ** 2 / np.where(psim > 0, N * psim, np.nan))
        fits.append(dict(pi0=pi0, h2_scale=sc, q=q, S=chi_gap + chi_cnt, chi_gap=chi_gap,
                         chi_cnt=chi_cnt,
                         sign_agree=float(np.mean(np.sign(gmean) == np.sign(rg))),
                         mad_gap=float(np.nanmean(np.abs(gmean - rg))),
                         tvd=float(0.5 * np.abs(psim - rn / N).sum())))
    fits = pd.DataFrame(fits).sort_values("S").reset_index(drop=True)
    n_terms = 2 * len(cells)
    s_min = fits["S"].iloc[0]
    tol = 7.81 * max(1.0, s_min / (n_terms - 3))
    fits["acceptable"] = fits["S"] <= s_min + tol
    fits.to_csv(OUT / "calibration_fit.csv", index=False)
    best = fits.iloc[0]
    acc = fits[fits["acceptable"]]
    print(f"   calibration: {len(grid)} priors, best S = {s_min:.1f} over {n_terms} terms, "
          f"tolerance {tol:.1f}, {len(acc)} acceptable")
    print("   best:", best[["pi0", "h2_scale", "q", "S", "chi_gap", "chi_cnt",
                          "sign_agree", "mad_gap", "tvd"]].to_dict())
    print(f"   acceptable pi0 range {acc['pi0'].min()}-{acc['pi0'].max()}, "
          f"h2_scale {sorted(acc['h2_scale'].unique())}, q {sorted(acc['q'].unique())}")
    calib_ok = (best["sign_agree"] >= 0.80 and best["mad_gap"] <= 0.010 and best["tvd"] <= 0.15)
    crit = [("sign agreement", best["sign_agree"], ">=", 0.80),
            ("mean |sim - real| gap", best["mad_gap"], "<=", 0.010),
            ("TVD of stratum shares", best["tvd"], "<=", 0.15)]
    print(f"   rule 1 (calibration): {'PASSES' if calib_ok else 'FAILS'}")
    for name, v, op, lim in crit:
        ok = v >= lim if op == ">=" else v <= lim
        print(f"     {'ok  ' if ok else 'FAIL'} {name} = {v:.4f} (needs {op} {lim})")

    # ---- per-stratum metrics: every stratum of the 5 x 5 grid ----------------
    all_cells = [f"{r}|{c}" for r in ORDER5 for c in H2_LABELS]
    full_codes = {c: i for i, c in enumerate(all_cells)}
    fi = d["cell"].map(full_codes).to_numpy()
    oh = np.zeros((len(d), len(all_cells)))
    oh[np.arange(len(d)), fi] = 1.0
    genes = d["gene_id"].astype("category").cat.codes.to_numpy()
    n_genes = genes.max() + 1
    rng = np.random.default_rng(20261028)
    boot = rng.poisson(1.0, size=(500, n_genes)).astype(float)[:, genes]   # 500 x rows

    def metrics(w):
        sw = w @ oh
        safe = np.where(sw > 0, sw, np.nan)
        out = {"ess": sw ** 2 / np.where((w ** 2) @ oh > 0, (w ** 2) @ oh, np.nan)}
        for name, col in (("win", "win"), ("win_resid", "win_resid"), ("mae_c", "eC"),
                          ("mae_g", "eG"), ("chews", "chews"), ("gcta", "gcta"),
                          ("truth", "h2"), ("gap", "gap")):
            out[name] = ((w * d[col].to_numpy()) @ oh) / safe
        if "gcta_uc" in d and d["gcta_uc"].notna().any():
            ok = d["gcta_uc"].notna().to_numpy()
            wu = w * ok
            su = wu @ oh
            out["gcta_uc"] = ((wu * np.nan_to_num(d["gcta_uc"].to_numpy())) @ oh) / np.where(su > 0, su, np.nan)
            euc = np.abs(np.nan_to_num(d["gcta_uc"].to_numpy()) - d["h2"].to_numpy())
            out["win_vs_uc"] = ((wu * (d["eC"].to_numpy() < euc)) @ oh) / np.where(su > 0, su, np.nan)
        bw = boot * w
        bwin = ((bw * d["win"].to_numpy()) @ oh) / np.where(bw @ oh > 0, bw @ oh, np.nan)
        out["win_lo"] = np.nanpercentile(bwin, 2.5, axis=0)
        out["win_hi"] = np.nanpercentile(bwin, 97.5, axis=0)
        return out

    per_prior = {}
    for _, f in acc.iterrows():
        per_prior[(f["pi0"], f["h2_scale"], f["q"])] = metrics(
            target_weights(f["pi0"], f["h2_scale"], f["q"]))
    bkey = (best["pi0"], best["h2_scale"], best["q"])
    bm = per_prior[bkey]
    rows = []
    for j, c in enumerate(all_cells):
        r, col = c.split("|")
        rates = np.array([m["win"][j] for m in per_prior.values()])
        los = np.array([m["win_lo"][j] for m in per_prior.values()])
        his = np.array([m["win_hi"][j] for m in per_prior.values()])
        if np.all(los > 0.5):
            rule = "2: closer to truth"
        elif np.all(his < 0.5):
            rule = "3: real loss"
        else:
            rule = "4: depends on prior"
        rows.append(dict(clump_bin=r, h2_bin=col, real_genes=int(real_cnt.loc[r, col]) if np.isfinite(real_cnt.loc[r, col]) else 0,
                         ess=bm["ess"][j], win=bm["win"][j], win_lo=bm["win_lo"][j],
                         win_hi=bm["win_hi"][j], win_min=np.nanmin(rates), win_max=np.nanmax(rates),
                         ci_min=np.nanmin(los), ci_max=np.nanmax(his),
                         win_resid=bm["win_resid"][j], rule=rule,
                         masked=bool(not (bm["ess"][j] >= 30) or not (real_cnt.loc[r, col] >= 30))))
    win_t = pd.DataFrame(rows)
    win_t.to_csv(OUT / "sim_strata_win.csv", index=False)
    bias_rows = []
    for j, c in enumerate(all_cells):
        r, col = c.split("|")
        row = dict(clump_bin=r, h2_bin=col, ess=bm["ess"][j], truth=bm["truth"][j],
                   mean_gcta=bm["gcta"][j], mean_chews=bm["chews"][j],
                   bias_gcta=bm["gcta"][j] - bm["truth"][j],
                   bias_chews=bm["chews"][j] - bm["truth"][j],
                   mae_gcta=bm["mae_g"][j], mae_chews=bm["mae_c"][j],
                   mae_diff=bm["mae_c"][j] - bm["mae_g"][j], gap=bm["gap"][j],
                   real_gap=real_gap.loc[r, col])
        if "gcta_uc" in bm:
            row.update(mean_gcta_uc=bm["gcta_uc"][j], bias_gcta_uc=bm["gcta_uc"][j] - bm["truth"][j],
                       win_vs_gcta_uc=bm["win_vs_uc"][j])
        bias_rows.append(row)
    bias_t = pd.DataFrame(bias_rows)
    bias_t.to_csv(OUT / "sim_strata_bias.csv", index=False)

    # win rate against pi0 in the PI's cells, other knobs at their best value
    curve = []
    for pi0 in np.round(np.arange(0, 1.0, 0.1), 1):
        m = metrics(target_weights(pi0, best["h2_scale"], best["q"]))
        for r in ("0", "1"):
            j = all_cells.index(f"{r}|{H2_LABELS[0]}")
            curve.append(dict(pi0=pi0, clump_bin=r, win=m["win"][j], lo=m["win_lo"][j],
                              hi=m["win_hi"][j], ess=m["ess"][j]))
    pd.DataFrame(curve).to_csv(OUT / "sim_win_vs_pi0.csv", index=False)

    with open(OUT / "calibration_summary.json", "w") as fh:
        json.dump(dict(best=dict(pi0=float(best["pi0"]), h2_scale=float(best["h2_scale"]),
                                 q=float(best["q"])), S_min=float(s_min), tolerance=float(tol),
                       n_acceptable=int(len(acc)), n_terms=int(n_terms),
                       calibration_passes=bool(calib_ok),
                       sign_agree=float(best["sign_agree"]), mad_gap=float(best["mad_gap"]),
                       tvd=float(best["tvd"])), fh, indent=2)

    show = win_t.assign(win=lambda t: (100 * t["win"]).round(0)).pivot(
        index="clump_bin", columns="h2_bin", values="win").reindex(ORDER5)[H2_LABELS]
    print("\n   best-prior win rate (% CHEWS closer to truth), half B by half-A stratum")
    print(show.to_string())
    print("\n   0-0.05 column:")
    print(win_t[win_t.h2_bin == H2_LABELS[0]][["clump_bin", "real_genes", "ess", "win", "win_lo",
                                                "win_hi", "win_min", "win_max", "ci_min", "ci_max",
                                                "win_resid", "rule"]].round(3).to_string(index=False))
    print(bias_t[bias_t.h2_bin == H2_LABELS[0]].round(4).to_string(index=False))
    return win_t, bias_t, fits


PARTS = {"a": part_a, "b": part_b, "c": part_c, "d": part_d, "e": part_e,
         "f": part_f, "g": part_g}


if __name__ == "__main__":
    import sys

    # no argument = the round-3 analyses exactly as before (A-D).  Round-4
    # parts are opt-in by letter:  python revision_analyses.py e f g
    for key in (sys.argv[1:] or ["a", "b", "c", "d"]):
        PARTS[key.lower()]()
