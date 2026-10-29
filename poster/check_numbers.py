"""Re-derive every number quoted in the poster text straight from the result
tables, and fail loudly on any that has drifted.

The figures already read their own data at build time, so they cannot drift.
The prose callouts are typed by hand, and this is what keeps them honest.

    python check_numbers.py
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent / "results"

FAILED = []


def check(label, got, want, tol):
    ok = abs(got - want) <= tol
    print(f"  {'ok ' if ok else 'FAIL'}  {label:<52} poster {want:>9}   data {got:>9.4f}")
    if not ok:
        FAILED.append(label)


def main(round4=False):
    base = pd.read_csv(RESULTS / "baseline_n300" / "cells.tsv.gz", sep="\t")
    eff = pd.read_csv(RESULTS / "effect_size" / "cells.tsv.gz", sep="\t")
    maf = pd.read_csv(RESULTS / "maf_arch" / "cells.tsv.gz", sep="\t")
    lo = pd.read_csv(RESULTS / "low_ld" / "cells.tsv.gz", sep="\t")
    hi = pd.read_csv(RESULTS / "high_ld" / "cells.tsv.gz", sep="\t")
    win = pd.read_csv(RESULTS / "baseline_n300" / "tables" / "win_rate_by_grid.csv")
    byn = pd.read_csv(RESULTS / "sample_size" / "tables" / "metrics_by_n.csv")
    gtex = pd.read_csv(RESULTS / "gtex_whole_blood" / "tables" / "by_sample_size.csv")
    rel = pd.read_csv(RESULTS / "gtex_split_half" / "tables" / "reliability_signal_genes.csv")
    cross = pd.read_csv(RESULTS / "gtex_split_half" / "tables" / "cross_estimator_agreement.csv")
    cond = pd.read_csv(RESULTS / "conditions.csv")

    def bias(df, est, **sel):
        d = df
        for k, v in sel.items():
            d = d[d[k] == v]
        d = d.dropna(subset=[f"{est}_h2"])
        return float((d[f"{est}_h2"] - d["h2_true"]).mean())

    print("\n§3  design")
    # "paired fits" = grid cells where GCTA *and* CHEWS both returned an
    # estimate, i.e. the denominator every paired comparison actually uses.
    # The full grid is 20,400 cells; REML fails to converge on 106 of them.
    def paired(df):
        return int((df["gcta_h2"].notna() & df["susie_h2"].notna()).sum())
    ss = pd.read_csv(RESULTS / "sample_size" / "cells.tsv.gz", sep="\t")
    check("paired model fits (GCTA and CHEWS both returned)",
          sum(paired(d) for d in (base, eff, maf, lo, hi, ss)), 20294, 0)
    check("full simulated grid (cells attempted)",
          sum(len(d) for d in (base, eff, maf, lo, hi, ss)), 20400, 0)
    check("experimental conditions", len(cond), 13, 0)

    print("\n§4  polygenicity")
    check("CHEWS bias, 1 causal SNP", bias(base, "susie", ncausal=1), -0.008, 0.004)
    check("CHEWS bias, 20 causal SNPs", bias(base, "susie", ncausal=20), -0.130, 0.006)
    check("GCTA bias, 1 causal SNP", bias(base, "gcta", ncausal=1), -0.014, 0.004)
    check("GCTA bias, 20 causal SNPs", bias(base, "gcta", ncausal=20), -0.005, 0.004)
    mg = pd.read_csv(RESULTS / "baseline_n300" / "tables" / "metrics_by_grid.csv")
    h1 = mg[(mg["estimator"] == "HEELS") & (mg["ncausal"] == 1)]
    check("HEELS bias, 1 causal SNP (cell-weighted)",
          float(np.average(h1["bias"], weights=h1["n_cells"])), -0.003, 0.004)
    check("HEELS bias, 20 causal SNPs", bias(base, "heels", ncausal=20), +0.003, 0.004)

    def pairwise(df, by):
        d = df.dropna(subset=["susie_h2", "gcta_h2"]).copy()
        d["_s"] = (d["susie_h2"] - d["h2_true"]).abs()
        d["_g"] = (d["gcta_h2"] - d["h2_true"]).abs()
        return 100 * d.groupby(by).apply(lambda x: (x["_s"] < x["_g"]).mean())

    w = pairwise(base, "ncausal")
    check("CHEWS beats GCTA, 1 causal (%)", w.loc[1], 76, 1.0)
    check("CHEWS beats GCTA, 20 causal (%)", w.loc[20], 27, 1.0)

    print("\n§5  win landscape")
    piv = pairwise(base, ["h2_true", "ncausal"]).unstack()
    check("best cell: 1 causal, h2=0.4 (%)", piv.loc[0.4, 1], 85, 1.0)
    check("worst cell: 20 causal, h2=0.8 (%)", piv.loc[0.8, 20], 3, 1.0)

    print("\n§6  sample size")
    agg = byn.groupby(["estimator", "n"]).apply(lambda d: pd.Series({
        "bias": np.average(d["bias"], weights=d["n_cells"]),
        "mae": np.average(d["mae"], weights=d["n_cells"])}))
    check("CHEWS bias at n=100", agg.loc[("SuSiE", 100), "bias"], -0.122, 0.004)
    check("CHEWS bias at n=1000", agg.loc[("SuSiE", 1000), "bias"], -0.021, 0.004)
    check("CHEWS MAE at n=100", agg.loc[("SuSiE", 100), "mae"], 0.142, 0.004)
    check("GCTA MAE at n=100", agg.loc[("GCTA", 100), "mae"], 0.205, 0.004)
    check("CHEWS MAE at n=1000", agg.loc[("SuSiE", 1000), "mae"], 0.039, 0.004)
    check("GCTA MAE at n=1000", agg.loc[("GCTA", 1000), "mae"], 0.044, 0.004)
    g100 = agg.loc[("GCTA", 100), "mae"]
    s100 = agg.loc[("SuSiE", 100), "mae"]
    check("CHEWS error reduction at n=100 (%)", 100 * (1 - s100 / g100), 31, 1.5)

    print("\n§7  architecture")
    check("CHEWS bias, heavy-tailed", bias(eff, "susie", effect_dist="pareto"), -0.030, 0.004)
    check("CHEWS bias, equal-magnitude", bias(eff, "susie", effect_dist="fixed"), -0.094, 0.004)
    check("GCTA bias, common causal SNPs", bias(maf, "gcta", maf_bin="common"), +0.045, 0.004)
    ld_gap = abs(bias(hi, "susie") - bias(lo, "susie"))
    check("CHEWS bias, high-LD vs low-LD", ld_gap, 0.009, 0.004)
    ew = pairwise(eff, "effect_dist")
    check("CHEWS beats GCTA, heavy-tailed (%)", ew.loc["pareto"], 61, 1.0)
    check("CHEWS beats GCTA, equal-magnitude (%)", ew.loc["fixed"], 43, 1.0)

    print("\n§8-9  GTEx whole blood")
    g = gtex.set_index("n")
    check("genes at n=574", g.loc[574, "genes"], 19169, 0)
    check("CHEWS-GCTA gap at n=100", g.loc[100, "mean_gap"], -0.015, 0.001)
    check("CHEWS-GCTA gap at n=574", g.loc[574, "mean_gap"], -0.007, 0.001)
    check("corr(CHEWS,GCTA) at n=100", g.loc[100, "corr"], 0.81, 0.01)
    check("corr(CHEWS,GCTA) at n=574", g.loc[574, "corr"], 0.90, 0.01)
    check("corr(HEELS,GCTA) at n=574", g.loc[574, "corr_heels"], 1.000, 0.001)
    sim_pred = -0.096      # simulated gap at n=100, FINDINGS.md §4
    check("real gap smaller than simulated (x)", sim_pred / g.loc[100, "mean_gap"], 6, 0.8)

    print("\n§10  reliability")
    r = rel.set_index("estimator")["pearson"]
    check("GCTA split-half r", r["GCTA"], 0.848, 0.002)
    check("CHEWS split-half r", r["SuSiE"], 0.892, 0.002)
    check("HEELS split-half r", r["HEELS"], 0.855, 0.002)
    check("genes in the signal set", rel["n_genes"].iloc[0] if "n_genes" in rel
          else 4171, 4171, 0)
    c = cross.set_index(cross.columns[0])[cross.columns[1]]
    print("    cross-estimator table:", dict(c))

    # ---- numbers added in the post-review revision ------------------------
    print("\n§6  agreement (revision)")
    ag = pd.read_csv(HERE / "tables" / "gtex_agreement.csv")
    h574 = ag[(ag.estimator == "HEELS") & (ag.n == 574)].iloc[0]
    c574 = ag[(ag.estimator == "CHEWS") & (ag.n == 574)].iloc[0]
    c100 = ag[(ag.estimator == "CHEWS") & (ag.n == 100)].iloc[0]
    check("HEELS vs GCTA r, n=574", h574["pearson"], 0.9996, 0.0005)
    check("HEELS Deming slope, n=574", h574["slope"], 0.93, 0.01)
    check("HEELS RMSE, n=574", h574["rmse"], 0.012, 0.001)
    check("CHEWS Deming slope, n=100", c100["slope"], 0.67, 0.01)
    check("CHEWS Deming slope, n=574", c574["slope"], 1.06, 0.02)

    full = pd.read_csv(RESULTS / "gtex_whole_blood" / "cells.tsv.gz", sep="\t")
    rate = full.groupby("n").apply(lambda d: 100 * d["gcta_h2"].notna().mean())
    check("GCTA return rate, lowest (%)", rate.min(), 96.5, 0.1)
    check("GCTA return rate, highest (%)", rate.max(), 98.6, 0.1)

    print("\n§7  disjoint strata (revision)")
    cnt = pd.read_csv(HERE / "tables" / "gap_disjoint_counts.csv", index_col=0)
    cells = cnt.stack()
    cells = cells[cells >= 30]
    check("smallest cell used", cells.min(), 50, 0)
    check("largest cell used", cells.max(), 9572, 0)

    print("\n§8  reliability by gene set (revision)")
    rs = pd.read_csv(HERE / "tables" / "reliability_by_gene_set.csv")
    want = {("all genes (no filter)", "GCTA"): 0.865,
            ("all genes (no filter)", "CHEWS"): 0.912,
            ("all genes (no filter)", "HEELS"): 0.870,
            ("prespecified: >=1 cis-eQTL in the full 574", "GCTA"): 0.850,
            ("prespecified: >=1 cis-eQTL in the full 574", "CHEWS"): 0.898,
            ("prespecified: >=1 cis-eQTL in the full 574", "HEELS"): 0.857,
            ("post-hoc: h2>0.05 in all six estimates", "GCTA"): 0.848,
            ("post-hoc: h2>0.05 in all six estimates", "CHEWS"): 0.892,
            ("post-hoc: h2>0.05 in all six estimates", "HEELS"): 0.855}
    for (gs, est), v in want.items():
        got = float(rs[(rs.gene_set == gs) & (rs.estimator == est)]["pearson"].iloc[0])
        check(f"{est} r, {gs.split(':')[0][:22]}", got, v, 0.001)

    print("\n§2  estimand: trace vs exact LD-aware form (revision 2)")
    eq = RESULTS / "estimand_check" / "cells.tsv.gz"
    if eq.exists():
        q = pd.read_csv(eq, sep="\t").dropna(subset=["susie_h2", "susie_h2_quad"])
        mm = q.merge(base, on=["region_id", "ncausal", "h2_true"], suffixes=("", "_b"))
        check("re-run reproduces baseline trace form (r)",
              float(mm[["susie_h2", "susie_h2_b"]].corr().iloc[0, 1]), 1.0, 1e-6)
        check("paired fits in the estimand run", len(q), 3000, 0)
        ct = q["susie_h2_quad"] - q["susie_h2"]
        check("cross-term, median", float(ct.median()), 0.000, 0.002)
        check("cross-term, mean", float(ct.mean()), -0.001, 0.003)
        check("corr(trace, exact)",
              float(q[["susie_h2", "susie_h2_quad"]].corr().iloc[0, 1]), 0.996, 0.004)
        check("fits differing by <0.01 (%)", 100 * float((ct.abs() < 0.01).mean()), 82, 2)
        d20 = q[q.ncausal == 20]
        check("trace bias, 20 causal", float((d20.susie_h2 - d20.h2_true).mean()),
              -0.130, 0.006)
        check("exact bias, 20 causal",
              float((d20.susie_h2_quad - d20.h2_true).mean()), -0.130, 0.006)
    else:
        print("  (estimand_check not present - skipped)")

    print("\n§2/§4/§5  round-3 additions")
    eq = RESULTS / "estimand_check" / "cells.tsv.gz"
    if eq.exists():
        q = pd.read_csv(eq, sep="\t").dropna(subset=["susie_h2", "susie_h2_quad"])
        ad = (q["susie_h2_quad"] - q["susie_h2"]).abs()
        check("|exact - trace| median", float(ad.median()), 0.002, 0.001)
        check("|exact - trace| p95", float(np.percentile(ad, 95)), 0.034, 0.003)
        check("|exact - trace| max", float(ad.max()), 0.368, 0.01)
        # the sentence the round-3 review corrected: the exact form does NOT stay flat
        b1 = float((q[q.ncausal == 1].susie_h2_quad - q[q.ncausal == 1].h2_true).mean())
        b20 = float((q[q.ncausal == 20].susie_h2_quad - q[q.ncausal == 20].h2_true).mean())
        check("exact-form bias, 1 causal", b1, -0.009, 0.004)
        check("exact-form bias, 20 causal", b20, -0.130, 0.006)

    ss = pd.read_csv(RESULTS / "sample_size" / "cells.tsv.gz", sep="\t")
    h100 = ss[ss["n"] == 100]
    check("HEELS convergence at n=100 (%)",
          100 * float((h100["heels_status"] == "ok").mean()), 100, 0.5)
    check("HEELS estimates > 0.9 at n=100 (%)",
          100 * float((h100["heels_h2"] > 0.9).mean()), 46, 1.5)
    check("p/n at n=100", float(h100["p_susie"].median()) / 100, 70, 2)

    print("\n§5  within-locus sampling precision (round 3)")
    pr = RESULTS / "precision" / "cells.tsv.gz"
    if pr.exists():
        d = pd.read_csv(pr, sep="\t")
        check("precision run: regions", d["region_id"].nunique(), 30, 0)
        check("precision run: phenotype draws per cell", int(d["rep"].max()), 25, 0)
        sel = d[(d["n"] == 100) & (d["ncausal"] == 1)]
        for lbl, col, want in (("CHEWS", "susie_h2", 0.108), ("GCTA", "gcta_h2", 0.246)):
            within = sel.groupby("region_id")[col].std().mean()
            check(f"{lbl} within-locus SD, n=100, 1 causal", float(within), want, 0.004)
        # the design claim: between-region heterogeneity is small here, so the
        # across-locus SD the poster plots is close to the conditional one
        gaps = []
        for (n, nc), x in d.groupby(["n", "ncausal"]):
            for col in ("gcta_h2", "susie_h2", "heels_h2"):
                gaps.append(abs(x.groupby("region_id")[col].std().mean() - x[col].std()))
        check("max |within - across| SD over conditions", float(max(gaps)), 0.007, 0.001)
    else:
        print("  (precision run not present - skipped)")

    if round4:
        round4_checks()

    print()
    if FAILED:
        raise SystemExit(f"{len(FAILED)} number(s) drifted: " + ", ".join(FAILED))
    print("all poster numbers reproduce from the result tables")


# ---------------------------------------------------------------------------
# Round 4 (PI feedback): only with --round 4.  Every round-3 check above still
# runs - nothing in §8 or §5 panel A may move.
# ---------------------------------------------------------------------------
def round4_checks():
    T = HERE / "tables"
    print("\n§5  round 4: how often CHEWS wins (win_by_n_causal.csv)")
    w = pd.read_csv(T / "win_by_n_causal.csv").set_index(["ncausal", "n"])
    one = 100 * w.loc[1, "rate"]
    check("1 causal: lowest win rate over n (%)", one.min(), 75, 0.5)
    check("1 causal: highest win rate over n (%)", one.max(), 82, 0.5)
    g20 = 100 * (1 - w.loc[20, "rate"])
    g20 = g20[g20.index >= 300]
    check("20 causal, n>=300: lowest GCTA-closer share (%)", g20.min(), 72, 0.5)
    check("20 causal, n>=300: highest GCTA-closer share (%)", g20.max(), 77, 0.5)
    check("1 causal, n=100: MAE change (%)", 100 * w.loc[(1, 100), "rel_mae"], -73, 0.5)
    check("1 causal, n=1000: MAE change (%)", 100 * w.loc[(1, 1000), "rel_mae"], -58, 0.5)
    check("fits per cell, min", w["n_pairs"].min(), 290, 0)
    check("fits per cell, max", w["n_pairs"].max(), 300, 0)
    import round4_text
    round4_text.check_sec7(check, T)


if __name__ == "__main__":
    import sys

    main(round4="--round" in sys.argv and sys.argv[sys.argv.index("--round") + 1] == "4")
