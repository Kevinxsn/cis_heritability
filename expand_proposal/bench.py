# Benchmark harness for the cis-heritability simulation pipeline.
# Reuses the exact functions from simulations_test.py and times each stage,
# so the numbers reflect the real homegrown code, not a re-implementation.
#
# It measures, for a single genomic region and sample size n:
#   - build_LD  : read PLINK, standardize, LD, regularize, Cholesky   (O(p^3))
#   - sim_geno  : sample MVN genotypes from L                          (O(n p^2))
#   - gcta      : build GRM + REML                                     (O(n^2 p) + REML)
#   - ols       : per-SNP OLS marginal sumstats                        (O(n p))
#   - plink_ld  : plink --r inter-chr pairwise LD                      (O(p^2))
#   - susie     : Rscript run_susie.R (IBSS fine-mapping)              (O(p^2 L))
# plus peak RSS for the Python process.

import sys, os, time, json, argparse, resource, gc
import numpy as np

# import the real pipeline
SIM_DIR = "/projects/ps-renlab2/sux002/cis_heritability/simulations"
sys.path.insert(0, SIM_DIR)
import simulations_test as sim   # noqa

def peak_rss_gb():
    # ru_maxrss is in kilobytes on Linux
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024.0 * 1024.0)

def time_region(region_prefix, n, ncausal, h2, work, plink_exe, gcta_exe, geno_full_bim_prefix, chrom):
    """Run one full grid-cell on one region, timing each stage. Returns dict."""
    rec = {"region": os.path.basename(region_prefix), "n": n, "ncausal": ncausal, "h2": h2}
    # p = number of SNPs in the region
    with open(region_prefix + ".bim") as fh:
        p = sum(1 for _ in fh)
    rec["p"] = p

    t0 = time.time()
    bim, L = sim.get_ld(region_prefix)          # includes Cholesky
    rec["t_build_ld"] = time.time() - t0

    t0 = time.time()
    Z = sim.sim_geno(L, n)
    rec["t_sim_geno"] = time.time() - t0

    t0 = time.time()
    b = sim.sim_beta(L, ncausal, h2)
    gexpr = sim.sim_trait(np.dot(Z, b), h2)[0]
    rec["t_sim_pheno"] = time.time() - t0

    curr = f"{os.path.basename(region_prefix)}_{ncausal}_{h2:.3f}_{n}"
    hsq_dir = os.path.join(work, "hsq"); os.makedirs(hsq_dir, exist_ok=True)
    t0 = time.time()
    h2g, hsq_se, hsq_p = sim.get_gcta_hsq(Z, gexpr, os.path.join(hsq_dir, curr), gcta_exe)
    rec["t_gcta"] = time.time() - t0
    rec["gcta_h2"] = float(h2g)

    # SuSiE lane, timed in sub-stages by replicating build_susie_sumstats internals
    import statsmodels.api as sm, subprocess, pandas as pd
    ss_dir = os.path.join(work, "sumstats"); os.makedirs(ss_dir, exist_ok=True)
    ld_dir = os.path.join(work, "ld"); os.makedirs(ld_dir, exist_ok=True)

    t0 = time.time()
    betas, pvals, std_errs = [], [], []
    for i in range(Z.shape[1]):
        x = sm.add_constant(Z[:, i].reshape(-1, 1))
        res = sm.OLS(gexpr, x).fit()
        betas.append(res.params[1]); pvals.append(res.pvalues[1]); std_errs.append(res.bse[1])
    ss = pd.DataFrame({"Gene": curr, "SNP": bim["snp"], "A1": bim["a0"], "A2": bim["a1"],
                       "BETA": betas, "SE": std_errs, "P": pvals})
    fn = os.path.join(ss_dir, curr + ".txt")
    ss.to_csv(fn, sep="\t", index=False, header=True)
    rec["t_ols"] = time.time() - t0

    ld_gene = os.path.join(ld_dir, curr)
    t0 = time.time()
    subprocess.run([plink_exe, "--bfile", region_prefix, "--r", "inter-chr", "--ld-window-r2", "0",
                    "--keep-allele-order", "--extract", fn, "--out", ld_gene],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    rec["t_plink_ld"] = time.time() - t0
    rec["ld_file_mb"] = os.path.getsize(ld_gene + ".ld") / 1e6 if os.path.exists(ld_gene + ".ld") else 0.0

    t0 = time.time()
    subprocess.run(["Rscript", "run_susie.R", "-c", fn, "-l", ld_gene + ".ld",
                    "-n", str(n), "-b", region_prefix + ".bim", "-o", fn],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    rec["t_susie"] = time.time() - t0

    rec["t_total_cell"] = sum(rec[k] for k in
                              ["t_build_ld", "t_sim_geno", "t_sim_pheno", "t_gcta",
                               "t_ols", "t_plink_ld", "t_susie"])
    rec["peak_rss_gb"] = peak_rss_gb()

    # cleanup per-cell files to avoid filling disk
    for f in [fn, ld_gene + ".ld", ld_gene + ".log", ld_gene + ".nosex"]:
        try: os.remove(f)
        except OSError: pass
    del Z, L, b, gexpr; gc.collect()
    return rec

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--regions", nargs="+", required=True, help="region plink prefixes")
    ap.add_argument("--n", type=int, default=500)
    ap.add_argument("--ncausal", type=int, default=5)
    ap.add_argument("--h2", type=float, default=0.2)
    ap.add_argument("--work", default="bench_work")
    ap.add_argument("--out", required=True)
    ap.add_argument("--plink", default="plink")
    ap.add_argument("--gcta", default="gcta64")
    ap.add_argument("--ns", nargs="*", type=int, help="optional list of n to sweep (overrides --n)")
    args = ap.parse_args()

    os.makedirs(args.work, exist_ok=True)
    ns = args.ns if args.ns else [args.n]
    rows = []
    for region in args.regions:
        for n in ns:
            t0 = time.time()
            rec = time_region(region, n, args.ncausal, args.h2, args.work,
                              args.plink, args.gcta, None, "22")
            rec["wall"] = time.time() - t0
            rows.append(rec)
            print(f"[done] p={rec['p']:>6} n={n:>4}  build_ld={rec['t_build_ld']:6.2f}s "
                  f"gcta={rec['t_gcta']:6.2f}s susie={rec['t_susie']:6.2f}s "
                  f"cell={rec['t_total_cell']:6.2f}s rss={rec['peak_rss_gb']:.2f}GB", flush=True)
            with open(args.out, "w") as fh:
                json.dump(rows, fh, indent=2)
    print("wrote", args.out)

if __name__ == "__main__":
    main()
