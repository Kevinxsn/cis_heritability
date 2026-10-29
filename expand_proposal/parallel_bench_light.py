# Node-level parallel-throughput test.
# The pipeline is embarrassingly parallel across genes/regions: each region is
# an independent task with no communication. We therefore run C identical
# single-threaded region-tasks concurrently, each pinned to its own core, and
# measure how throughput (regions/hour) scales as C goes 1 -> many cores on ONE
# Expanse-class node (2x AMD EPYC, 128 physical cores).
#
# Reported: wall time for the batch, throughput, speedup, parallel efficiency.

import sys, os, time, json, subprocess

SIM_DIR = "/projects/ps-renlab2/sux002/cis_heritability/simulations"
REGION  = "bench_regions/reg_22_kb100"   # p=3567
N, NCAUSAL, H2 = 500, 5, 0.2
PLINK = "/projects/ps-renlab2/sunan/anaconda3/envs/dsc180/bin/plink"
GCTA  = "/projects/ps-renlab2/sunan/anaconda3/envs/dsc180/bin/gcta64"

WORKER = r'''
import sys, os, time, json
sys.path.insert(0, "%s")
os.environ.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
                  MKL_NUM_THREADS="1", NUMEXPR_NUM_THREADS="1")
import numpy as np, statsmodels.api as sm, subprocess, pandas as pd
import simulations_test as sim
region, n, ncausal, h2, work, wid = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4]), sys.argv[5], sys.argv[6]
os.makedirs(work, exist_ok=True)
t0=time.time()
bim,L=sim.get_ld(region)
Z=sim.sim_geno(L,n)
b=sim.sim_beta(L,ncausal,h2)
gexpr=sim.sim_trait(np.dot(Z,b),h2)[0]
curr=f"w{wid}"
sim.get_gcta_hsq(Z,gexpr,os.path.join(work,curr),"%s")
betas=[];pvals=[];ses=[]
for i in range(Z.shape[1]):
    x=sm.add_constant(Z[:,i].reshape(-1,1)); r=sm.OLS(gexpr,x).fit()
    betas.append(r.params[1]);pvals.append(r.pvalues[1]);ses.append(r.bse[1])
ss=pd.DataFrame({"Gene":curr,"SNP":bim["snp"],"A1":bim["a0"],"A2":bim["a1"],"BETA":betas,"SE":ses,"P":pvals})
fn=os.path.join(work,curr+".txt"); ss.to_csv(fn,sep="\t",index=False)
ldg=os.path.join(work,curr+"_ld")
subprocess.run(["%s","--bfile",region,"--r","inter-chr","--ld-window-r2","0","--keep-allele-order","--extract",fn,"--out",ldg],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
subprocess.run(["Rscript","run_susie.R","-c",fn,"-l",ldg+".ld","-n",str(n),"-b",region+".bim","-o",fn],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
for f in [fn, ldg+".ld", ldg+".log"]:
    try: os.remove(f)
    except OSError: pass
print(json.dumps({"wid":wid,"cell":time.time()-t0}))
''' % (SIM_DIR, GCTA, PLINK)

def run_batch(C, ncore_phys=128):
    """Launch C worker processes concurrently, each pinned to a distinct core."""
    procs = []
    work = f"bench_work/par_C{C}"
    os.makedirs(work, exist_ok=True)
    t0 = time.time()
    for i in range(C):
        core = i % ncore_phys
        cmd = ["taskset", "-c", str(core), "python", "-c", WORKER,
               REGION, str(N), str(NCAUSAL), str(H2), work, str(i)]
        procs.append(subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True))
    cells = []
    for p in procs:
        out, _ = p.communicate()
        for line in out.strip().splitlines():
            try: cells.append(json.loads(line)["cell"])
            except Exception: pass
    wall = time.time() - t0
    return {"C": C, "wall": wall, "n_done": len(cells),
            "mean_cell": sum(cells)/len(cells) if cells else None,
            "max_cell": max(cells) if cells else None,
            "throughput_per_hr": C / wall * 3600}

def main():
    concurrencies = [int(x) for x in sys.argv[1].split(",")] if len(sys.argv) > 1 else [1,2,4,8,16,32,64,128]
    rows = []
    base = None
    for C in concurrencies:
        r = run_batch(C)
        if base is None: base = r["wall"]        # C=1 wall = single-core time for one cell
        r["speedup"] = (base * C) / r["wall"]     # ideal = C
        r["efficiency"] = r["speedup"] / C
        rows.append(r)
        print(f"C={C:>4}  wall={r['wall']:7.1f}s  mean_cell={r['mean_cell']:6.1f}s  "
              f"thru={r['throughput_per_hr']:8.1f}/hr  speedup={r['speedup']:6.1f}x  eff={r['efficiency']:.2f}", flush=True)
        with open("results/parallel_scaling_light.json", "w") as fh:
            json.dump(rows, fh, indent=2)
    print("wrote results/parallel_scaling_light.json")

if __name__ == "__main__":
    main()
