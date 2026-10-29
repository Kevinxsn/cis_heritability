# Cost model + budget in Expanse shared-partition core-hours, with memory-equivalent
# billing (nodes are 2 GB/core -> a job bills max(cores, ceil(mem_GB/2)) * wall_hours).
import json, numpy as np, math

ps  = json.load(open("results/psweep_n500.json"))
fg  = None
# ground-truth full-grid (from fullgrid_bench): p=1544, 30 cells, grid_s measured
GT = {"p":1544, "ncells":30, "grid_s":218.64, "per_cell_s":7.288}

p    = np.array([r["p"] for r in ps], float)
cell = np.array([r["t_total_cell"] for r in ps], float)
rss  = np.array([r["peak_rss_gb"] for r in ps], float)

# power-law fits
b_t, a_t = np.polyfit(np.log(p), np.log(cell), 1); a_t = math.exp(a_t)
b_m, a_m = np.polyfit(np.log(p), np.log(rss),  1); a_m = math.exp(a_m)
def wall_s(pp):  return a_t*pp**b_t         # single-core seconds per model
def mem_gb(pp):  return a_m*pp**b_m
def mem_cores(pp): return max(1, math.ceil(mem_gb(pp)/2.0))   # memory-equivalent cores
def corehr_per_model(pp): return mem_cores(pp)*wall_s(pp)/3600.0

print("== fits ==")
print(f"wall_s(p)  = {a_t:.4g} * p^{b_t:.2f}")
print(f"mem_gb(p)  = {a_m:.4g} * p^{b_m:.2f}")

print("\n== cross-validation against ground-truth full grid (p=1544, 30 cells) ==")
pred = 30*wall_s(1544)
print(f"predicted grid wall = {pred:.1f}s ; observed = {GT['grid_s']:.1f}s ; error = {(pred-GT['grid_s'])/GT['grid_s']*100:+.1f}%")

print("\n== per-model cost + memory billing by p ==")
print(f"{'p':>6} {'wall_s':>8} {'mem_GB':>7} {'mem_cores':>9} {'core-hr/model':>13}")
for pp in [1000,1544,3000,3567,5000,8000,8227,17287]:
    print(f"{pp:>6} {wall_s(pp):>8.1f} {mem_gb(pp):>7.2f} {mem_cores(pp):>9d} {corehr_per_model(pp):>13.4f}")

# ---------- campaign definitions ----------
P_SIM = 3000     # representative simulated cis-region
P_REAL= 3000     # representative real cis-window (1000G, ~1 Mb)
# simulation grid per locus (fully crossed upper design)
GRID = dict(n=5, ncausal=5, h2=6, dist=3, LD=3, MAF=2)
models_per_locus = np.prod(list(GRID.values()))    # 2700
LOCI = 1000
def sim_campaign_corehr(p=P_SIM, loci=LOCI, models_per=models_per_locus):
    return loci*models_per*corehr_per_model(p)

print("\n== SIMULATION campaign (1 full sweep) ==")
one_sweep = sim_campaign_corehr()
print(f"models/locus={models_per_locus}, loci={LOCI} -> {LOCI*models_per_locus:,} models ; {one_sweep:,.0f} core-hours/sweep")

# ---------- annual plan ----------
rows=[]
def add(name, models, p, note):
    ch = models*corehr_per_model(p)
    rows.append((name, models, p, mem_cores(p), ch, note)); return ch

# 1) simulation method development: several full sweeps + many partial
add("Simulation: full architecture sweeps (method development)", 6*LOCI*models_per_locus, P_SIM, "6 sweeps as code/estimators evolve")
add("Simulation: targeted single-axis campaigns", 40*LOCI*30, P_SIM, "40 focused campaigns x 30-cell grids x 1000 loci")
add("Simulation: high-memory large-window stress tests", 5*LOCI*30, 8227, "wide 1 Mb+ windows; memory-billed at ceil(3.3/2)=2 cores")
# 2) real data GTEx x methods
GTEX = 49*20000
add("Real data: GTEx cis-h2 (GCTA)", GTEX, P_REAL, "49 tissues x ~20k genes")
add("Real data: GTEx cis-h2 (SuSiE)", GTEX, P_REAL, "49 tissues x ~20k genes")
add("Real data: GTEx cis-h2 (HEELS, 3rd method)", GTEX, P_REAL, "planned benchmark method")
# 3) multi-ancestry reference panels
add("Multi-ancestry replication (5 panels, sim+real subset)", 5*(LOCI*30 + 20000), P_REAL, "EUR/AFR/EAS/SAS/AMR")
# 4) GWAS-trait fine-mapping
add("GWAS-trait locus fine-mapping", 100*1500, P_REAL, "100 traits x ~1,500 loci")
# 5) extrapolation to large scRNA-seq atlas (inflation target)
add("Large scRNA-seq atlas extension (e.g. TenK10K-scale)", 30*20000, P_REAL, "~30 cell states/contexts x 20k genes")

subtotal = sum(r[4] for r in rows)
retries  = 0.10*subtotal   # convergence retries
margin   = 0.20*subtotal   # safety margin
total    = subtotal+retries+margin

print("\n== ANNUAL BUDGET (core-hours, memory-equivalent billing) ==")
for name,models,pp,mc,ch,note in rows:
    print(f"  {ch:>10,.0f} ch | {models:>10,} models @ p={pp} ({mc} core-eq) | {name}")
print(f"  {retries:>10,.0f} ch | convergence retries (+10%)")
print(f"  {margin:>10,.0f} ch | safety margin (+20%)")
print(f"  {'-'*10}")
print(f"  {total:>10,.0f} ch  TOTAL REQUEST")

json.dump({
 "fit":{"wall_a":a_t,"wall_b":round(b_t,2),"mem_a":a_m,"mem_b":round(b_m,2)},
 "xval":{"pred_grid_s":round(pred,1),"obs_grid_s":GT["grid_s"],"err_pct":round((pred-GT['grid_s'])/GT['grid_s']*100,1)},
 "per_model":{str(pp):{"wall_s":round(wall_s(pp),1),"mem_gb":round(mem_gb(pp),2),
              "mem_cores":mem_cores(pp),"corehr":round(corehr_per_model(pp),4)}
              for pp in [1000,1544,3000,3567,5000,8000,8227,17287]},
 "models_per_locus":int(models_per_locus),
 "one_full_sweep_corehr":round(one_sweep),
 "budget":[{"name":n,"models":int(m),"p":pp,"core_eq":mc,"core_hours":round(ch)} for n,m,pp,mc,ch,note in rows],
 "retries_ch":round(retries),"margin_ch":round(margin),"total_ch":round(total),
}, open("results/budget.json","w"), indent=2)
print("\nwrote results/budget.json")
