# Turn measured benchmark JSON into proposal figures + an SU cost model.
# Design: Okabe-Ito colorblind-safe categorical palette, thin marks, one axis,
# legend + selective direct labels, log-log where power-laws are the story.

import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

os.makedirs("figures", exist_ok=True)

# ---- Okabe-Ito colorblind-safe palette ----
OI = {"blue":"#0072B2","orange":"#E69F00","green":"#009E73","vermillion":"#D55E00",
      "purple":"#CC79A7","sky":"#56B4E9","yellow":"#F0E442","black":"#222222"}
plt.rcParams.update({
    "figure.dpi":130, "savefig.dpi":200, "font.size":11,
    "axes.spines.top":False, "axes.spines.right":False,
    "axes.grid":True, "grid.alpha":0.25, "grid.linewidth":0.6,
    "axes.titleweight":"bold", "axes.titlesize":12, "legend.frameon":False,
    "font.family":"DejaVu Sans",
})

def load(f): return json.load(open(f))
ps   = load("results/psweep_n500.json")
ns   = load("results/nsweep_p3567.json")
par  = load("results/parallel_scaling.json")
parL = load("results/parallel_scaling_light.json")

def powfit(x, y):
    x=np.asarray(x,float); y=np.asarray(y,float)
    b,a=np.polyfit(np.log(x),np.log(y),1); return b, np.exp(a)

# =====================================================================
# FIG 1: runtime + memory vs #SNPs (p)  -- the dominant cost driver
# =====================================================================
p   = [r["p"] for r in ps]
tot = [r["t_total_cell"] for r in ps]
rss = [r["peak_rss_gb"] for r in ps]
stages = [("t_build_ld","Build LD + Cholesky",OI["blue"]),
          ("t_susie","SuSiE fine-mapping",OI["vermillion"]),
          ("t_plink_ld","PLINK pairwise LD",OI["orange"]),
          ("t_gcta","GCTA REML",OI["green"]),
          ("t_ols","Per-SNP OLS",OI["purple"])]

fig, ax = plt.subplots(1,2, figsize=(11,4.4))
# left: per-stage runtime vs p (log-log)
for key,label,c in stages:
    ax[0].plot(p,[r[key] for r in ps],"o-",color=c,lw=2,ms=6,label=label)
b_tot,_=powfit(p,tot)
ax[0].plot(p,tot,"s--",color=OI["black"],lw=2,ms=7,label=f"Total per model (∝ p^{b_tot:.1f})")
ax[0].set_xscale("log"); ax[0].set_yscale("log")
ax[0].set_xlabel("SNPs in cis-region, p"); ax[0].set_ylabel("wall-clock seconds (1 core)")
ax[0].set_title("Runtime is driven by region size p")
ax[0].legend(fontsize=6.5, loc="upper left", handlelength=1.4, handletextpad=0.4,
             labelspacing=0.3, borderpad=0.3, framealpha=0.85, edgecolor="none")
for a in ax:
    a.set_xticks([1544,3567,8227,17287])
    a.set_xticklabels(["1.5k","3.6k","8.2k","17k"])
    a.xaxis.set_minor_formatter(plt.NullFormatter())
# right: memory vs p
b_rss,a_rss=powfit(p,rss)
ax[1].plot(p,rss,"o-",color=OI["blue"],lw=2,ms=7,label=f"Peak RAM (∝ p^{b_rss:.1f})")
xx=np.linspace(min(p),max(p),50)
ax[1].plot(xx,a_rss*xx**b_rss,":",color=OI["black"],lw=1.2,alpha=0.6)
for xi,yi in zip(p,rss):
    ax[1].annotate(f"{yi:.1f} GB",(xi,yi),textcoords="offset points",xytext=(6,-10),fontsize=8,color=OI["black"])
ax[1].set_xscale("log"); ax[1].set_yscale("log")
ax[1].set_xlabel("SNPs in cis-region, p"); ax[1].set_ylabel("peak resident memory (GB)")
ax[1].set_title("Memory footprint ≈ O(p²)")
ax[1].legend(fontsize=9, loc="upper left")
fig.suptitle("Per-model cost scales with the number of SNPs fine-mapped (n=500 individuals)", fontsize=12, y=1.02, fontweight="bold")
ax[0].text(0.5,-0.30,"A",transform=ax[0].transAxes,ha="center",va="top",fontsize=17,fontweight="bold")
ax[1].text(0.5,-0.30,"B",transform=ax[1].transAxes,ha="center",va="top",fontsize=17,fontweight="bold")
fig.tight_layout(); fig.savefig("figures/fig1_scaling_p.png", bbox_inches="tight"); plt.close(fig)

# =====================================================================
# FIG 2: runtime vs sample size n (fixed region) -- n is NOT the driver
# =====================================================================
nn=[r["n"] for r in ns]
fig,ax=plt.subplots(figsize=(6.2,4.4))
ax.plot(nn,[r["t_total_cell"] for r in ns],"s--",color=OI["black"],lw=2,ms=7,label="Total per model")
ax.plot(nn,[r["t_gcta"] for r in ns],"o-",color=OI["green"],lw=2,ms=6,label="GCTA REML (∝ n²)")
ax.plot(nn,[r["t_susie"] for r in ns],"o-",color=OI["vermillion"],lw=2,ms=6,label="SuSiE (n-independent)")
ax.plot(nn,[r["t_build_ld"] for r in ns],"o-",color=OI["blue"],lw=2,ms=6,label="Build LD (n-independent)")
ax.set_xlabel("sample size, n (individuals)"); ax.set_ylabel("wall-clock seconds (1 core)")
ax.set_title("At fixed region (p=3,567), sample size\nbarely moves total cost")
ax.legend(fontsize=9); ax.set_ylim(0,None)
fig.tight_layout(); fig.savefig("figures/fig2_scaling_n.png", bbox_inches="tight"); plt.close(fig)

# =====================================================================
# FIG 3: within-node parallel throughput + efficiency
# =====================================================================
def col(o,k): return [r[k] for r in o]
fig,ax=plt.subplots(1,2,figsize=(11,4.4))
ax[0].plot(col(parL,"C"),col(parL,"throughput_per_hr"),"o-",color=OI["green"],lw=2,ms=7,label="Light task (p=1,544, 0.25 GB)")
ax[0].plot(col(par,"C"),col(par,"throughput_per_hr"),"o-",color=OI["vermillion"],lw=2,ms=7,label="Heavy task (p=3,567, 0.8 GB)")
ax[0].set_xscale("log",base=2); ax[0].set_xlabel("concurrent tasks per node (cores used)")
ax[0].set_ylabel("throughput (models / hour / node)")
ax[0].set_title("More cores → more throughput,\nwith memory-bandwidth diminishing returns")
ax[0].legend(fontsize=9, loc="upper left")
ax[0].xaxis.set_major_formatter(ScalarFormatter()); ax[0].xaxis.set_minor_formatter(plt.NullFormatter())
# efficiency
C=col(parL,"C")
ax[1].plot(C,[1]*len(C),":",color=OI["black"],lw=1.2,label="ideal linear")
ax[1].plot(col(parL,"C"),col(parL,"efficiency"),"o-",color=OI["green"],lw=2,ms=7,label="Light task")
ax[1].plot(col(par,"C"),col(par,"efficiency"),"o-",color=OI["vermillion"],lw=2,ms=7,label="Heavy task")
ax[1].set_xscale("log",base=2); ax[1].set_xlabel("concurrent tasks per node (cores used)")
ax[1].set_ylabel("parallel efficiency (speedup / cores)")
ax[1].set_title("Efficiency depends on per-task memory footprint")
ax[1].legend(fontsize=9); ax[1].set_ylim(0,1.6)
ax[1].xaxis.set_major_formatter(ScalarFormatter()); ax[1].xaxis.set_minor_formatter(plt.NullFormatter())
ax[1].annotate("super-linear:\nboost clock + warm cache",(2,1.48),textcoords="offset points",
               xytext=(14,-6),fontsize=7.5,color=OI["black"],alpha=0.8)
fig.suptitle("Single Expanse-class node (2× AMD EPYC, 128 cores): embarrassingly-parallel across genes",
             fontsize=12,y=1.03,fontweight="bold")
fig.tight_layout(); fig.savefig("figures/fig3_parallel.png", bbox_inches="tight"); plt.close(fig)

# =====================================================================
# SU COST MODEL  (Expanse: 1 SU = 1 core-hour; standard node = 128 cores)
# =====================================================================
# per-model (one grid cell) single-core seconds as a function of p, from measured data
p_arr=np.array([r["p"] for r in ps]); cell_arr=np.array([r["t_total_cell"] for r in ps])
b_cell,a_cell=powfit(p_arr,cell_arr)
def cell_sec(pp): return a_cell*pp**b_cell
def su_per_model(pp): return cell_sec(pp)/3600.0   # core-hours = SU

# grid definition
GRID={"n (sample size)":[100,300,500,800,1000],
      "ncausal (causal SNPs)":[1,3,5,10,20],
      "h2 (heritability)":[0.05,0.1,0.2,0.4,0.6,0.8],
      "effect-size dist":["normal","fixed","pareto"],
      "LD architecture":["regular","high-LD","low-LD"],
      "MAF bin":["common","rare"]}

P_TYP=3000     # representative cis-region SNP count
scenarios=[]
# 1) one baseline campaign: 1000 regions x (5 ncausal x 6 h2) cells, n=500
n_models_baseline=1000*5*6
scenarios.append(("Baseline sim campaign\n(1000 regions × 5 ncausal × 6 h²)",n_models_baseline))
# 2) full factorial-ish: 1000 regions x ncausal(5) x h2(6) x dist(3) x LD(3) x MAF(2)
n_models_full=1000*5*6*3*3*2
scenarios.append(("Full architecture sweep\n(1000 regions × all parameter axes)",n_models_full))
# 3) real data GTEx: 49 tissues x ~20000 genes (one model each)
n_models_gtex=49*20000
scenarios.append(("GTEx real-data pass\n(49 tissues × 20k genes)",n_models_gtex))
# 4) 100 GWAS traits fine-mapped genome-wide (~1500 loci each)
n_models_gwas=100*1500
scenarios.append(("100 GWAS traits\n(~1,500 loci each)",n_models_gwas))

labels=[s[0] for s in scenarios]
su=[m*su_per_model(P_TYP) for _,m in scenarios]
node_days=[s/128/24 for s in su]

fig,ax=plt.subplots(figsize=(8.6,4.6))
colors=[OI["sky"],OI["blue"],OI["green"],OI["orange"]]
bars=ax.barh(range(len(labels)),su,color=colors,height=0.62)
ax.set_yticks(range(len(labels))); ax.set_yticklabels(labels,fontsize=9)
ax.invert_yaxis(); ax.set_xscale("log")
ax.set_xlabel("Service Units (core-hours) — log scale")
for i,(bar,s,nd) in enumerate(zip(bars,su,node_days)):
    ax.text(s*1.05,i,f"{s:,.0f} SU  ({nd:.1f} node-days)",va="center",fontsize=8.5,color=OI["black"])
ax.set_title("Estimated Expanse Service Units per campaign\n(measured per-model cost × number of models, p≈3,000)",fontsize=11)
ax.set_xlim(min(su)*0.5, max(su)*6)
fig.tight_layout(); fig.savefig("figures/fig4_su_cost.png", bbox_inches="tight"); plt.close(fig)

# =====================================================================
# FIG 5 (combined): A = runtime vs n ; B = Service Units per campaign
# =====================================================================
fig,ax=plt.subplots(1,2,figsize=(11,6.2))
# --- Panel A: runtime vs sample size n ---
ax[0].plot(nn,[r["t_total_cell"] for r in ns],"s--",color=OI["black"],lw=2,ms=7,label="Total per model")
ax[0].plot(nn,[r["t_gcta"] for r in ns],"o-",color=OI["green"],lw=2,ms=6,label="GCTA REML (∝ n²)")
ax[0].plot(nn,[r["t_susie"] for r in ns],"o-",color=OI["vermillion"],lw=2,ms=6,label="SuSiE (n-independent)")
ax[0].plot(nn,[r["t_build_ld"] for r in ns],"o-",color=OI["blue"],lw=2,ms=6,label="Build LD (n-independent)")
ax[0].set_xlabel("sample size, n (individuals)"); ax[0].set_ylabel("wall-clock seconds (1 core)")
ax[0].set_title("At fixed region (p=3,567),\nsample size barely moves total cost")
ax[0].set_ylim(0,None)
ax[0].legend(fontsize=7.5, loc="upper center", bbox_to_anchor=(0.5,-0.30), ncol=2,
             handlelength=1.5, handletextpad=0.4, columnspacing=1.2, labelspacing=0.3,
             borderpad=0.3, framealpha=0, edgecolor="none")
# --- Panel B: Service Units per campaign ---
short=["Baseline sim\n(1k loci × grid)","Full sweep\n(all axes)",
       "GTEx pass\n(49 tissues)","100 GWAS\ntraits"]
bars=ax[1].barh(range(len(short)),su,color=colors,height=0.62)
ax[1].set_yticks(range(len(short))); ax[1].set_yticklabels(short,fontsize=8.5)
ax[1].invert_yaxis(); ax[1].set_xscale("log")
ax[1].set_xlabel("Service Units (core-hours) — log scale")
for i,(s,nd) in enumerate(zip(su,node_days)):
    ax[1].text(s*1.06,i,f"{s:,.0f} SU",va="center",fontsize=8,color=OI["black"])
ax[1].set_title("Estimated Expanse Service Units\nper campaign (p≈3,000)")
ax[1].set_xlim(min(su)*0.5, max(su)*7)
ax[0].text(0.5,-0.62,"A",transform=ax[0].transAxes,ha="center",va="top",fontsize=17,fontweight="bold")
ax[1].text(0.5,-0.62,"B",transform=ax[1].transAxes,ha="center",va="top",fontsize=17,fontweight="bold")
fig.suptitle("Sample size is not the cost driver; total cost is set by number of models run",
             fontsize=12,y=1.02,fontweight="bold")
fig.tight_layout(); fig.savefig("figures/fig5_n_and_su.png", bbox_inches="tight"); plt.close(fig)

# =====================================================================
# dump numbers for the writeup
# =====================================================================
out={
 "fit_total_cell_exponent":round(b_cell,2),
 "fit_total_cell_coef":a_cell,
 "fit_memory_exponent":round(powfit(p_arr,np.array([r['peak_rss_gb'] for r in ps]))[0],2),
 "cell_sec_at_p":{str(pp):round(cell_sec(pp),1) for pp in [1000,3000,5000,8000]},
 "su_per_model_at_p":{str(pp):round(su_per_model(pp),5) for pp in [1000,3000,5000,8000]},
 "scenarios":[{"name":l.replace("\n"," "),"n_models":m,"SU":round(s), "node_days":round(nd,2)}
              for (l,(_,m)),s,nd in zip(zip(labels,scenarios),su,node_days)],
 "parallel_heavy":par,"parallel_light":parL,
 "grid":GRID,
 "ld_file_mb":{str(r["p"]):round(r["ld_file_mb"],1) for r in ps},
}
json.dump(out,open("results/numbers.json","w"),indent=2)
print("Total-per-model exponent p^%.2f ; memory exponent p^%.2f"%(b_cell,out["fit_memory_exponent"]))
print("SU per model @ p=3000:",out["su_per_model_at_p"]["3000"])
for sc in out["scenarios"]:
    print(f"  {sc['n_models']:>10,} models  {sc['SU']:>10,} SU  {sc['node_days']:>7.1f} node-days  | {sc['name']}")
print("wrote figures/ + results/numbers.json")
