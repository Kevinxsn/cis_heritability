# ACCESS Maximize — Project Documentation
## Fine-mapping-based estimation of cis-heritability (Amariuta Lab)
**Contributor:** Sunan Xu · **Prepared:** July 2026 · all numbers below are *measured*, not estimated, on SDSC AMD-EPYC hardware (the same Zen2/Zen3 generation as Expanse compute nodes). Reproduce with the scripts in this folder (`bench.py`, `parallel_bench.py`, `analyze.py`).

> **How to read this doc:** everything the PI's checklist asked for is a labeled section below. Tables and the four figures in `figures/` are drop-in ready for the proposal. Section 10 gives the Service-Unit (SU) budget built by the "small-experiment × multiplier" method ACCESS recommends.

---

## 1. One-paragraph project description
We benchmark **two estimators of *cis*-heritability of gene expression** — GCTA/REML (individual-level, the field-standard benchmark) and a **SuSiE fine-mapping-based estimator** we implemented (summary-statistic based) — to determine *when* fine-mapping gives more reliable *cis*-h² than GCTA. Because the ground-truth h² is never known in real data, the core of the project is a **simulation engine**: it draws a realistic local LD structure from 1000 Genomes, simulates gene expression with a *known* true h² under many genetic architectures, runs both estimators, and compares their error. We then confirm the simulation-derived patterns on real GTEx expression data. This is a **method-evaluation project**, so the compute cost is dominated by running the same pipeline across a large grid of architectures and thousands of replicate loci to obtain confidence intervals on bias, power, and type-I error.

---

## 2. "Homegrown" code we developed, and why it was necessary
Roughly **6,000 lines across 32 scripts** (Python + R + shell). It is small in size but is the scientific core of the project — no existing tool performs this end-to-end simulate-and-compare workflow, which is why it had to be written in-house.

| Component | Lang | What it does | Why off-the-shelf tools don't cover it |
|---|---|---|---|
| Simulation engine (`simulations_test.py`, `simulation_effect_size.py`, `simulation_high_al_freq.py`, high/low-LD variants) | Python | Draws a real *cis* region, builds & regularizes the LD matrix, Cholesky-samples genotypes, simulates causal effects to a **known** target h², runs both estimators, records error | No public tool simulates expression from a *fixed known cis-h²* **and** benchmarks GCTA vs. fine-mapping in one loop |
| SuSiE→h² estimator (`run_susie.R`) | R | Runs summary-statistic fine-mapping (`susie_rss`), performs allele QC / flip harmonization, and converts posterior effects into an h² estimate `Σ(μ² + σ²)` with an analytic SE | SuSiE returns fine-mapping posteriors, **not** a heritability estimate or its SE — that derivation is ours |
| GCTA harness | Python | Builds the GRM directly from the simulated genotype matrix, writes GCTA-format inputs, runs REML, parses results | Gluing simulated matrices into GCTA's file formats at scale is custom |
| Grid + retry/■job orchestration (`simulation_*_job*.sh`) | shell | Sweeps the parameter grid, manages temp dirs, retries non-converged replicates, cleans up | Bespoke to our parameter design |

The design principle we had to enforce ourselves (see lab-meeting slides 19): **the LD used to *simulate* genotypes must be the exact same ordered LD given to the fine-mapper**, or the comparison is uninterpretable. That alignment constraint is baked into our code and is not something a generic pipeline provides.

---

## 3. Third-party software required
| Software | Version | Role | License / source |
|---|---|---|---|
| **GCTA** | 1.94.1 | REML variance-component h² (benchmark estimator) | bioconda |
| **SuSiE (susieR)** | CRAN | summary-statistic fine-mapping (`susie_rss`, IBSS) | R / CRAN |
| **PLINK** | 1.90 | extract *cis* regions, compute pairwise LD | open source |
| **Python scientific stack** | NumPy, SciPy, pandas, statsmodels, scikit-learn, pandas-plink | LD/Cholesky, MVN genotype sampling, per-SNP OLS, I/O | conda-forge |
| **R stack** | Rfast, data.table, dplyr, optparse, ggplot2 | fast matrix ops + I/O inside `run_susie.R` | CRAN |

All are installed in a single conda environment; nothing requires a GPU. The workload is **CPU + memory-bandwidth bound** (dense linear algebra on p×p LD matrices), which maps directly onto Expanse's standard compute nodes.

---

## 4. Storage needs, code size, and file-set profile
| Category | Size | File-set shape | Notes |
|---|---|---|---|
| Homegrown code | ~0.5 MB (6k lines) | few small text files | version-controlled |
| Reference genotypes (1000G EUR, chr1–22) | **14 GB** | 22 × (large `.bed` + `.bim`), few large files | read-only, reused by every job |
| Per-run result tables | ~0.1–0.3 MB each | many small append-only text files | one summary row per (region × architecture) |
| **Transient LD matrices** (`*.ld`) | **0.1–12 GB *each*** | O(p²), written then deleted per model | **the dominant I/O object** — see §8 |
| Transient GRM / sumstats / phenotype | KB–MB each | many small files, created & deleted per model | churns the metadata layer |

**Data movement:** the 14 GB reference is staged once to node-local NVMe (`/scratch/$USER/$SLURM_JOB_ID` on Expanse) and read many times; all heavy transient files (LD, GRM, sumstats) are written to and deleted from **node-local NVMe**, never the parallel filesystem, so scaling to hundreds of concurrent tasks does not hammer Lustre. Only the small final result table is copied back to project storage. Steady-state persistent storage is modest (tens of GB); the pressure is on *transient scratch* and on *metadata* (file create/delete rate), not on total capacity.

---

## 5. Code performance & scaling on Expanse-class hardware
All timings are single-core wall-clock on 2× AMD EPYC (Zen3, 128-core) nodes with threads pinned, i.e. the per-core rate you get inside a Slurm allocation on Expanse.

### 5a. The dominant cost driver is the number of SNPs fine-mapped (p), *not* sample size
**→ Figure 1 (`figures/fig1_scaling_p.png`).** For a fixed n = 500, cost rises steeply with the number of SNPs *p* in the *cis*-region. Measured per-model breakdown:

| p (SNPs in region) | Build LD+Cholesky | SuSiE | PLINK LD | GCTA | OLS | **Total / model** | **Peak RAM** | LD file written |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1,544 | 0.2 s | 4.9 s | 0.9 s | 1.4 s | 0.6 s | **8.0 s** | **0.25 GB** | 93 MB |
| 3,567 | 0.9 s | 12.5 s | 4.3 s | 1.4 s | 1.2 s | **20.6 s** | **0.81 GB** | 496 MB |
| 8,227 | 6.6 s | 59.8 s | 23.7 s | 1.8 s | 2.7 s | **96.0 s** | **3.27 GB** | 2.8 GB |
| 17,287 | 45.2 s | 253.0 s | 104.4 s | 3.4 s | 5.6 s | **417.8 s** | **12.9 GB** | 12.3 GB |

Empirical fits: **runtime ∝ p¹·⁷**, **memory ∝ p¹·⁶** (SuSiE's dense p×p operations and the Cholesky dominate; Cholesky alone approaches the theoretical O(p³) at large p). Practical consequence for the allocation: analyses that widen the *cis*-window or move to denser reference panels get expensive **super-linearly**, and per-task memory sets how many tasks fit on a 256 GB node.

### 5b. Sample size (n) barely moves the cost
**→ Figure 2 (`figures/fig2_scaling_n.png`).** At fixed region (p = 3,567), growing n from 100 → 1,000 changes total per-model time only from 19.9 s → 25.0 s. GCTA's REML grows ~O(n²) but is a small fraction; SuSiE and the LD build are **n-independent** (summary-statistic method). This justifies scaling experiments along *p* / number-of-loci rather than along sample size.

### 5c. Scaling from 1 → 1,000 simulations (and why 1,000)
The engine is a **linear loop over independent replicate loci**, so wall-time scales linearly in the number of simulations. Measured ground truth for one full architecture grid (5 ncausal × 6 h² = 30 models) on one region: **3.65 min/region** at p=1,544 (7.3 s/model). Therefore:
- 1 simulated locus (full grid) → ~4 min on 1 core
- **1,000 loci → ~61 core-hours** at p≈1,500, or **~165 core-hours (SU)** at a representative p≈3,000 — a fraction of one node-day (Figure 4).

**Why 1,000 replicates is scientifically necessary:** power, bias, and type-I error are estimated as means/proportions across replicate loci, so their confidence intervals shrink as 1/√N. Monte-Carlo SE of a proportion (e.g. power or type-I error) is ≈ √(p(1−p)/N): **N = 100 → ±5%** CI half-width, **N = 1,000 → ±1.6%**. To claim, e.g., that SuSiE's type-I error is controlled at 5% (not 7%), or to resolve bias differences of a few percent between methods across the h² grid, ~1,000 replicates per cell is the standard in the fine-mapping/TWAS-simulation literature (e.g. twas_sim, the SuSiE papers). Fewer replicates leave the method-comparison conclusions inside the noise.

### 5d. Node-level parallel scaling (1 → many cores)
**→ Figure 3 (`figures/fig3_parallel.png`).** The workload is **embarrassingly parallel**: each locus is an independent Slurm-array task with zero inter-task communication, so scaling *across nodes is linear by construction*. **Within** a single 128-core node we measured throughput as concurrent tasks go 1 → 128:

| Cores | Throughput, light task p=1.5k (models/hr) | Throughput, heavy task p=3.6k (models/hr) |
|---:|---:|---:|
| 1 | 263 | 131 |
| 8 | 1,866 | 532 |
| 32 | 3,029 | 718 |
| 64 | 3,541 | 753 |
| 128 | — | 789 |

Throughput keeps rising but **sub-linearly**, and the gap is explained by **per-task memory footprint**: light 0.25 GB tasks reach ~13.5× speedup at 64 cores, while heavy 0.8 GB tasks saturate near ~6× because the dense p×p linear algebra becomes **memory-bandwidth bound** once many tasks share a socket. This is the key operational finding: the efficient strategy is a *moderate* number of tasks per node (below the bandwidth knee) spread **across many nodes** — which is exactly what a Maximize allocation enables. (Measurements were taken on a shared SDSC node; dedicated Expanse compute nodes will match or exceed the light-task curve.)

### 5e. Memory footprint & other eval metrics
- **Peak RAM** per task = 0.25 GB (p=1.5k) → 12.9 GB (p=17k), scaling ∝ p¹·⁶ (Figure 1, right). At the 2 GB/core Expanse ratio, typical p≈3k tasks (0.8 GB) are core-bound; large-window p≥8k tasks are memory-bound and should be given multiple cores' worth of RAM.
- **Transient disk per task:** the LD matrix write is 0.1–12 GB (§8).
- **Convergence:** the engine retries non-converged SuSiE replicates (up to 5×); a small fraction of high-LD/high-h² cells need retries, which we account for in the budget.

---

## 6. Parameters swept, and their scientific meaning
Each axis isolates one genetic-architecture factor that is known to affect heritability estimation. **These are the variable parameters explicitly considered:**

| Parameter | Levels (measured) | Why it is scientifically meaningful |
|---|---|---|
| True *cis*-h² | 0.05, 0.1, 0.2, 0.4, 0.6, 0.8 (6) | Spans realistic gene-expression heritabilities; tests bias across the full range |
| Number of causal SNPs | 1, 3, 5, 10, 20 (5) | Sparse vs. polygenic *cis* architecture — the regime where fine-mapping should win or fail |
| Sample size *n* | 100, 300, 500, 800, 1,000 (5) | Matches real eQTL study sizes (e.g. GTEx tissues) |
| Effect-size distribution | normal, fixed, heavy-tailed Pareto (3) | Tests robustness to the (unknown) true effect-size prior |
| LD architecture | regular, high-LD, low-LD (3) | LD is *the* confounder fine-mapping is meant to resolve |
| Allele-frequency bin | common vs. rare/MAF bins (2–3) | Rare-variant regime stresses both estimators |
| Replicate loci | up to 1,000 | Confidence intervals on bias/power/type-I error (§5c) |

A single fully-crossed sweep of one replicate = 5×6×5×3×3×2 ≈ **2,700 models per locus**; with replicate loci this is the workhorse campaign. Scientific end-targets: **~20,000 genes × up to 49 GTEx tissues** (real-data validation), extensible to **~100 GWAS traits** for downstream TWAS interpretation, and a planned **third estimator (HEELS)** that multiplies the method axis.

---

## 7. Efficiency work already done, and I/O minimization
- **All heavy transient I/O on node-local NVMe**, not Lustre — LD/GRM/sumstats are written to and deleted from local scratch so thousands of array tasks don't overwhelm the parallel filesystem's metadata server.
- **Aggressive per-model cleanup** (`clean_up()`, `rm` of `.grm.gz/.grm.id/.phen/.ld`) keeps the transient-file count bounded even across a full grid.
- **Cholesky pre-factorization**: the LD matrix is factored once per region and reused to sample all genotype replicates (O(p³) done once, not per replicate).
- **Retry-with-reseed** instead of failing a whole job when SuSiE doesn't converge, so a single bad replicate never wastes a node-hour.
- **Append-only single result table** per campaign (one summary row per model) rather than thousands of small output files.
- **Long-running jobs past the 2-day wall limit:** campaigns are checkpointable by design (independent loci; the result table is resumable — the header/skip logic already appends only new rows), so a campaign is split into Slurm-array chunks and/or auto-requeued, letting a 1,000-locus sweep run past Expanse's 48-hour limit without losing completed work. *(This mirrors the requeue approach the team has used to keep jobs alive past the wall-clock cap.)*

---

## 8. Example bottlenecks (filesystem)
1. **The pairwise-LD text file is O(p²) and huge.** Measured: **93 MB (p=1.5k) → 2.8 GB (p=8k) → 12.3 GB (p=17k)** per model — consistent with the 2.4 GB `temp_ld.ld` still on disk in the repo. Writing/reading multi-GB text LD files is the single biggest I/O cost and must stay on local NVMe.
2. **Redundant LD recomputation.** In the current code the LD file is recomputed for **every** grid cell of a region even though the region's genotypes (and thus its LD) are identical across all 30 cells. Caching LD **once per region** removes ~29/30 of the PLINK-LD I/O and, at p=8k, saves ~23 s **and 2.8 GB of write** per cell — a concrete optimization we will land. *(Flagged here because ACCESS asks for known bottlenecks and planned efficiencies.)*
3. **Metadata churn from many small files.** Each model creates/deletes ~5 small files (GRM, id, phen, sumstats, ld). At hundreds of concurrent tasks this stresses the filesystem's create/delete rate far more than its capacity — another reason all of it lives on node-local NVMe.

---

## 9. Independence from NAIRR
This CPU simulation-and-fine-mapping workload is entirely separate from any NAIRR-pilot (GPU/LLM) usage. It needs standard CPU compute nodes, no accelerators, and the resources described here are requested solely for the *cis*-heritability benchmarking project.

---

## 10. Service-Unit budget (small-experiment × multiplier)
**Method (per ACCESS guidance):** we measured the true single-core cost of *one* model, then multiply by the number of models — no wasteful full-scale scaling runs. On Expanse, **1 SU = 1 core-hour**; a standard node-hour = 128 SU.

Measured single-core cost per model (full pipeline) at representative region sizes:

| p (SNPs) | seconds/model | SU/model (core-hr) |
|---:|---:|---:|
| 1,000 | 3.2 | 0.0009 |
| 3,000 | 19.8 | 0.0055 |
| 5,000 | 46.2 | 0.0128 |
| 8,000 | 100.5 | 0.0279 |

**Per-campaign estimates (Figure 4, `figures/fig4_su_cost.png`, at representative p≈3,000):**

| Campaign | # models | **SU (core-hr)** | node-days |
|---|---:|---:|---:|
| One baseline sim campaign (1,000 loci × 5 ncausal × 6 h²) | 30,000 | **165** | 0.1 |
| Full architecture sweep (1,000 loci × all parameter axes) | 540,000 | **2,971** | 1.0 |
| GTEx real-data pass (49 tissues × 20k genes) | 980,000 | **5,393** | 1.8 |
| 100 GWAS traits (~1,500 loci each) | 150,000 | **825** | 0.3 |

**Annualized request (transparent multipliers the PI can adjust):** method development requires re-running the full sweep many times as the estimator/code evolves; real-data validation runs at *realistic* window sizes (p≈3k–8k, i.e. 1×–5× the per-model cost above) and across multiple estimators (GCTA, SuSiE, + planned HEELS). A defensible annual envelope:

| Activity | Assumption | SU/yr |
|---|---|---:|
| Simulation method development | ~40 full sweeps × ~3,000 SU | ~120,000 |
| Real-data GTEx passes | ~10 passes × ~5,400 SU × ~3 (realistic p, 3 methods) | ~160,000 |
| GWAS-trait fine-mapping | ~5 passes × ~825 SU × ~3 | ~12,000 |
| Replicate/CI + retries + safety margin | ~30% | ~90,000 |
| **Total requested** | | **≈ 380,000 SU/yr** |

Every base number in this table is measured on this hardware; only the *counts* of campaigns are policy choices the PI can tune to the allocation ceiling.

---

### Reproducibility
`bench.py` (per-stage timing + peak RSS), `parallel_bench.py` (node throughput), `fullgrid_bench.py` (ground-truth per-region), `analyze.py` (figures + SU model). Raw measurements in `results/*.json`; figures in `figures/`. Environment: conda env `dsc180` (GCTA 1.94.1, susieR, PLINK 1.90, Python 3 scientific stack) on `mediator.sdsc.edu` (2× AMD EPYC 7713).
