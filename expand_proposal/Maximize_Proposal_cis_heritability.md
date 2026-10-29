X: Fine-mapping-based Estimation of cis-Heritability of Gene Expression

*(Section letter "X" is a placeholder — set it to match your slot in the combined lab document. All timing, memory, storage, and cost figures below are measured on SDSC AMD-EPYC hardware in the lab's `dsc180` environment, not estimated; they are reproducible with `bench.py`, `parallel_bench.py`, `fullgrid_bench.py`, and `analyze.py`/`compute_budget.py` in this folder.)*

# 1. Resource Usage Plan to Achieve the Research Objectives

This project determines when fine-mapping-based estimation of the *cis*-heritability of gene expression is more reliable than the field-standard GCTA/REML estimator. *cis*-heritability — the fraction of a gene's expression variance explained by local genetic variation — is a central quantity for transcriptome-wide association studies (TWAS) and for prioritizing genes and tissues in complex-trait interpretation. We compare two estimators: GCTA/REML, which uses individual-level genotype and expression data, and a **SuSiE fine-mapping-based estimator we developed**, which uses only summary statistics and a linkage-disequilibrium (LD) matrix and summarizes the inferred sparse posterior effects as heritability, `h² = Σ(μ² + σ²)`. Because ground-truth heritability is never observed in real data, the computational core of the project is a **simulation engine** that draws a realistic local LD structure from 1000 Genomes reference genotypes, simulates gene expression from a *known* true *cis*-h² under a controlled grid of genetic architectures, runs both estimators, and measures their error; the simulation-derived conclusions are then confirmed on real GTEx expression data.

Custom software is necessary because no existing tool performs this simulate-and-compare workflow. Generic fine-mapping and heritability tools each cover only one step: SuSiE returns fine-mapping posteriors but not a heritability estimate or its standard error (that derivation is ours), and no public tool simulates expression from a *fixed known cis-h²* while benchmarking GCTA against fine-mapping in a single loop. The repository contains approximately **6,000 lines of homegrown Python, R, and shell** across 32 scripts, excluding reference data and third-party binaries. Internal components provide: real-region LD construction with regularization and Cholesky factorization; multivariate-normal genotype sampling; causal-effect and expression simulation calibrated to a target heritability; a GCTA harness that builds the genetic relationship matrix from simulated genotypes and runs REML; the SuSiE→h² estimator (`run_susie.R`, including allele-flip quality control and the posterior-to-heritability derivation with an analytic standard error); and shell-based grid orchestration with per-replicate retry and temporary-file cleanup. A design constraint we enforce in code is that the LD used to *simulate* genotypes must be the exact same ordered LD handed to the fine-mapper, without which the comparison is uninterpretable. Required third-party software is **GCTA 1.94.1** (REML benchmark), **SuSiE/susieR** (summary-statistic fine-mapping), **PLINK 1.90** (region extraction and pairwise LD), the Python scientific stack (NumPy, SciPy, pandas, statsmodels, scikit-learn, pandas-plink), and the R stack (Rfast, data.table, dplyr, optparse). The workload uses no GPUs; it is CPU- and memory-bandwidth-bound dense linear algebra over p×p LD matrices.

All computation is performed as **CPU-only jobs on Expanse's shared partition**, with no GPU resources and no exclusive whole-node allocations, consistent with the embarrassingly parallel structure of the analysis: each simulated locus, each real gene, and each architecture cell is an independent task with no inter-task communication. The scientific workload is organized as a staged plan rather than a single monolithic sweep:

**(a) Simulation campaigns.** Each campaign draws up to 1,000 replicate loci and, for each locus, sweeps a grid of genetic architectures. The variable parameters, and why each is scientifically meaningful, are:

| Parameter | Levels | Scientific meaning |
|---|---|---|
| True *cis*-h² | 0.05, 0.1, 0.2, 0.4, 0.6, 0.8 (6) | Spans realistic expression heritabilities; tests bias across the range |
| Number of causal SNPs | 1, 3, 5, 10, 20 (5) | Sparse vs. polygenic *cis* architecture — the regime where fine-mapping should win or fail |
| Sample size *n* | 100, 300, 500, 800, 1,000 (5) | Matches real eQTL study sizes (e.g. GTEx tissues) |
| Effect-size distribution | normal, fixed, heavy-tailed Pareto (3) | Robustness to the unknown true effect-size prior |
| LD architecture | regular, high-LD, low-LD (3) | LD is the confounder fine-mapping is meant to resolve |
| Allele-frequency bin | common, rare (2) | Rare-variant regime stresses both estimators |
| Replicate loci | up to 1,000 | Confidence intervals on bias, power, and type-I error |

A fully crossed design is 5 × 6 × 5 × 3 × 3 × 2 = **2,700 architecture cells per locus**. One thousand replicate loci is required because power, bias, and type-I error are estimated as means and proportions across replicates, whose confidence intervals shrink as 1/√N: the Monte-Carlo standard error of a proportion is ≈ √(p(1−p)/N), giving a ±5% half-width at N = 100 but ±1.6% at N = 1,000. Resolving whether the fine-mapping estimator controls type-I error at 5% versus 7%, or resolving few-percent bias differences across the h² grid, requires ~1,000 replicates, consistent with the fine-mapping and TWAS-simulation literature (e.g. twas_sim, the SuSiE papers).

**(b) Real-data analysis.** We apply both estimators (and a planned third summary-statistic method, HEELS) to GTEx expression across up to **49 tissues × ~20,000 genes**, testing whether the simulation-derived reliability patterns hold in real data.

**(c) Downstream and extensions.** Fine-mapping a variable number of loci from GWAS summary statistics for ~100 immune and complex traits; multi-ancestry replication across additional 1000 Genomes reference panels; and, if warranted, extension to larger single-cell expression atlases (e.g. TenK10K-scale) to estimate *cis*-h² across many more cell states.

Job granularity is set empirically rather than assumed (Section 3): individual architecture cells are batched into jobs of a moderate size that keep per-job memory near the shared-partition's 2 GB-per-core billing floor and stay within Expanse's walltime limits, with many such jobs submitted as SLURM arrays that the scheduler spreads across nodes.

# 2. Justification of the Allocation Amounts for All Resources and Resource Types

**Cost model and cross-validation.** The dominant cost driver is the number of SNPs *p* in the *cis*-region, not the sample size. A single-core scaling test across regions of 1,544–17,287 SNPs (n = 500, threads pinned) fit the model **wall_seconds ≈ 3.46 × 10⁻⁵ × p^1.66** per architecture cell (the full GCTA + SuSiE pipeline). We cross-validated this against an independent completed full-grid production run of 30 architecture cells on one region (p = 1,544): the model predicts 197.9 s of compute versus 218.6 s observed, an error of −9.5% (slightly conservative). By contrast, holding the region fixed and varying sample size from n = 100 to 1,000 changes per-cell time only from 19.9 s to 25.0 s, because SuSiE and the LD build are sample-size-independent and only GCTA's REML grows (~O(n²)); we therefore scale experiments along *p* and the number of loci, not along *n*.

Measured per-stage cost at n = 500 (SuSiE fine-mapping and the Cholesky dominate):

| p (SNPs in region) | Build LD + Cholesky | SuSiE | PLINK LD | GCTA | Total / cell | Peak RAM |
|---:|---:|---:|---:|---:|---:|---:|
| 1,544 | 0.2 s | 4.9 s | 0.9 s | 1.4 s | **8.0 s** | 0.25 GB |
| 3,567 | 0.9 s | 12.5 s | 4.3 s | 1.4 s | **20.6 s** | 0.81 GB |
| 8,227 | 6.6 s | 59.8 s | 23.7 s | 1.8 s | **96.0 s** | 3.27 GB |
| 17,287 | 45.2 s | 253.0 s | 104.4 s | 3.4 s | **417.8 s** | 12.9 GB |

**Memory-equivalent billing.** Because Expanse standard compute nodes provide ~2 GB of memory per core, a shared-partition job is billed for the larger of its requested core count or its memory-equivalent core count, `ceil(requested_GB / 2)`. Our peak memory scales as **≈ 1.46 × 10⁻⁶ × p^1.63**, so per-cell memory ranges from 0.25 GB (p ≈ 1.5k) to 12.9 GB (p ≈ 17k). This directly affects effective cost: at a representative *cis*-window (p ≈ 3,000, 0.68 GB) a cell bills at 1 core-equivalent, but a wide 1 Mb+ window (p ≈ 8,000, 3.4 GB) bills at 2 core-equivalents and a very large window (p ≈ 17,000, 11.8 GB) at 6 core-equivalents — up to a sixfold hidden cost multiplier. We right-size requests to ≤ 2 GB per core wherever the region permits, and the estimates below assume this correction. Combined per-cell cost in Expanse core-hours (memory-equivalent billing applied):

| p (SNPs) | wall s/cell | peak GB | core-equiv billed | core-hours / cell |
|---:|---:|---:|---:|---:|
| 1,000 | 3.2 | 0.11 | 1 | 0.0009 |
| 3,000 | 19.8 | 0.68 | 1 | 0.0055 |
| 5,000 | 46.2 | 1.56 | 1 | 0.0128 |
| 8,000 | 100.5 | 3.37 | 2 | 0.0558 |
| 17,287 | 360.0 | 11.82 | 6 | 0.6000 |

**Campaign and annual estimates.** At a representative p ≈ 3,000, one architecture cell costs 0.0055 core-hours, so one *full* simulation sweep (1,000 loci × 2,700 cells = 2.7 million cells) costs ≈ **14,900 core-hours**, and a targeted single-axis campaign (1,000 loci × 30 cells) costs ≈ **165 core-hours**. One GTEx real-data pass for a single estimator (49 tissues × 20,000 genes = 980,000 gene models) costs ≈ **5,400 core-hours**. Aggregating the staged plan:

| Workstream | Models | Rep. p (core-eq) | Core-hours |
|---|---:|---:|---:|
| Simulation: full architecture sweeps (method development, ~6 sweeps) | 16,200,000 | 3,000 (1) | 89,100 |
| Simulation: targeted single-axis campaigns (~40) | 1,200,000 | 3,000 (1) | 6,600 |
| Simulation: large-window / high-memory stress tests (~5) | 150,000 | 8,227 (2) | 8,800 |
| Real data: GTEx *cis*-h² — GCTA | 980,000 | 3,000 (1) | 5,400 |
| Real data: GTEx *cis*-h² — SuSiE | 980,000 | 3,000 (1) | 5,400 |
| Real data: GTEx *cis*-h² — HEELS (3rd method) | 980,000 | 3,000 (1) | 5,400 |
| Multi-ancestry replication (5 reference panels) | 250,000 | 3,000 (1) | 1,400 |
| GWAS-trait locus fine-mapping (100 traits × ~1,500 loci) | 150,000 | 3,000 (1) | 800 |
| Large scRNA-seq atlas extension (e.g. TenK10K-scale) | 600,000 | 3,000 (1) | 3,300 |
| Convergence retries (+10%) | — | — | 12,600 |
| Safety margin (+20%) | — | — | 25,200 |
| **Total requested** | | | **≈ 164,000 core-hours** |

The simulation sweeps dominate because the full architecture grid is combinatorial and each configuration must be replicated ~1,000× for valid confidence intervals; the real-data passes scale with the ~980,000 gene-tissue models per estimator. This request is one component of the lab's shared ACCESS allocation and is not intended to independently approach the lab-wide target; the campaign *counts* are tunable, while the per-model costs are measured.

**Storage.** Persistent storage is modest. The read-only 1000 Genomes EUR reference (chr 1–22) is **14 GB** in 22 large PLINK filesets and is reused by every job. Final results are append-only summary tables of one row per model (≈ 0.1–0.3 MB per campaign). The pressure is instead on **transient scratch and filesystem metadata**: for each model the pipeline writes and deletes a pairwise-LD text file whose size is O(p²) — measured at **93 MB at p ≈ 1.5k, 496 MB at p ≈ 3.5k, 2.8 GB at p ≈ 8k, and 12.3 GB at p ≈ 17k** (consistent with a 2.4 GB `temp_ld.ld` still present in the working directory) — plus ~4 small GRM, phenotype, and sumstats files. A full simulation sweep therefore creates and deletes on the order of ten million transient files over its lifetime. All of these transient objects are written to and removed from **node-local NVMe scratch**, never the parallel filesystem, so the campaign never persists more than a few tens of GB and does not accumulate a large file count on shared storage.

# 3. Resource Appropriateness

All computation uses **CPU-only jobs on Expanse's shared partition**, with no GPU and no exclusive whole-node allocations. This matches the embarrassingly parallel structure of the workload, in which computation divides naturally across replicate loci, genes, tissues, and architecture cells, each an independent task with no communication. Across nodes the workload therefore scales linearly by construction.

**Job granularity was tuned empirically rather than assumed.** We ran a controlled node-level concurrency test, launching an increasing number of identical single-core model-tasks (1, 2, 4, 8, 16, 32, 64, 128) pinned to distinct cores on one 128-core node. Aggregate throughput rose but sub-linearly, and — importantly — the efficiency depended on per-task memory footprint: light 0.25 GB tasks (p ≈ 1.5k) reached ~13.5× speedup at 64 concurrent tasks, whereas heavy 0.8 GB tasks (p ≈ 3.5k) saturated near ~6× at 64 tasks and ~6× at 128, because the dense p×p linear algebra becomes **memory-bandwidth bound** once many tasks share a socket. (These measurements were taken on a shared node and are a lower bound; dedicated compute nodes track the better, light-task curve.) The operational conclusion is that packing a node to its full core count is counter-productive for large-region tasks; instead we batch a moderate number of loci per job and submit many independent array jobs that the scheduler spreads across nodes, capturing the linear cross-node scaling while staying below the per-node memory-bandwidth knee.

**Memory requests are right-sized to control effective cost.** Under the shared partition's memory-driven billing (`ceil(GB / 2)` core-equivalents), an over-requested job silently multiplies its credit cost; our measured per-task footprints (0.25–3.4 GB for typical windows) let us cap most jobs at 2 GB per core so they bill at a single core-equivalent, and we reserve multi-core-equivalent billing only for the deliberately wide-window stress tests where p ≳ 8,000 forces it. Where regions are small, several genes are batched per 2 GB core so that the memory floor is not wasted.

**Filesystem strain is a known risk and is managed.** The workload produces both large objects (multi-GB O(p²) LD matrices) and many small files (per-model GRM, phenotype, and summary-statistic files), a pattern known to strain Expanse's Lustre metadata server. We mitigate this by (i) keeping all heavy and high-count transient I/O on node-local NVMe rather than the parallel filesystem; (ii) aggressive per-model cleanup that bounds the live transient-file count; (iii) an append-only single result table per campaign rather than thousands of small output files; and (iv) a planned consolidation of the redundant LD computation — the LD matrix is currently recomputed for every architecture cell of a region even though a region's genotypes, and thus its LD, are identical across cells, so caching LD once per region eliminates ~29/30 of that I/O and, at p ≈ 8k, saves ~24 s and 2.8 GB of write per cell. Additional efficiency measures already in place include factoring the Cholesky once per region and reusing it across genotype replicates (paying the O(p³) cost once, not per replicate), and retry-with-reseed on SuSiE non-convergence so a single non-converging replicate never wastes a job.

Because loci are independent and the result table skips already-completed rows, campaigns are checkpointable and resumable: a long sweep is split into requeued SLURM-array chunks that survive Expanse's 48-hour walltime limit without losing completed work. The main future scaling dimensions are the number of SNPs per region (the super-linear p^1.66 cost driver, controlled by *cis*-window width and reference-panel density), the number of replicate loci and architecture cells, the number of estimators benchmarked (GCTA, SuSiE, HEELS), and the number of tissues, traits, and ancestries analyzed. This CPU-only simulation-and-fine-mapping workload uses no accelerators and is entirely independent of any NAIRR-pilot usage.
