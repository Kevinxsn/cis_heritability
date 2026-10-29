# cish2 — benchmarking estimators of *cis*-heritability

> For the project as a whole — the scientific question, what was found, design
> decisions and known traps — start at [`../PROJECT.md`](../PROJECT.md).
> This file documents the pipeline itself.

Three estimators of *cis*-heritability, compared on simulated gene expression
where the truth is known and on real GTEx data where it is not:

| estimator | data it needs | what it assumes |
|---|---|---|
| **GCTA / REML** | individual-level genotypes | infinitesimal mixed model |
| **SuSiE** | summary statistics + LD | a **sparse** set of causal effects, `h² = Σ(μ² + σ²)` |
| **HEELS** | summary statistics + LD | infinitesimal mixed model |

HEELS shares SuSiE's inputs and GCTA's model, which is what lets a
SuSiE-versus-GCTA difference be attributed to the sparsity assumption rather
than to working from summary statistics.

Genotypes come from the **1000 Genomes EUR** reference panel; the REML
benchmark is **GCTA 1.94.1**; HEELS is implemented here from Li, Mazumder &
Lin, *Nat Commun* 14:7954 (2023). Everything runs in the `dsc180` conda
environment.

```
cish2/                      code   (this directory)
results/<campaign>/         output (one directory per campaign)
```

---

## Quick start

```bash
cd /projects/ps-renlab2/sux002/cis_heritability/cish2
source scripts/env.sh

python -m cish2 plan  configs/baseline.yaml     # what would this run?
python -m cish2 run   configs/baseline.yaml     # run it (resumable)
python -m cish2 analyze ../results/baseline_n300  # tables + figures + REPORT.md
```

or in one step, `scripts/run_campaign.sh configs/baseline.yaml`.

```bash
scripts/analyze_all.sh                        # every campaign + overview + write-up
python -m analysis.overview ../results        # cross-condition tables and figures
python -m analysis.findings ../results        # the narrative FINDINGS.md
python -m analysis.compare_legacy \
    --new ../results/baseline_n300 \
    --legacy ../simulations/result_eur_num_people_100_300_regular_LD.txt
```

Read the results starting at **`results/FINDINGS.md`** (simulations, narrative
with figures), then **`results/gtex_whole_blood/GTEX_FINDINGS.md`** (whether
those patterns hold in real GTEx data), then `results/OVERVIEW.md` (all
simulated conditions in one table), then a campaign's own `REPORT.md`.

### HEELS implementation notes

The published algorithm alternates two variance-component updates to a fixed
point. We solve for that fixed point directly instead: heritability depends on
the components only through their ratio `a = s2e/s2g` (`h² = 1/(1+a)`), and for
fixed `a` the residual-variance update *is* the profile maximiser of the scale,
so one bracketed 1-D search over `log a` lands on the same answer. This is not
an approximation — it agrees with the iteration to ~1e-8 — and it matters,
because the iteration crawls near the `h² → 0` boundary (one locus still read
0.098 after 200 steps against a converged 0.067).

The eigenbasis comes from an SVD of the `n × p` genotype matrix rather than an
eigendecomposition of the `p × p` LD: same decomposition, `O(n²p)` instead of
`O(p³)`, and it is the in-sample LD the method requires. When `p > n` the
`p − n` null-space directions carry no signal but do contribute `(p−rank)/a` to
`tr(W⁻¹)`, and are handled exactly.

`HeelsSolver.solve_em` keeps the published iteration for validation; a test
asserts the two agree.

### Real data

```bash
python -m cish2 gtex configs/gtex_whole_blood.yaml   # ~17 min on 48 workers
python -m analysis.gtex_report                        # prediction tests + report
```

The GTEx pass fits all three estimators for every gene at five **nested**
sample-size subsets, turning each simulation result into a falsifiable
prediction. Real data has no ground truth, so it measures the *disagreement*
`estimator − GCTA`, not bias. Here both summary-statistic estimators get
correctly specified in-sample LD, so SuSiE-versus-HEELS is a clean test of the
sparsity assumption in a way the simulation comparison is not.
Covariates are handled per estimator, each the way it expects: **REML gets them
as fixed effects** (`--qcovar`) on the unadjusted phenotype, while the
summary-statistic estimators work on residualised expression and genotypes with
`n_eff = n − k`. The previous GTEx script accepted a `--covariate_file` and
never used it at all; an earlier version of *this* pipeline residualised for
everyone, which hid the lost degrees of freedom from REML and inflated it badly
at small *n*. The inferred-factor count follows GTEx's own sample-size rule.

**Semi-synthetic and permutation modes (round 4, default off).** The GTEx
runner takes `phenotype_mode: real | simulated | permuted` (default `real`,
which is every campaign above, unchanged and regression-tested to 1e-8):

- `simulated` keeps the donors, the split, the covariates and every estimator,
  and replaces each gene's expression with `Zβ + ε` on that gene's **real**
  genotypes (`simulate.sim_beta_empirical`, `var(Zβ) = h²` exactly). `(h², k)`
  is drawn per gene and draw from `sim_design`; both halves share one β. Rows
  gain `draw, h2_target, ncausal_target, h2_sample, h2_resid, sim_beta_hash`.
- `permuted` keeps real expression and covariates and permutes genotype rows
  within each subset, so true cis-h² = 0 for every gene.
- `gcta_also_unconstrained: true` adds a second, unconstrained REML fit
  (`gcta_h2_uc`, `gcta_h2_uc_se`, `gcta_uc_pval`, `gcta_uc_status`).

RNG streams are `(seed, crc32(gene_id), draw)`, so results do not depend on
chunking or worker order. Configs: `gtex_sim_split_half.yaml`,
`gtex_perm_split_half.yaml`.

### Conditions, not campaigns

A campaign like `sample_size` sweeps five values of *n* inside one directory.
Reporting it as one row averages the sweep away, so the reader cannot tell
whether "sample size" meant more or less of it. `analysis/conditions.py` splits
every campaign into one **condition per swept level**, each carrying its
direction relative to the baseline (*n* = 300, normal effects, unrestricted
causal MAF, unselected regions). All cross-campaign output is per condition.

One naming caution, encoded in that module: the `effect_dist` axis does **not**
vary effect magnitude. `sim_beta` rescales every draw so `β'Rβ = h²` exactly, so
all three arms carry identical total effect — what changes is how *concentrated*
it is across causal SNPs (`fixed` spread evenly → `normal` → `pareto`
concentrated). Effect magnitude is already the `h2_true` axis (total) and the
`ncausal` axis (per-SNP, `√(h²/k)`).

Everything is driven by a YAML config; there is **one** simulation engine, and
LD architecture / MAF bin / effect-size prior / sample size / causal count /
true h² are all axes of its grid rather than separate scripts.

```bash
python -m cish2 run configs/baseline.yaml -o n_regions=10 -o compute.workers=8
python -m cish2 run configs/baseline.yaml --name baseline_matchedLD \
    -o method.ld_for_susie=sim
```

---

## What each campaign config covers

| Config | Axis under study | Replaces |
|---|---|---|
| `baseline.yaml` | causal count × true h², n = 300 | `simulations_test.py` |
| `sample_size.yaml` | n = 100 … 1000 | `simulations.py` n-sweep runs |
| `effect_size.yaml` | normal / fixed / Pareto effect prior | `simulation_effect_size.py` |
| `maf_arch.yaml` | common vs. rare causal SNPs | `simulation_high_al_freq.py` |
| `high_ld.yaml` | high-LD regions | `simulation_high_ld_test.py` |
| `low_ld.yaml` | low-LD regions | `simulation_low_ld_test.py` |
| `smoke.yaml` | end-to-end check, ~30 s | — |

---

## How the pipeline works

For each **region slot** (the unit of parallelism, one independent task):

1. **Draw a region.** A random variant, ± `window_bp`, taken directly out of
   the chromosome `.bed` by byte offset. High/low-LD campaigns accept-reject on
   a subsampled LD score, so a rejected region costs a slice read, not a p×p
   matrix.
2. **Build LD once.** `R_panel` (raw reference correlation) and `R_sim`
   (ridge-regularised) plus the Cholesky factor of `R_sim` — one GEMM, reused
   by every cell of the grid.
3. **Simulate once per sample size.** Genotypes `Z ~ N(0, R_sim)` via the
   Cholesky factor; the GRM built from `Z` is identical across the grid, so
   it is written once.
4. **Per cell:** draw causal effects rescaled to hit the target h² exactly,
   add noise, run REML, and accumulate marginal sumstats.
5. **HEELS per cell**, against an eigenbasis built once per (region, n) from
   the SVD of the simulated genotypes — 0.015 s per fit.
6. **One SuSiE call per region.** All cells' sumstats go to R in one batch
   against the LD matrix that is already in R's memory.
7. **Write a shard**, append to the manifest, delete the scratch.

### Which LD the fine-mapper sees

`method.ld_for_susie` picks between two scientifically different setups:

- `panel` — the raw reference-panel correlation matrix. This is what the
  legacy PLINK `--r` pipeline handed SuSiE, so it is the setting for
  reproducing historical results. Note it is rank-deficient whenever
  *p* > 503 (the panel size), which is the normal case for a ±500 kb window.
- `sim` — the same regularised matrix the genotypes were sampled from. This is
  the design principle stated in the proposal ("the LD used to simulate must be
  the LD given to the fine-mapper") and removes LD mismatch as a confounder.

**HEELS does not use that setting — it always gets in-sample LD**, because it
has no choice. Its residual-variance update stays non-negative only through
`S' R⁺ S = y' P_X y ≤ y'y`, which holds when R is the LD of the same
individuals the statistics came from. With a reference panel the bound fails
outright (measured: `S' R⁺ S` = 1383 against `y'y` = 300 on a p = 893, n = 300
region) and the iteration diverges. The consequence is an asymmetry worth
remembering: in the `panel` configuration **HEELS sees better-specified LD than
SuSiE does**, so part of the difference between them is the input, not the
model. On the GTEx pass this asymmetry vanishes — both get in-sample LD.

`method.legacy_ambiguous_filter` reproduces the legacy allele QC, which dropped
strand-ambiguous SNPs (A/T, C/G) before fine-mapping. Since our sumstats and
reference panel are the same `.bim`, there is nothing to disambiguate and the
filter only removes usable SNPs — it exists to match historical runs.

---

## Efficiency

The engine is CPU- and memory-bandwidth-bound dense linear algebra. The
optimisations that matter, relative to the previous script family:

| Change | Effect |
|---|---|
| LD passed to R as raw float64, not a PLINK long-format text file | at p = 8k: one 520 MB binary write per **region** instead of a 2.8 GB text write + parse per **cell** |
| One `Rscript` per region instead of per cell | R startup and LD load paid once per 30 cells, not 30 times |
| Region extraction by direct `.bed` byte-range read | no PLINK subprocess, no per-region temp filesets |
| Vectorised marginal OLS (one BLAS call) | replaces p `statsmodels.OLS` fits per cell |
| GRM built once per (region, n), written in GCTA **binary** format | ~15× less CPU and ~4× less disk than gzipped text, and 1× instead of 30× |
| Transient I/O auto-routed to node-local NVMe | the project filesystem here is NFS at ~117 MB/s vs. ~1.5 GB/s local; this alone roughly halved wall time at 48 workers |
| Multiprocessing over regions, 1 BLAS thread per worker | linear scaling below the memory-bandwidth knee |

Measured on `mediator.sdsc.edu` (2× AMD EPYC 7713): **≈ 18 s of single-core
work per grid cell** at a median region size of ~8,200 SNPs, so a 100-region ×
30-cell baseline campaign is ~15 core-hours and finishes in ~23 minutes on 48
workers. Cost scales as p^2.0 in the SNPs per window and is essentially flat in
sample size.

**Sizing note.** Peak memory per worker is dominated by three p×p float64
matrices plus R's copy — about 2.5 GB at p ≈ 8k and ~4× that at p ≈ 16k. Set
`compute.workers` from available RAM, not just core count, and cap
`data.max_snps` if regions get very wide.

---

## Results layout

```
results/<campaign>/
  config.yaml          exact config this run used
  manifest.jsonl       one line per finished region slot -> resumability
  regions.tsv          one row per drawn region
  shards/*.tsv.gz      per-slot results, written as each slot completes
  cells.tsv.gz         consolidated long table, one row per grid cell
  tables/*.csv         metrics with Monte-Carlo standard errors
  figures/*.png|pdf    400 dpi review copies + vector masters
  REPORT.md            headline numbers, attrition, cost, figures
  logs/run.log

results/                 (across campaigns)
  FINDINGS.md          narrative walk-through with figures - start here
  OVERVIEW.md          every condition in one table
  conditions.csv       per-condition metrics (machine readable)
  overview.csv         per-campaign provenance
  figures/             fig8 condition grid, fig9 effect-of-each-change
```

`cells.tsv.gz` carries region identity, the design cell, the realised
sample heritability, **both** estimators' point estimates, standard errors and
diagnostics (SuSiE credible sets, max PIP, iterations; GCTA's V(G), V(e), Vp,
p-value), the convergence status of each, and per-stage timings.

A campaign is **resumable and idempotent**: rerunning it computes only the
region slots missing from `manifest.jsonl`, so an interrupted or extended run
never duplicates or loses work.

---

## Figures

Two families: the **metric** figures (means with Monte-Carlo intervals, which
answer "how biased?") and the **distribution** figures (`figD*`, box-and-whisker
over replicate loci, which answer "how spread out, and how skewed?"). Means and
spread tell different stories here — fine-mapping is tighter than REML in the
sparse corner and wider in the polygenic one — so both are generated for every
campaign.

Each is written as a 400 dpi PNG and a vector PDF (fonts embedded, editable in
Illustrator). Colors come from a validated colorblind-safe palette: one fixed
hue per estimator, a single-hue ramp for ordinal encodings, and a two-pole
diverging map with a neutral midpoint for "who wins".

| Figure | Question it answers |
|---|---|
| `fig1_accuracy_grid` | Does each estimator track the truth, and where does it depart? |
| `fig2_bias_by_polygenicity` | How does bias respond to the number of causal SNPs? |
| `fig3_rmse_landscape` | RMSE per estimator over the (causal count × h²) grid, plus which one wins each cell |
| `fig4_calibration` | Are the reported standard errors honest, and do the 95 % intervals cover? |
| `fig5_head_to_head` | Each summary-statistic estimator against GCTA: agreement, Bland-Altman, win shares |
| `fig6_axis_<axis>` | What does the campaign's extra axis (n / prior / MAF bin) change? |
| `fig7_compute_profile` | How does cost scale with region size, and where does the time go? |
| `figD1_estimate_distributions` | Box-and-whisker of every estimate against the truth trace |
| `figD2_error_distributions` | The same, centred on zero — spread as well as bias |
| `figD3_relative_error` | Percent error, the scale that matters at low true h² |
| `figD4_boxplots_<axis>` | Error distributions across a campaign's extra axis |
| `gtex_figD1_distributions` | Real-data spread, and the paired per-gene disagreement |
| `fig8_condition_grid` | Every swept level side by side, direction labelled (overview) |
| `fig9_condition_effects` | What each change does relative to the baseline (overview) |
| `fig0_legacy_reproduction` | Does the rebuilt engine match the old pipeline? |

---

## Metrics

All metrics are means or proportions over replicate loci and are reported with
their Monte-Carlo standard error (`sd/√N` for a mean, `√(p(1−p)/N)` for a
proportion) — that is what makes a few-percent bias difference claimable rather
than eyeballed.

`bias`, `rel_bias`, `rmse`, `mae`, `median_abs_err`, `sd_h2_hat`, `coverage95`
(95 % interval coverage of the truth), `corr_truth`, `power05` (estimate
significantly > 0), plus a paired `susie_win_rate` and a Wilcoxon signed-rank
test on absolute error — paired because both estimators see the *same*
simulated phenotype.

Two robustness diagnostics are reported beside them, because the two estimators
fail in opposite ways. `frac_impossible` counts estimates outside [0, 1]:
GCTA under `--reml-no-constrain` produces small negatives *often* (22 % at
n = 100), while SuSiE on a rank-deficient panel LD produces runaway estimates
*rarely* (<1 %) but 3-4x out of range, reporting clean convergence as it does
so. `max_abs_err` exposes those. **Prefer MAE to RMSE** in any summary that
spans them — a single runaway fit can supply 70 % of a squared-error total.

---

## Reproduction check

`analysis/compare_legacy.py` puts a legacy result table and a new campaign on
the same footing and compares *cell means over replicate loci* (the two
pipelines draw different random regions, so agreement is expected within
Monte-Carlo error, not row by row). Against
`simulations/result_eur_num_people_100_300_regular_LD.txt`, the rebuilt engine
gives cell-mean correlations of **0.999** for both estimators, with 93 % (GCTA)
and 100 % (SuSiE) of the 30 design cells inside 2 Monte-Carlo standard errors —
i.e. indistinguishable from sampling noise.

---

## Tests

```bash
source scripts/env.sh && python tests/test_pipeline.py
```

Covers the parts that fail silently: the `.bed` reader against
`pandas_plink`, vectorised OLS against `statsmodels`, effect rescaling hitting
the target genetic variance, sampled genotypes reproducing their LD, the
GCTA binary-GRM layout, `.hsq` parsing, and the Python→R matrix round-trip.

---

## Notes on faithfulness to the previous code

Deliberate differences, all documented in `cish2/config.py`:

- **Causal SNPs are drawn without replacement.** The original baseline script
  used `np.random.choice(..., ncausal)` (with replacement), so a "5-causal"
  model occasionally had fewer than five distinct causal SNPs. The later
  variant scripts had already fixed this.
- **GCTA h² is reported as V(G)/Vp** in `gcta_h2`, with the raw V(G) kept in
  `gcta_vg`. The legacy tables recorded V(G) and called it h²; the two agree
  closely only because the phenotype is standardised.
- **Failures are recorded, not silently coerced.** A REML failure is `NaN` with
  its status string rather than `(0, 10000, 1)`, and the attrition table counts
  them.
- **The high/low-LD thresholds were re-derived from measurement.** The legacy
  scripts accepted a region as "high LD" when the fraction of sampled SNP pairs
  with r² ≥ 0.8 exceeded 0.02. Measured over 300 random ±500 kb 1000G EUR
  windows that statistic has median 0.0024 and a maximum of 0.0109 — the cut is
  unreachable and rejects every region. `high_ld.yaml` and `low_ld.yaml`
  instead use the empirical **top and bottom deciles** (≥ 0.005 and ≤ 0.0014),
  which is the contrast the axis was meant to create.
- Not yet ported: the HEELS estimator and the GTEx real-data pass.
