# *cis*-heritability of gene expression — project guide

**Question.** When is a **SuSiE fine-mapping** estimator of *cis*-h² more
reliable than the field-standard **GCTA/REML**, and how much of any difference
is the *data* each needs versus the *model* each assumes?

| estimator | data | assumption |
|---|---|---|
| **GCTA / REML** | individual-level genotypes | infinitesimal mixed model |
| **SuSiE** | summary statistics + LD | a **sparse** set of causal effects, `ĥ² = Σ(μ² + σ²)` |
| **HEELS** | summary statistics + LD | infinitesimal mixed model |

The third arm is what makes the first two interpretable. HEELS (Li, Mazumder &
Lin, *Nat Commun* 2023) shares SuSiE's inputs and GCTA's model, so a
SuSiE-versus-GCTA difference can finally be attributed to one or the other.

Ground truth is never observable in real data, so the work has two halves that
answer different things:

| half | what it can measure | where |
|---|---|---|
| **Simulation** | **bias** — the truth is set by construction | [`results/FINDINGS.md`](results/FINDINGS.md) |
| **Real data (GTEx)** | **disagreement**, and whether it behaves as predicted | [`results/gtex_whole_blood/GTEX_FINDINGS.md`](results/gtex_whole_blood/GTEX_FINDINGS.md) |
| **Real data, reliability** | **precision** — test-retest on disjoint donor halves | [`results/gtex_split_half/RELIABILITY.md`](results/gtex_split_half/RELIABILITY.md) |

The simulation half produces predictions; the real-data half tests them. Neither
alone is worth much — that pairing is the design of the project.

---

## Where everything lives

```
cis_heritability/
  PROJECT.md              this file — start here
  cish2/                  all code (~6,400 lines)
    README.md             pipeline documentation, how-to, efficiency notes
    cish2/                simulation + GTEx engine
    analysis/             metrics, figures, distributions, report generators
    configs/*.yaml        one file per campaign
    R/                    susie_batch.R (simulation), susie_genes.R (GTEx)
    scripts/              env.sh, run_campaign.sh, analyze_all.sh
    tests/                self-checks for the things that fail silently
  results/
    FINDINGS.md           simulation write-up  ← read first
    OVERVIEW.md           all 13 simulated conditions in one table
    conditions.csv        the same, machine readable
    figures/              cross-condition figures (fig8, fig9)
    <campaign>/           per-campaign REPORT.md, tables/, figures/, cells.tsv.gz
    gtex_whole_blood/     GTEX_FINDINGS.md + the real-data results
    gtex_split_half/      RELIABILITY.md — test-retest on disjoint donor halves
  poster/                 IGVF 2026 poster (36x48 in) + its own figure set
    README.md             how to rebuild it, and the design decisions
    XuS_IGVF2026_poster.pptx      the deliverable
    make_figures.py       poster-scale redraws of the key panels
    build_poster.py       content + layout engine + PNG proof
    check_numbers.py      re-derives every quoted number from the tables
  simulations/            SUPERSEDED (see note below)
  simulation_gtex/        SUPERSEDED
  expand_proposal/        ACCESS allocation proposal + its benchmarks
```

**Data** (read-only, not in this tree):
`/projects/ps-renlab2/sux002/data/` — 1000G EUR reference panel (chr 1–22) for
the simulations, and `data/gtex/` — GTEx v8 EUR HapMap3 genotypes, Whole Blood
expression, covariates.

### On the superseded folders

`simulations/` and `simulation_gtex/` hold the original script family (a
separate `simulation_*.py` per axis). Everything they did is now one
config-driven engine, and they can be ignored for new work. **Do not delete
`simulations/` outright**: `analysis/compare_legacy.py` reads
`simulations/result_eur_num_people_100_300_regular_LD.txt` for the reproduction
check documented below. Nothing else in either folder is referenced.

---

## What we found

### Simulations — 20,294 model fits, 350 regions, 13 conditions

**Polygenicity is the switch.** With one causal SNP the fine-mapping estimator
is essentially unbiased (−0.008) and closer to the truth on 76% of loci. At 20
causal SNPs it loses 0.13 of h² and wins 27%. GCTA is flat across the range.
The mechanism is structural: SuSiE only counts variance it can localise, and
effect spread thinly across variants in LD never enters the sum.

**Sample size is the strongest modifier — and it runs opposite to usefulness.**
Bias goes −0.122 (n=100) → −0.021 (n=1,000), monotone. But at n=100, where
SuSiE is *most* biased, its MAE is 31% *lower* than GCTA's, because REML is so
noisy there. By n=1,000 the edge is 11%. Fine-mapping is most worth using
exactly where it is most biased.

**Concentration matters as much as count.** At 20 causal SNPs the bias is
−0.222 with equal-magnitude effects but −0.069 with heavy-tailed effects. Note
the effect-size axis does *not* vary magnitude — `sim_beta` rescales so
`β'Rβ = h²` exactly — it varies how concentrated a fixed total effect is.

**GCTA breaks in exactly one place.** Restricting causal variants to a MAF bin
turns it **upward**-biased (+0.045 for common variants), the only sign flip in
the study — the known GREML MAF-stratification issue. So "SuSiE differs from
GCTA on real data" cannot be read as "SuSiE is wrong".

**LD architecture matters far less than expected.** Top vs bottom decile
separates the bias by only 0.009 — an order of magnitude less than the causal
count. What defeats fine-mapping is signal being *diffuse*, not *correlated*.

**HEELS answers the question the other two could not.** It shares SuSiE's
inputs and GCTA's model, and it lands on GCTA: bias +0.003 against −0.006,
coverage 89% against 87%. So the polygenicity failure is **not** the price of
working from summary statistics — it is the price of assuming sparsity. This
near-identity is expected rather than remarkable: `(X'y, X'X)` are sufficient
statistics for the individual-level likelihood. The useful corollary is that
**nothing is lost by sharing summary statistics instead of genotypes, provided
the LD is in-sample.** HEELS also costs 0.015 s per fit, ~800× less than SuSiE.

**Opposite failure profiles.** GCTA leaves [0,1] often but mildly (22% at
n=100, small negatives from `--reml-no-constrain`); SuSiE rarely (<1%) but
catastrophically (up to 4.3), while reporting clean convergence; HEELS never,
by construction. One runaway fit supplied 70% of a squared-error total —
**report MAE beside RMSE.**

**Fine-mapping trades variance for bias — and at high polygenicity loses both.**
At one causal SNP SuSiE's inter-quartile error spread is 0.056 against GCTA's
0.158, which is why it wins so often there. By 20 causal SNPs its box has both
slid down *and* widened to 0.149, overtaking GCTA's 0.120. The `figD*`
box-and-whisker figures show this directly; the mean-based figures cannot.

**Neither estimator's error bars are trustworthy.** 95% intervals cover 87%
(GCTA) and 68% (SuSiE). SuSiE's analytic SE quantifies uncertainty about the
effects it *found*, not the effects it missed — and the missed effects are the
whole problem.

### GTEx Whole Blood — ~19,200 genes × 5 nested sample sizes = 97,160 fits

Real data has no truth, so this measures the disagreement `estimator − GCTA`
and asks whether it behaves as predicted. The one tissue available can't give a
cross-tissue sample-size axis, so it comes from **nested subsampling**
(n=100 ⊂ 200 ⊂ 300 ⊂ 450 ⊂ 574), which holds the biology fixed.

| prediction | verdict |
|---|---|
| gap closes as n grows | **direction only** — −0.015 → −0.007, i.e. ~6× smaller than simulated |
| gap is largest where signal is diffuse | **confirmed after conditioning on h²** — see below |
| REML fails often/mildly, SuSiE rarely/wildly | **not testable** — constrained REML cannot go negative |
| HEELS tracks GCTA where SuSiE need not | **confirmed emphatically** — correlation **0.997–1.000** over ~19,000 genes |

**An earlier version of this analysis was wrong and is superseded.** It
pre-residualised expression on covariates and handed REML the residual, hiding
the degrees of freedom the covariates consume. That inflated GCTA at small *n*
and manufactured a large gap that appeared to close with sample size. REML now
receives covariates as fixed effects (`--qcovar`); the gap is ~10× smaller and
nearly flat. Doing it properly forced a second choice: unconstrained REML then
fails on 60–75% of genes, so the GTEx pass uses constrained REML (~2%
attrition), which is why GCTA can no longer produce a negative estimate.

The residual story is that **the simulated failure mode is real but rarely
triggered in Whole Blood** — most *cis*-windows here have zero or one detectable
independent signal, exactly the corner where the simulations say fine-mapping is
fine.

**The architecture test nearly came out backwards.** Stratified by LD-clumped
independent signals, the *marginal* table says the gap grows with signal count.
That is confounded — genes with more detectable signals also have far more
heritability (0.02 → 0.37). Conditioning on h² reverses it: within *every*
stratum the gap **shrinks** as signals accumulate. Textbook Simpson's paradox,
and the conditional direction is the one matching the simulations. Clump count
measures how **localisable** a gene's signal is, not how polygenic it is; the
worst case is a gene with real heritability and no detectable eQTL. Both views
are in the figure so the confounding is on the record.

### Practical summary

| situation | use | why |
|---|---|---|
| no individual-level data | **HEELS** | matches REML exactly (corr 1.000 on real genes) — if LD is in-sample |
| summary stats + *reference-panel* LD | **SuSiE** | HEELS diverges outright on mismatched LD when p > n |
| small eQTL study (n ≲ 300) | **SuSiE**, expect downward bias | REML variance dominates; ~30% lower MAE |
| large study (n ≳ 800) | **either**, GCTA if unbiasedness matters | they converge; SuSiE keeps a small residual bias |
| sparse architecture | **SuSiE** | 76% win rate at one causal variant |
| polygenic architecture | **GCTA or HEELS** | both flat across polygenicity; SuSiE loses ~0.13 of h² at 20 causal SNPs |
| frequency-stratified | **neither unmodified** | GCTA turns upward-biased; use a MAF-stratified GRM |
| gene with h² but no detectable eQTL | **distrust fine-mapping h²** | the largest disagreement in real data |
| any reported h² standard error | **do not trust as-is** | coverage is 68–87%, not 95% |
| any fine-mapping h² in production | **screen against [0,1]** | runaway fits report clean convergence |

---

## Running it

```bash
cd /projects/ps-renlab2/sux002/cis_heritability/cish2
source scripts/env.sh                      # dsc180 env + BLAS pinned to 1 thread

python -m cish2 plan  configs/baseline.yaml   # what would this run?
python -m cish2 run   configs/baseline.yaml   # ~23 min on 48 workers, resumable
python -m cish2 gtex  configs/gtex_whole_blood.yaml   # ~17 min

scripts/analyze_all.sh                     # every campaign + OVERVIEW + FINDINGS
python -m analysis.gtex_report             # GTEx prediction tests
python tests/test_pipeline.py              # self-checks
```

Any campaign is **resumable and idempotent** — rerunning computes only the
region slots missing from `manifest.jsonl`. Override anything from the command
line: `-o n_regions=10 -o compute.workers=8`.

**Cost.** ~18 s of single-core work per simulated model fit at a median region
of ~8,200 SNPs, scaling as p^2.0 and essentially flat in sample size. The whole
simulation study is ~101 core-hours (~3.5 h wall on 48 workers); the GTEx pass
is ~13 core-hours (17 min).

---

## Design decisions worth knowing

**One engine, many axes.** LD architecture, MAF bin, effect-size prior, sample
size, causal count and true h² are axes of one grid, not separate scripts.
Adding an axis is a config change.

**Conditions, not campaigns.** A campaign that sweeps five values of *n* lives
in one directory; reporting it as one row would average the sweep away and hide
the direction. `analysis/conditions.py` splits every campaign into one
condition per swept level, each labelled with its direction against the
baseline. All cross-campaign output is per condition.

**HEELS is given in-sample LD, always.** Not a preference — a requirement. Its
residual-variance update stays non-negative only through
`S' R⁺ S = y' P_X y ≤ y'y`, which holds when R is the LD of the same people the
statistics came from. With a reference panel the bound fails outright (measured:
`S' R⁺ S` = 1383 against `y'y` = 300 at p = 893, n = 300) and the iteration
diverges in a few steps. Consequence: in the `panel` configuration **HEELS sees
better-specified LD than SuSiE**, so part of the difference between the two
summary-statistic methods is input rather than model. On the GTEx pass the
asymmetry vanishes — both get in-sample LD — which is why that comparison is the
clean one.

**Which LD the fine-mapper sees** (`method.ld_for_susie`) is a real scientific
choice, not plumbing:
- `panel` — raw reference-panel correlation, what the legacy PLINK pipeline
  used. Rank-deficient whenever p > 503, which is normal for a ±500 kb window.
  **All simulation results reported so far use this.**
- `sim` — the same regularised LD the genotypes were sampled from. Removes LD
  mismatch as a confounder. **Not yet run.**

**`legacy_ambiguous_filter`** reproduces the old allele QC, which dropped
strand-ambiguous SNPs before fine-mapping. Since sumstats and reference come
from the same `.bim`, there is nothing to disambiguate — it is a ~15% unearned
handicap to SuSiE, kept on only to match historical runs.

**Efficiency choices that mattered** (details in `cish2/README.md`): LD crosses
into R as raw float64 rather than a PLINK text file; one `Rscript` per region
rather than per cell; direct `.bed` byte-range reads instead of PLINK
subprocesses; vectorised marginal OLS; the GRM built once per (region, n) in
binary format; and **transient I/O routed to node-local NVMe** — the project
filesystem is NFS at ~117 MB/s against ~1.5 GB/s local, and that one change
roughly halved wall time at 48 workers.

**Metrics carry Monte-Carlo standard errors** (`sd/√N` for means,
`√(p(1−p)/N)` for proportions). That is what makes a few-percent bias
difference claimable rather than eyeballed.

---

## Traps

Things that cost real debugging time and would silently corrupt results if
reintroduced. All are now covered by `tests/test_pipeline.py` or by explicit
comments at the site.

1. **`numpy.ndarray.tofile` always writes C order**, ignoring the array's own
   memory layout. Matrices handed to R (column-major) must be transposed before
   writing — marking them Fortran-ordered does nothing. This silently garbled
   every SuSiE input until caught.
2. **`susie_rss`'s default `check_prior = TRUE` is load-bearing.** With
   rank-deficient panel LD the IBSS fit runs away and returns h² in the tens of
   thousands. Never pass `check_prior = FALSE` here.
3. **Residualising propagates NaN across a whole column.** `X − Q(Q'X)` turns a
   genotype column holding one missing call into an all-NaN column, silently
   dropping that SNP. Impute *before* residualising.
4. **Covariate count must scale with sample size.** At n=100, 68 covariates
   plus an intercept consume 69 of 100 df and REML fails on every gene. The
   GTEx pass follows GTEx's own rule (15/30/45/60 inferred factors by n).
5. **RMSE is fragile here.** A single runaway fine-mapping fit supplied 70% of
   a squared-error total, creating a fake reversal in the sample-size trend.
   Always read MAE beside it; `frac_impossible` and `max_abs_err` are in every
   metrics table.
6. **Marginal stratification on real data is confounded.** The architecture
   test reverses direction unless conditioned on heritability.
7. **HEELS diverges on reference-panel LD when p > n**, silently and fast —
   `s2e` goes negative, `a → 0`, and the iteration blows up. Always feed it
   in-sample LD. Its published setting is n ≫ p, where the problem cannot arise;
   ours is p/n ≈ 20–30, where it always does.
8. **HEELS's published iteration crawls near the h² → 0 boundary.** One locus
   still read 0.098 after 200 steps against a converged 0.067. We maximise the
   profile likelihood directly instead — same fixed point, verified to ~1e-8,
   a few hundred evaluations instead of 2×10⁵ iterations.
9. **Pre-residualising a phenotype hides degrees of freedom from whatever
   consumes it next.** Residualising on k covariates leaves n−k, and both REML
   and HEELS inflate badly if told n. Give REML the covariates via `--qcovar`;
   give HEELS `n_eff = n − k`. This one silently corrupted a headline GTEx
   result before it was caught.
10. **GCTA needs a full-rank covariate block.** GTEx's ships rank-deficient (a
   constant column), and REML then fails with "the X^t * V^-1 * X matrix is not
   invertible" on every gene.
11. **`np.linalg.svd` fails to converge on ~12% of real genotype chunks.** The
   default LAPACK driver (`gesdd`) is not robust here; fall back to `gesvd` and
   then to an eigendecomposition of the Gram matrix.
12. **The high-LD threshold from the legacy scripts is unreachable.** It
   required a pair-fraction ≥ 0.02; over 300 random ±500 kb 1000G windows the
   statistic has median 0.0024 and maximum 0.0109. The configs now use measured
   top/bottom deciles.
13. **The SuSiE h² is the TRACE form; the exact LD-aware form needs BOTH cross
   terms.** `susie_batch.R` reports `sum(mu^2 + sd^2)` = `tr(E[bb'])`, but the
   simulated truth is `beta' R beta`. SuSiE's posterior is **not** diagonal:
   components are independent across `l`, but within one single effect the SNPs
   are mutually exclusive. With `m_l = alpha_l * mu_l` and
   `s_l = sum_j alpha_lj * mu2_lj`,

   ```
   E[bb'] = sum_l diag(alpha_l * mu2_l) + sum_{l != l'} m_l m_l''
   h2     = tr(R E[bb']) = mu'R mu - sum_l m_l' R m_l + sum_l s_l
   ```

   **The middle term is not optional.** An earlier attempt computed
   `mu'R mu + sum(sd^2)`, which keeps the positive between-effect term and drops
   the negative within-credible-set one; it over-counted badly (a spurious
   `+0.065` mean, and an apparent 60% "explanation" of the polygenicity bias
   that was entirely an artefact of the omission). Verified against a
   brute-force `tr(R E[bb'])` in `tests/test_susie_estimand.R`, and the `R = I`
   reduction is asserted on every fit at run time.

   **Result:** with both terms the two forms agree — median difference 0.000,
   mean −0.005, r = 0.995 over 3,000 paired fits. So the trace form is a good
   approximation here, and the polygenicity shortfall in `FINDINGS.md` §1 is a
   property of the sparse-effect model, not of the moment approximation.
   `susie_h2_quad` ships in every new run; older result files lack it.
14. **Dispersion across regions is not the same quantity as sampling SD.** The
   main grid runs ONE phenotype per (region, architecture), so the SD of
   `h2hat` over that grid mixes phenotype sampling variability with
   between-region LD heterogeneity. `grid.replicates` (see
   `configs/precision.yaml`) draws repeated phenotypes at a fixed region and
   architecture, which makes a conditional sampling SD estimable. Measured on
   30 regions x 25 draws: the two agree to <= 0.007 everywhere tested, so the
   across-region number is a good stand-in *here* — but that is an empirical
   finding, not a licence to conflate them elsewhere.
15. **`susie_win_rate` is a THREE-way share, not "SuSiE beats GCTA".** It is an
   alias for `win_SuSiE` from `analysis/metrics.py::win_rate`, i.e. the fraction
   of cells where SuSiE is closest of GCTA / SuSiE / HEELS — the three columns
   sum to 1. The head-to-head number is `pairwise_win_rate`, and the two differ
   materially: at one causal SNP the baseline grid gives 70% three-way against
   **76%** pairwise. The submitted abstract quotes the pairwise 76%, so anything
   phrased as "better than GCTA" must use `pairwise_win_rate`. The poster
   computes it directly and `poster/check_numbers.py` pins every value.

Deliberate departures from the previous code, all documented in
`cish2/config.py`: causal SNPs drawn *without* replacement; `gcta_h2` is
V(G)/Vp with raw V(G) kept separately (the old tables recorded V(G) and called
it h²); failures recorded as NaN with a status string rather than coerced to
sentinel values.

---

## Verification

Against the previous pipeline's `result_eur_num_people_100_300_regular_LD.txt`,
the rebuilt engine reproduces cell means with correlation **0.9987** (GCTA) and
**0.9995** (SuSiE); 93% and 100% of the 30 design cells fall within two
Monte-Carlo standard errors — the difference is sampling noise. Figure:
`results/baseline_n300/figures/fig0_legacy_reproduction.png`.

`tests/test_pipeline.py` covers the `.bed` reader against `pandas_plink`,
vectorised OLS against `statsmodels`, effect rescaling hitting the target
genetic variance, sampled genotypes reproducing their LD, the GCTA binary-GRM
layout, `.hsq` parsing, and the Python→R matrix round-trip.

---

## Not done yet

**Two one-line config changes, ~25 min each**, that would sharpen the main
result:

1. `-o method.ld_for_susie=sim` — separates "fine-mapping is biased" from "the
   LD matrix was rank-deficient". This is the single most valuable missing run.
2. `-o method.legacy_ambiguous_filter=false` — removes the ~15% SNP handicap.

**Larger:**

- **A fourth estimator** would need a fourth palette slot; three is what the
  validated categorical palette supports all-pairs, so a fourth arm means
  faceting rather than a fourth colour.

- **More GTEx tissues.** Would replace subsampling with a genuine sample-size
  axis and add a tissue dimension. The GTEx loader takes any
  expression/covariate pair.
- **Denser genotypes for the real-data pass.** HapMap3 gives ~365 SNPs per
  window against ~8,200 in the simulations, so real and simulated data are not
  in the same p-regime. Matching them would make the quantitative agreement
  meaningful rather than suggestive.
- **A calibrated SuSiE standard error.** The current analytic SE is the
  weakest part of the estimator, and §6 of `FINDINGS.md` says why.

---

## Environment

Runs on `mediator.sdsc.edu`: 2× AMD EPYC 7713 (256 cores), ~4 TB RAM, A100s
(unused — the workload is CPU and memory-bandwidth bound), **no scheduler**, so
jobs are launched with `nohup` and concurrency is self-limited (48 workers is a
polite share).

Conda env `dsc180` (`/projects/ps-renlab2/sunan/anaconda3/envs/dsc180`) carries
GCTA 1.94.1, PLINK 1.90, R with susieR/data.table/Rfast, and the Python stack.
It is Python 3.8 / pandas 2.0.3 — **no `pyarrow`, no `jsonlite` in R**, and
system `node` is v12. Result tables are therefore gzipped TSV, and data crosses
into R as raw binary rather than JSON.

Storage: `/projects/ps-renlab2` is NFS (~117 MB/s); `/data` is local NVMe
(~1.5 GB/s) and is where all transient LD/GRM traffic goes.
