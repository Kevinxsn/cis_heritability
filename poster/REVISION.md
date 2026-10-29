# Response to PI review — what changed, and what did not

Point-by-point against `XuS_IGVF2026_poster_feedback.md`. Numbers cited here are
re-derived by `check_numbers.py` (54 checks, all passing).

---

## Priority 0

### 1. GTEx values vs the submitted abstract — **resolved, abstract was wrong**

Not two analyses to choose between: the abstract's numbers came from a run with a
degrees-of-freedom bug. It pre-residualised expression on covariates and handed
REML the residual, which hides the df the covariates consume. With 24 covariates
at *n* = 100 that is a quarter of the sample, so GCTA was inflated at small *n*
and the apparent gap that "closed with sample size" was largely an artefact of
its own construction.

Covariates now go to REML as fixed effects (`--qcovar`). The poster carries an
explicit **"Updated analysis"** callout in §6 naming the cause in one sentence.

| quantity | abstract | poster (corrected) |
|---|---:|---:|
| mean CHEWS − GCTA at *n* = 100 | −0.151 | **−0.015** |
| mean CHEWS − GCTA at *n* = 574 | −0.023 | **−0.007** |
| CHEWS–GCTA correlation | 0.61 → 0.88 | **0.81 → 0.90** |

Gene counts reconciled and now stated on the poster: 19,432 genes attempted at
every subset; GCTA returns on 96.5–98.6% (CHEWS and HEELS on all); the 19,169 is
the *n* = 574 complete-case set; the 4,171 is the post-hoc reliability set, now
shown alongside an all-genes and a prespecified set (§8).

### 2. The CHEWS estimand — **stated exactly, and measured**

The review is correct and this was worth catching. The implementation is
literally `sum(mu^2) + sum(sd^2)` (`R/susie_batch.R`), i.e. the **trace** form.
The target is

```
E[β'Rβ | data] = μ'Rμ + tr(RΣ)
```

With a diagonal posterior covariance and unit-diagonal `R`, `tr(RΣ) = Σⱼσⱼ²`
exactly — so the diagonal form and the quadratic form differ by precisely
`μ'Rμ − μ'μ`: the cross terms between credible sets that are themselves in LD.
These vanish for a single resolved effect or mutually uncorrelated sets, and not
otherwise. §2 of the poster now says this.

Rather than assert the term is small, `susie_batch.R` was extended to emit both
(one extra matrix-vector product against an `R` already in memory) and the
baseline grid was re-run (`configs/estimand_check.yaml`, same seed, same
phenotypes). **Complete, 3,000 paired fits.** The trace column reproduces the
original run at r = 1.000000, so the two forms are compared on identical
phenotypes.

| causal SNPs | bias, trace | bias, quadratic | omitted term |
|---:|---:|---:|---:|
| 1 | −0.008 | **+0.031** | +0.039 |
| 3 | −0.028 | +0.029 | +0.056 |
| 5 | −0.044 | +0.026 | +0.070 |
| 10 | −0.081 | +0.003 | +0.084 |
| 20 | **−0.130** | **−0.052** | +0.079 |
| overall | −0.058 | +0.007 | **+0.065** |

**The term is not negligible, and this qualifies the poster's central mechanism.**
It accounts for ~60% of the trace form's shortfall at 20 causal SNPs. So a large
part of what we had been attributing to the sparsity assumption is the trace
approximation in the estimator's implementation.

Two things keep this from being a simple "use the quadratic form" fix:

1. Restoring the term **over-corrects at low polygenicity** (+0.031 at one causal
   SNP against −0.008 for the trace form), so overall MAE is a wash: 0.074 trace
   against 0.076 quadratic. What changes is the composition — the quadratic form
   is much better at 20 causal SNPs (MAE 0.091 vs 0.135) and worse at 1 (0.058
   vs 0.037).
2. `R` here is the **rank-deficient reference panel** (503 individuals,
   ~8,200 SNPs), so `μ'Rμ` is evaluated against a mis-specified matrix. The
   upward bias is plausibly that rather than a property of the estimand.

Neither form is therefore the right estimator as it stands, and the poster says
so. §4's figure now shows all four curves — GCTA, HEELS, CHEWS trace, CHEWS
quadratic — so a viewer can see exactly how much of the shortfall each
explanation carries. **This is the most important thing the review surfaced.**

### 3. HEELS is not a clean sparsity control — **conceded on the poster**

The two statements were indeed inconsistent. §5 now carries a callout saying
plainly that HEELS received in-sample LD while CHEWS received reference-panel LD,
so they differ in model **and** LD source, and that the factorial run is in
progress. The claim "attributable to the sparsity assumption alone" is gone.

HEELS's *n* = 100 blow-up (+0.18 bias, MAE 0.30) is left in the figure and
labelled: at p/n ≈ 82 it is far outside its published n ≫ p regime.

### 4. Accuracy / precision / agreement / reliability — **defined and separated**

§1 opens with an operational definition callout. The title is unchanged (it is
the accepted abstract), but "precision" is now defined where a reader meets it.

Correlation is no longer the agreement claim. New `P5_agreement` reports mean
difference, Deming calibration slope and RMSE:

| | r | slope | RMSE |
|---|---:|---:|---:|
| HEELS vs GCTA, *n* = 574 | 0.9996 | **0.932** | 0.0116 |
| CHEWS vs GCTA, *n* = 100 | 0.815 | **0.666** | 0.079 |
| CHEWS vs GCTA, *n* = 574 | 0.897 | **1.064** | 0.058 |

"HEELS and GCTA are one method" is gone; the poster says HEELS tracks REML
closely **but not exactly**, slope 0.93.

The "a constant estimator scores a perfect ICC" line is **removed**. It was
wrong — a constant has no between-gene variance, so ICC is undefined, not
perfect. Replaced with: reliability rewards low sampling noise *relative to
between-gene variance*, which shrinkage delivers.

### 5. Sparsity inference — **softened to the reviewer's wording**

"Whole Blood cis-architecture is mostly sparse" is gone. §6 now reads: consistent
with concentrated cis effects under the CHEWS model, but does **not** establish
that whole blood is sparse.

### 6. Heatmap conditioning — **the concern was right; the claim is withdrawn**

Re-ran the stratification with strata defined on a **disjoint half of donors**
(half A) and the gap measured in the other half (half B). No shared donors, so
no shared noise. `P6_strata_compare` shows both side by side, and the pattern
largely does not survive:

- the positive low-h² cells (+0.026, +0.032) become −0.004, −0.013;
- the −0.099 extreme becomes +0.003;
- mean |coupled − disjoint| = **0.023**, larger than most of the cell values.

**The Simpson's-paradox claim and the "diffusely regulated undetected-eQTL genes"
interpretation are withdrawn from the poster.** What survives is a mild gap that
scales with h². Per-cell gene counts are in `tables/gap_disjoint_counts.csv` and
the range is quoted in the caption.

### 7. Split-half selection bias — **both analyses now shown**

Reliability is reported on three gene sets, one of which is chosen entirely
independently of the split-half estimates (≥1 cis-eQTL in the full 574 donors):

| gene set | GCTA | CHEWS | HEELS |
|---|---:|---:|---:|
| all genes (18.5–19.4k) | 0.865 | 0.912 | 0.870 |
| prespecified, ≥1 cis-eQTL (9,665) | 0.850 | 0.898 | 0.857 |
| post-hoc, ĥ² > 0.05 ×6 (4,171) | 0.848 | 0.892 | 0.855 |

The ordering is stable, so it is not a filter artefact. The bare 0.835-vs-0.848
comparison is replaced by a paired bootstrap: CHEWS − GCTA split-half r =
**+0.050, 95% CI +0.038 to +0.063** over 2,000 resamples. Still one split; the
poster says so.

---

## Priority 1

- **Fit accounting** (§3): 20,400 cells attempted, 20,294 paired, unconstrained
  REML failed on 106. CHEWS and HEELS returned on all; CHEWS left [0,1] on <0.3%.
  GCTA's constrained/unconstrained status is stated in both sections.
- **Out-of-range values**: the poster says *screen*, and reports the failures
  rather than presenting silent removal as routine.
- **Standard errors**: narrowed to "nominal 95% Wald intervals under-covered
  here (68–87%, by estimator and condition) — report method-specific coverage."
- **Decision guide** (§9) is now conditional on architecture as well as *n*, and
  the MAF row says to validate on discrete genotypes first.
- **Simulation caveats** are in §10, including that regions recur across grid
  cells so Monte-Carlo intervals are not independent.

**Not done:** region-clustered bootstrap; validating the MAF result on sampled
or haplotype genotypes; many random half-splits. These are stated as limitations
rather than quietly omitted.

---

## Priority 2 — design

- **12 sections → 10, in 3 clear modules** (question/estimand/benchmark ·
  simulation · real data + practice). Prose cut from ~1,240 to ~640 words.
- **Type scale raised throughout**: title 72 pt (was 52), section headers 33 (27),
  body 25 (19–23), captions 21 (17), tables 22 (19–20).
- **Every figure is now drawn at the column width (11 in).** This was the real
  cause of unreadable annotations: a figure drawn 22 in wide and placed in an
  11 in column has all of its type halved. `P5` and `P6` annotations were landing
  at ~8 pt on the printed sheet. They are now true-size.
- **Full-width header** carries a one-sentence take-home and a CHEWS workflow
  strip (summary stats + LD → SuSiE posterior → posterior genetic variance → h²).
- **Neutral figure titles**; direct labelling instead of legend boxes where a box
  would cover data; ALL-CAPS emphasis replaced with restrained bold.
- **QR code** in the footer (placeholder URL — see below).
- Detailed robustness figure and the MAF/LD sensitivity panels moved off the
  poster, to be linked from the QR supplement.

---

## Before printing

1. **Wait for `estimand_check`** and fold the measured cross-term size into §2.
2. Point the QR at the real repository and regenerate `figures/qr.png`.
3. Fill the acknowledgements placeholder (grant numbers) in §10.
4. Export to PDF at 100% and check fonts and equations.
