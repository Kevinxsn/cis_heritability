# PI feedback, round 4 — design doc

**Status:** proposed · **Written:** 2026-10-28 · **Owner:** Sunan Xu
**Audience:** Sunan, and the coding agent who will implement it.
**Scope:** `poster/` (analysis, figures, layout, number pins) and one new
`cish2` feature: a semi-synthetic / permutation mode for the GTEx runner.

Earlier rounds: `poster/REVISION.md`, `REVISION2.md`, `REVISION3.md`. Read
`PROJECT.md` (especially traps 5, 6, 13–15) and `MIGRATION.md` before touching
estimators.

---

## 0. The feedback, mapped to the poster

| # | PI comment | Poster element | Action |
|---|---|---|---|
| 1 | "Really like the comparison of standard errors in simulations. Chews is much lower." | §5 panel A, *Dispersion of ĥ²* (`P3_sample_size`, left) | **Keep as is.** |
| 2 | "MAE difference is mostly blue but makes the differences seem marginal (except the lower left corner). Show how *frequently* Chews improves. % difference? Or does that not help because true cis-h² is constant?" | §5 panel B, *Mean absolute error difference* (`P3_sample_size`, right) | **Redesign** → §2 of this doc. |
| 3 | "Really love the split-half consistency." | §8 table | **Keep as is.** |
| 4 | "For stratification, show visually when Chews is performing *better* (not just differently) than GCTA. Is the underestimation in the bottom-left 0–0.05 stratum good or bad?" | §7, *A conditioning artefact* (`P6_strata_compare`) | **New analysis + redesign** → §3 of this doc. |

Items 1 and 3 are constraints: nothing below may change those two panels or
their numbers. `check_numbers.py` already pins them; it must still pass.

> **Ambiguity to confirm with the PI.** "Bottom left 0–0.05 stratum": in the
> current figure the 0–0.05 *column* is on the left, but the 0-signal *row* is
> at the top, and the bottom-left cell (4+ signals, 0–0.05) is n/a. This doc
> reads it as **the whole 0–0.05 column of panel B (disjoint strata)**, with the
> 0-signal cell as the main case. That cell holds 9,708 of the ~19,400 genes,
> so it is where most genes live, not a corner case.

---

## 1. Summary of decisions

1. **§5 panel B becomes a win-rate heatmap** (n × causal SNPs, "% of loci where
   CHEWS is closer to truth"), using the same colour scale as §4's `P2`. The
   relative MAE change (%) is the magnitude companion. No new runs; about an hour
   of work.
2. **Answer to the PI's % question:** h² is *not* constant in that panel (it
   pools six true h² values, 0.05–0.8). Per-fit "% of truth" error reweights
   toward h² = 0.05 and can **flip the sign** of a cell, so we don't use it.
   "% change in MAE" (MAE_CHEWS / MAE_GCTA − 1) is scale-free and safe. Win rate
   does not depend on any scaling at all.
3. **§7 needs ground truth, so we make some.** Run the GTEx split-half pipeline
   unchanged, except that the expression is replaced by phenotypes simulated
   on the **real GTEx genotypes** (same donors, same halves, same covariates,
   same SNP density). Strata are built exactly as in the real analysis. For
   each stratum we then report how often CHEWS is closer to the truth.
4. **A permutation null on real expression** anchors the null corner: it
   measures where each estimator's estimates sit when h² is truly 0.
5. **§7 figure:** panel A = the observed disjoint-donor gap (today's panel B);
   panel B = the simulated "% CHEWS closer to truth" in the same strata. The
   coupled-strata panel moves to the supplement and keeps one sentence in the
   text.

---

## 2. Item A — §5 panel B: show how *often* CHEWS wins

### 2.1 Why the current panel looks marginal

Numbers below come from `results/sample_size/cells.tsv.gz` (quick check run
2026-10-28; the script is at
`/expanse/lustre/scratch/sux002/temp_project/pi_feedback_metrics/metrics.py`.
That is purgeable scratch, so re-implement it properly in step A1).

- **Absolute differences shrink with the error itself.** MAE falls about 4×
  from n = 100 to 1000, so a constant *proportional* advantage looks like
  −0.17 at n = 100 and −0.03 at n = 1000.
- **The colour limit is set by the single −0.17 cell**, so every other cell is
  drawn in the pale centre of the ramp.
- **Mean differences hide how often CHEWS wins.** At one causal SNP, CHEWS is
  closer to truth on **75–82 % of loci at every n**. What shrinks with n is the
  size of the win, not how often it happens. The current panel cannot show this.

| causal SNPs ↓ / n → | 100 | 300 | 500 | 800 | 1000 |
|---|---:|---:|---:|---:|---:|
| **win rate (CHEWS closer to truth, paired)** | | | | | |
| 1  | 82 % | 81 % | 78 % | 78 % | 75 % |
| 3  | 65 % | 64 % | 62 % | 62 % | 64 % |
| 5  | 64 % | 57 % | 54 % | 61 % | 57 % |
| 10 | 53 % | 41 % | 37 % | 40 % | 43 % |
| 20 | 40 % | 28 % | 24 % | 23 % | 24 % |
| **relative MAE change, MAE_C / MAE_G − 1** | | | | | |
| 1  | −73 % | −65 % | −61 % | −57 % | −58 % |
| 3  | −49 % | −42 % | −32 % | −33 % | −36 % |
| 5  | −35 % | −16 % | −11 % | +4 % | −19 % |
| 10 | −5 % | +19 % | +31 % | +29 % | +14 % |
| 20 | +24 % | +91 % | +112 % | +94 % | +91 % |
| **absolute MAE difference (current panel)** | | | | | |
| 1  | −0.171 | −0.064 | −0.042 | −0.036 | −0.034 |
| 20 | +0.042 | +0.063 | +0.053 | +0.034 | +0.031 |

~290–300 paired fits per cell (50 regions × 6 h² values). There are no exact
ties.

**Be honest about what this shows.** The frequency view makes the wins at
≤ 5 causal SNPs clearer. It also makes the losses at 20 causal SNPs clearer:
GCTA is closer on 72–77 % of loci for n ≥ 300. That is the trade-off the poster
already states, so this is a clarity gain, not cherry-picking. Do not drop the
10/20 rows to make the panel look better.

### 2.2 Answering the PI's "% difference" question directly

Two different things could be meant by "% difference":

| metric | definition | verdict |
|---|---|---|
| per-fit relative error | mean of \|ĥ² − h²\| / h² | **Don't use.** Within one true h² it just rescales absolute error, so it adds nothing; the PI's intuition is right. Across the h² values pooled in a cell, it multiplies the h² = 0.05 fits by 20 and the 0.8 fits by 1.25, so a few low-h² cells dominate. Concrete failure: at 20 causal SNPs and n = 100, MAE says CHEWS is **worse** (+0.042), but mean relative error says **better** (−0.34), only because CHEWS wins 76 % of the h² = 0.05 fits. |
| relative MAE change | MAE_CHEWS / MAE_GCTA − 1 | **Use as the magnitude.** Scale-free: it divides by the competitor's error, not by h², so n = 100 and n = 1000 become comparable and the h² weighting is left alone. |
| win rate | P(\|ĥ²_C − h²\| < \|ĥ²_G − h²\|), paired | **Use as the headline.** It is exactly "how frequently Chews results in an improvement", and no rescaling, by h² or anything else, can change it. |

One sentence for the PI: *"h² isn't constant in that panel (it pools 0.05–0.8),
and dividing by it lets the h² = 0.05 fits dominate. Win rate answers 'how
often' directly, and % change in MAE answers 'by how much' without that
weighting."*

### 2.3 Design options considered

| option | encoding | pros | cons |
|---|---|---|---|
| A1 relative-MAE heatmap | colour = MAE_C/MAE_G − 1 | answers the % question literally | asymmetric scale (−73 % vs +112 %), needs a log-ratio colour map; still a magnitude, not a frequency |
| **A2 win-rate heatmap** ✅ | colour + label = win rate, centred at 50 % | literally answers "how frequently"; same encoding and colour scale as §4 `P2`, so the poster reads as one system | loses magnitude |
| A3 win rate + % MAE, two-line label | colour = win rate; label "82 %" over "−73 %" | both answers in one panel | two lines per cell in a 5×5 grid on an 11 in column may fall below the 20 pt annotation floor |
| A4 paired-difference violins | distribution of \|e_C\| − \|e_G\| per cell | most complete | 25 distributions in one column cannot be read from 1.5 m |

**Recommendation: A3 if it passes the legibility check, otherwise A2.**
Prototype A3 at true column width (`COLW_IN`). Accept it only if the win-rate
label is ≥ 19 pt bold and the MAE line is ≥ 16 pt. Otherwise ship A2 and put the
relative-MAE range in the caption. Only one cell has the two metrics pointing
in opposite directions (5 causal, n = 800: win rate 61 %, MAE +4 %), so A2 loses
very little.

### 2.4 Specification (for the coding agent)

**A1. Table — `poster/revision_analyses.py`, new `part_e()`**

- Input: `results/sample_size/cells.tsv.gz`, rows with both `susie_h2` and
  `gcta_h2` non-NaN. Same filter as `make_figures.pairwise_win`; reuse that
  function. Do **not** use the pipeline's three-way `susie_win_rate`
  (trap 15).
- For each `(ncausal, n)`: `win_rate`, `mae_chews`, `mae_gcta`,
  `rel_mae = mae_chews / mae_gcta − 1`, `n_pairs`.
- **95 % CI by region-clustered bootstrap** (2,000 resamples, resample
  `region_id`), reusing `make_figures.region_bootstrap_ci` (move it into a
  shared helper if importing it is awkward). The same region supplies every
  h² and causal count, so a binomial SE is too small.
- Write `poster/tables/win_by_n_causal.csv`, long format, with
  `rate, lo, hi, rel_mae, rel_mae_lo, rel_mae_hi, n_pairs`.
- Pooling over h²: average per-fit indicators with equal weight per fit. That
  is what the table above does and it matches the current panel's pooling.
  State it in the caption.

**A2. Figure — `poster/make_figures.py::fig_sample_size`, panel B only**

- Panel A: **unchanged**, byte-for-byte if possible.
- Panel B: `imshow` of win rate (%), rows = causal SNPs (20 at top, as now),
  cols = n. Use **the same colour map and limits as `fig_landscape` (`P2`)**
  so a given colour means the same win rate in §4 and §5. Read the limits from
  that function rather than hard-coding them. Centre the map at 50 %.
- Cell labels: `f"{rate:.0f}%"`, using the existing `cell_text` helper
  (dark/white switching). For A3, add a second line `f"{rel:+.0%}"` at the
  smaller size and check the legibility rule in §2.3.
- Optional: a thin outline on cells whose bootstrap CI excludes 50 %. This is
  the "how often is it real" cue. Leave it out if it clutters.
- Colour-bar label: **"% of loci where CHEWS is closer to truth"**.
- Panel title: **"How often CHEWS has the lower error"** (neutral wording,
  following round 3's style).

**A3. Poster text — `poster/build_poster.py`, §5 block**

Suggested replacement for the second half of the §5 paragraph (check the word
count; round 3 ended at 964 words):

> With one causal SNP CHEWS is closer to truth on 75–82 % of loci **at every
> sample size**; what shrinks with n is the margin (MAE −73 % at n = 100,
> −58 % at n = 1000). At twenty causal SNPs the ordering reverses: GCTA is
> closer on 72–77 % of loci for n ≥ 300.

Caption: *"(B) % of paired fits where CHEWS is closer to truth; pooled over
h² = 0.05–0.8, ~300 fits per cell, region-clustered bootstrap CIs."*

**A4. Number pins — `poster/check_numbers.py`, §5 block**

Pin every number quoted in the text: the 75 and 82 % endpoints at one causal
SNP, the 72–77 % GCTA share at 20 causal SNPs for n ≥ 300, and −73 % / −58 %.
Keep the existing MAE pins (0.142 / 0.205 / 0.039 / 0.044) if those numbers
still appear anywhere; delete pins only for text that was removed.

**Acceptance for item A**

- [ ] `check_numbers.py` passes, with the new pins added.
- [ ] Panel A pixels unchanged (`cmp` the panel crop, or at least eyeball it).
- [ ] Smallest type on panel B ≥ 16 pt at print size.
- [ ] The P2 and P3B colour bars are identical.

---

## 3. Item B — §7: *better*, not just *different*

### 3.1 What panel B can and cannot say today

Panel B stratifies genes by **half-A** GCTA ĥ² and half-A independent-signal
count, then plots mean(CHEWS − GCTA) in **half B**. Because no donor is in
both halves, half B's noise does not depend on the stratum. So if GCTA were
unbiased, E[GCTA_B | stratum] would equal E[h² | stratum], and the gap would
read directly as CHEWS's bias. **Two things break that reading in exactly
the 0–0.05 column:**

1. **GTEx GCTA is *constrained* REML** (`gcta_reml_no_constrain: false`;
   unconstrained REML failed on 60–75 % of genes, see `PROJECT.md`). A
   non-negative estimator is biased **upward** when the true h² is near 0, so
   GCTA is not a neutral reference in this column. The simulations used
   *unconstrained* REML, so "GCTA is unbiased" from §4 does not carry over.
2. **"Better" is about per-gene error, not mean bias.** A mean gap of −0.004
   is consistent with CHEWS being closer on most genes, or on few. Only
   per-gene truth can tell these apart.

Real-data facts for the 0–0.05 column (half B, disjoint strata; same quick run):

| half-A signals | genes | mean GCTA_B | mean CHEWS_B | CHEWS_B < 0.005 | GCTA_B at 0 bound |
|---|---:|---:|---:|---:|---:|
| 0 | 9,708 | 0.0257 | 0.0213 | 28 % | 3 % |
| 1 | 1,205 | 0.0467 | 0.0429 | 14 % | 2 % |
| 2 | 50 | 0.0744 | 0.0614 | 14 % | 2 % |

For context, the median per-gene GCTA SE in a half is 0.049, about **twice the
stratum's mean estimate**. CHEWS sets 28 % of these genes to ~0. Whether that is
right depends on **how many of those 9,708 genes are truly null** (CHEWS is
correct, GCTA's floor is too high) or truly carry small, diffuse h² (CHEWS is
throwing real signal away). Real data cannot settle this, because any real-data
"reference" is just another estimator.

**Rejected proxies, and why:**

- Full-sample (n = 574) GCTA as truth: it shares donors with half B, which
  re-creates the coupling we just removed. It is also still constrained.
- Out-of-sample prediction R²: a lower bound that depends on the fine-mapping
  weights, which tilts it toward CHEWS.
- Agreement with HEELS: HEELS tracks GCTA at r ≈ 1, so it is not independent.

### 3.2 Approach: three pieces, one question

| piece | what it answers | new compute |
|---|---|---|
| **B1 Semi-synthetic GTEx (main)** | Per stratum: how often is CHEWS closer to truth, and what are the bias and MAE of each method? | 1 GTEx-scale run (~15–30 min on one node) |
| **B2 Permutation null** | Where does each estimator sit when h² = 0 on *real* expression and covariates? | 1 GTEx-scale run (~7 min) |
| **B3 Calibration by reweighting** | Which mixes of null / sparse / diffuse genes reproduce the *observed* panel-B pattern, and do the B1 conclusions hold across all of them? | none (post-hoc) |

#### B1. Semi-synthetic GTEx

Keep **everything** about `configs/gtex_split_half.yaml`: donors, the
`disjoint_subsamples` split (same seed), covariate handling, residualisation,
in-sample LD, clumping, constrained REML. Change only the phenotype: for each
gene window, simulate

```
Z_full = standardize(impute(raw_geno))                    # all 574 donors
beta   = sparse effects on HM3 SNPs in the window;
         rescaled so var(Z_full @ beta) = h2_target
y_full = Z_full @ beta + N(0, 1 - h2_target)              # then standardised
```

then hand `y_full` to the unchanged `_prepare_cell` in place of
`raw_expr`. Both halves see **the same β**: truth is a property of the
gene, as in the real-data question.

**Design prior** (deliberately broad, so B3 can reweight to any target prior):

- `h2_target ∈ {0, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.4}`. The low end is
  dense on purpose: that is the column in question.
- `ncausal ∈ {1, 2, 3, 5, 10, 20}` (ignored when h² = 0).
- Null mass 25 %; the rest uniform over the 7 × 6 non-null grid.
- `effect_dist: normal`, causal SNPs drawn without replacement (matching the
  main simulations).
- **K = 2 independent draws per gene** → ~39k simulated genes × 2 halves
  (~900 per non-null design cell). Increase K if B3's effective sample size per
  stratum falls below 30 (same masking rule as the real heatmap).

**Truth definition:** primary = `h2_target`, the population parameter, which
matches the question being asked. Also record, per half,
`h2_sample = var(g_half)/var(y_half)` and
`h2_resid = var(resid(g_half))/var(resid(y_half))`. The genotype PCs in the
covariate block absorb some genetic variance, and the summary-statistic
estimators only ever see the residualised scale. Report the primary metric
against `h2_target` and check that conclusions don't change with `h2_resid`.

**GCTA, both ways:** for the simulated run, fit REML **constrained and
unconstrained** (two `GctaRunner` instances; columns `gcta_h2` and
`gcta_h2_uc`). This separates "the constraint drives the 0–0.05 gap" from
"GCTA is noisy there". Expect unconstrained failures; record them as NaN with a
status, as elsewhere.

Known limitation to state on the poster: causal variants are drawn **among
the genotyped HM3 SNPs**, so tagging is better than in reality. This helps both
methods, and probably CHEWS more.

#### B2. Permutation null

Real expression and real covariates, with **genotype rows permuted within
each half** (seeded per gene). The expression–covariate relationship stays
intact; every gene's true cis-h² becomes 0. Output: the null distribution of
ĥ² for GCTA (constrained) and CHEWS: mean, median, and % at the floor, overall
and within the permuted half-A strata. That answers "where is each method's
floor" on real expression distributions, without relying on any simulation
model.

It also gives an **empirical null for GCTA's LRT p-value**, which yields a
Storey-type estimate of the null fraction π₀ in each real stratum:
π̂₀ = P̂_real(p > λ) / P̂_perm(p > λ), λ = 0.5. Use it as a prior anchor in B3.
(Constrained-REML p-values have a point mass under the null, which is why we
use the permutation null instead of Uniform(0,1).)

#### B3. Calibration by importance reweighting

B1 is run once under the broad design prior. Any target prior
π(h², ncausal) on the same grid can then be evaluated **without re-running**,
by weighting each simulated gene with
`w = π_target(h², k) / π_design(h², k)`.

Parameterise the target prior with three knobs:

- `π0`: null fraction (grid 0–0.9; anchored by B2's π̂₀),
- `h2_scale`: mean of a discretised exponential over the non-null h² grid,
- `sparse_weight`: geometric weight on small `ncausal`.

Score each prior on how well its **simulated** half-A stratum gene counts
(`tables/gap_disjoint_counts.csv`) and half-B gap means
(`tables/gap_disjoint_strata.csv`, SEMs from `gap_disjoint_sem.csv`) match the
real ones: a weighted squared error across the 25 cells, with cells masked
exactly as in the real figure. Keep the **set of acceptable priors** (within a
χ²-style tolerance of the best fit), not just the single best one.

Report for every stratum: the win rate, MAE_C − MAE_G, the bias of each method,
and the effective sample size (Kish), as **best-fit value + range over the
acceptable set**.

### 3.3 Decision rules, fixed before looking at B1

Write these into the analysis script's docstring before running, so the
interpretation can't drift after we see the results.

| outcome in the 0–0.05 column | conclusion we may print |
|---|---|
| Calibration fails: no prior reproduces the observed gap pattern (cell-wise sign agreement < 80 % or the gap means badly missed) | **Do not put B1 on the poster.** Say "simulation could not reproduce the real pattern" in the limitations and keep §7 as it is. |
| Win rate > 50 % with CI excluding 50 % across the whole acceptable set | "CHEWS's lower values there are **closer to truth**: GCTA's constrained floor sits above it." |
| Win rate < 50 % with CI excluding 50 % across the whole set | "CHEWS's underestimation there is a real **loss**: it zeroes genes with small diffuse h²." |
| The acceptable set straddles 50 % | "**Depends on the null fraction**, which the data do not pin down": show the win rate as a function of π₀ (B4 below), and give B2's π̂₀. |

Apply the same rules to the other strata, but the poster sentence should be
about the 0–0.05 column (the PI's question) plus one about the high-h² /
many-signal corner.

### 3.4 Figure specification — `P6_strata_compare` → `P6_better_or_different`

Two stacked panels, same 5 × 5 grid, same axes and labels as today:

- **A — "Observed: CHEWS − GCTA (strata from a disjoint half)"**: today's panel
  B, unchanged (diverging map, ±0.06).
- **B — "Simulated on GTEx genotypes: % genes where CHEWS is closer to
  truth"**: the best-fit-prior win rate per stratum, with the **same colour map
  and limits as `P2` and the new P3B** (a third use of one encoding, so the
  reader learns it once). Label `f"{rate:.0f}%"`. **Hatch** cells where the
  acceptable-prior range crosses 50 %. Mask cells whose effective sample size
  is < 30.
- The coupled-strata panel goes to the QR supplement. The §7 text keeps one
  sentence: *"Defining strata on the same donors manufactures much of the
  pattern (mean |coupled − disjoint| = 0.023); panel A uses disjoint donors."*

Optional **inset or panel C** if space allows. It answers the PI's question
for one column:

- **B4. "The 0–0.05 stratum, with truth"**: for the 0-signal and 1-signal
  cells, dots for mean truth, mean GCTA, mean CHEWS (simulated, best-fit
  prior), plus the B2 null floor as a dashed line per method. Or, if §3.3's
  third row applies: the win rate against π₀ as a line, with B2's π̂₀ marked.

Section title: **"7 Better, or just different?"** (replaces "A conditioning
artefact").

Text template (fill from the tables; every number pinned in `check_numbers.py`):

> In genes GCTA places at h² < 0.05, CHEWS reads lower (−0.004). Simulating
> phenotypes on the same GTEx genotypes, under null/sparse/diffuse mixes that
> reproduce panel A, CHEWS is closer to truth in **X %** of those genes
> (range Y–Z). Under permuted genotypes (true h² = 0), constrained GCTA averages
> **a** and CHEWS **b**.

### 3.5 Implementation — `cish2` changes

All changes are additive and default-off. A run with an existing config must
produce identical output.

**Config (`cish2/gtex_campaign.py` defaults):**

```yaml
phenotype_mode: real          # real | simulated | permuted
sim_design:                   # used only when phenotype_mode == simulated
  h2_grid: [0, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.4]
  ncausal_grid: [1, 2, 3, 5, 10, 20]
  null_mass: 0.25
  draws_per_gene: 2
  effect_dist: normal
gcta_also_unconstrained: false
```

New configs: `configs/gtex_sim_split_half.yaml` and
`configs/gtex_perm_split_half.yaml`. Each copies `gtex_split_half.yaml` and
sets the mode. **The split seed must be identical** (`seed: 20260804`) so the
halves are the same donors.

**Runner (`cish2/gtex_runner.py::GeneWorker._run_chunk`):** the hook point is
right after

```python
raw_geno = self.data.genotypes(win)
raw_expr = self.data.expr[win.expr_row, :]
```

- `simulated`: for `d in range(draws_per_gene)`, draw `(h2, k)` from the design
  prior and β via a new helper `simulate.sim_beta_empirical(Z_full, k, h2, rng)`.
  It reuses `sim_beta`'s effect-drawing logic but rescales empirically with
  `var(Z_full @ beta)`, because there is no Cholesky factor here and no `L`.
  Build `y_full`, then loop over subsets as usual with `raw_expr=y_full`. Label
  cells `f"{label}|d{d}"` so the two draws stay distinct, and add the columns
  `draw, h2_target, ncausal_target, h2_sample, h2_resid`.
- `permuted`: per subset, a seeded permutation `perm` of `idx` is applied to
  the **genotype rows only**. That needs `_prepare_cell` to accept a separate
  genotype index; add an optional `geno_idx=None` argument that defaults to
  `idx`.
- **RNG seeding:** `np.random.default_rng([seed, crc32(gene_id), d])`. Results
  must not depend on chunking or worker order (same principle as the main
  campaign).
- `gcta_also_unconstrained`: a second `GctaRunner(reml_no_constrain=True)`; its
  results go into `gcta_h2_uc`, `gcta_h2_uc_se`, `gcta_uc_status`.

**Tests (`cish2/tests/test_pipeline.py`):**

1. `sim_beta_empirical`: for h² = 0.3 and k ∈ {1, 5}, `var(Zβ)` = 0.3 to 1e-10,
   and exactly k non-zero entries.
2. Simulated mode: the same β is used for both halves (check via the stored
   causal indices or a hash).
3. Permuted mode: the expression vector handed to REML is bit-identical to
   real mode, and only the genotype rows differ.
4. **Regression:** `phenotype_mode: real` on one chunk reproduces the existing
   `results/gtex_split_half` rows for those genes to 1e-8. The GTEx genotypes
   did not change in the migration, so this should hold exactly; if it doesn't,
   stop and investigate before running anything else.

**Running on Expanse** (see `cish2/scripts/expanse/config.sh`; account
`sdp202`):

- **Nothing heavy on the login node.** `import pandas` alone hung for > 100 s
  there on 2026-10-28. Even the poster scripts go through `sbatch` (the
  `debug` partition, 30 min, is enough for those).
- Add `cish2/scripts/expanse/gtex_sim.sbatch` and `gtex_perm.sbatch`,
  modelled on the existing scripts. The real split-half pass took **6.4 min on
  48 workers** for 38,864 cells; the simulated run is ~2× (K = 2) plus a second
  REML fit, so request one node for 1 h, with ~1 h of slack. GTEx inputs are
  the read-only symlinks under `$DATA/gtex`; write results under
  `results/gtex_sim_split_half/` and `results/gtex_perm_split_half/`.
- Smoke test first: `max_genes: 200` on `debug`.

**Analysis (`poster/revision_analyses.py`):**

- `part_f()`: B2 null summary + π̂₀ per real stratum → `tables/perm_null.csv`,
  `tables/pi0_by_stratum.csv`.
- `part_g()`: B1 + B3 reweighting, calibration grid, per-stratum metrics →
  `tables/sim_strata_win.csv` (best fit + acceptable range + ESS),
  `tables/sim_strata_bias.csv`, `tables/calibration_fit.csv`.
  Reuse `H2_BINS`, `H2_LABELS` and `clump_bin` from `part_b` **exactly**, so the
  strata match the real ones.

### 3.6 Acceptance for item B

- [ ] Tests 1–4 pass; the real-mode regression matches to 1e-8.
- [ ] The simulated run reproduces the real per-stratum gene counts to within
      the calibration tolerance (if not, follow §3.3 row 1).
- [ ] The decision rule in §3.3 is written into `part_g`'s docstring
      **before** the B1 output is read.
- [ ] Every number on the new §7 is pinned in `check_numbers.py`.
- [ ] `REVISION4.md` written in the style of `REVISION3.md`: point by point,
      including any result that goes against CHEWS.

---

## 4. Layout and budget

- Panel A of §5, §8, and the rest of the poster: untouched.
- §7's figure keeps its footprint (two stacked 5×5 panels). An optional B4
  inset only goes in if the packer in `build_poster.py` still fits without
  shrinking type. It fails loudly otherwise; don't fight it, drop B4 to the
  supplement.
- Word budget: round 3 ended at 964. §5 text is roughly neutral; §7 grows by
  ~40 words. Trim §7's existing coupling paragraph to one sentence to pay for
  it.
- Rebuild order is unchanged:
  `revision_analyses.py → make_figures.py → build_poster.py → check_numbers.py`,
  all via `sbatch` on Expanse.

## 5. Execution order

| phase | work | depends on | rough effort |
|---|---|---|---|
| 1 | Item A: A1–A4, rebuild, send the PI the new §5 alone | — | ~1 h |
| 2 | `cish2` modes + tests + smoke test | — | ~half a day |
| 3 | B2 permutation run, then B1 simulated run | 2 | ~1 h wall each, mostly queue |
| 4 | `part_f`, `part_g`, calibration | 3 | ~half a day |
| 5 | §7 figure + text + pins + `REVISION4.md` | 4 | ~2 h |

Phase 1 is independent and quick; ship it to the PI first.

## 6. Open questions (Sunan, or the PI)

1. **Which "bottom-left" cell did the PI mean?** (See §0.) The plan covers the
   whole 0–0.05 column, so it's robust either way, but the poster sentence
   should name the cell they had in mind.
2. **Is this for a reprint of the poster or for follow-up slides?** Slides have
   room for B4 and for the coupled panel; the poster does not.
3. **Is it OK to drop the coupled panel from the poster?** It was the subject
   of round 1's correction, so check the PI doesn't want it kept for that
   reason.
4. **How much prior dependence will the PI accept on the poster?** If §3.3
   ends in "depends on π₀", is a π₀ curve acceptable, or should that stay in
   the supplement?
5. **Should the unconstrained-REML GCTA (simulated run only) also appear?**
   It is the cleanest test of "the constraint makes GCTA's floor too high", but
   it adds a fourth series.
