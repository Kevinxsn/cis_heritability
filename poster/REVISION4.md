# Response to the PI review (round 4)

Point by point against the four comments mapped in
[`../docs/pi_feedback_round4_design.md`](../docs/pi_feedback_round4_design.md) §0.
Earlier rounds: [`REVISION.md`](REVISION.md), [`REVISION2.md`](REVISION2.md),
[`REVISION3.md`](REVISION3.md). Every number below is re-derived by
`check_numbers.py --round 4` (all 76 round-3 checks plus 15 new ones, all passing).

The round-4 poster is a separate deck, `XuS_IGVF2026_poster_round4.pptx`; the
round-3 deck is unchanged and still rebuilds with the round-3 commands.

---

## 1. "Really like the comparison of standard errors" — **kept as is**

§5 panel A is drawn by the same code as round 3, at the same position: the
left 40% of the new figure is pixel-identical to the round-3 figure (0
differing pixels).

## 2. "Show how frequently CHEWS improves. % difference?" — **redesigned**

**Panel B now shows how often CHEWS is closer to truth**, by causal SNPs × n,
on the same colour scale as the §4 win-rate heatmap, so one colour means one
win rate across the poster.

| causal SNPs ↓ / n → | 100 | 300 | 500 | 800 | 1000 |
|---|---:|---:|---:|---:|---:|
| 1  | 82% | 81% | 78% | 78% | 75% |
| 3  | 65% | 64% | 62% | 62% | 64% |
| 5  | 64% | 57% | 54% | 61% | 57% |
| 10 | 53% | 41% | 37% | 40% | 43% |
| 20 | 40% | 28% | 24% | 23% | 24% |

290–300 paired fits per cell, pooled over the six true h² values (0.05–0.8)
with equal weight per fit. 95% CIs resample whole regions (2,000 draws); 23 of
25 cells exclude 50%, and the two that do not (10 causal at n = 100, 5 causal at
n = 500) carry a dotted frame.

**On "% difference": h² is not constant in that panel.** It pools 0.05–0.8, and
dividing each fit's error by its h² lets the h² = 0.05 fits dominate. At 20
causal SNPs and n = 100 that would even flip the sign of the cell. Win rate
answers "how often" with no rescaling at all. The scale-free magnitude is the
change in MAE, MAE_CHEWS / MAE_GCTA − 1, which runs from −73% (1 causal SNP,
n = 100) to +112% (20 causal SNPs, n = 500). It is in
`tables/win_by_n_causal.csv`, with CIs.

New §5 text:

> With one causal SNP CHEWS is closer to truth on **75–82% of loci at every
> sample size**; what shrinks with n is the margin (MAE −73% at n = 100, −58%
> at n = 1,000). At twenty causal SNPs GCTA is closer on 72–77% of loci for
> n ≥ 300.

**This view also makes the loss clearer, and it stays on the poster.** At
20 causal SNPs GCTA wins three loci in four once n ≥ 300. That is the
trade-off the poster already states, now shown as a frequency.

We also tried putting the MAE change under each win rate (design option A3).
The two-line labels collided across cells at the required ≥ 16 pt, so the
poster ships win rate only (A2). The rejected prototype is
`figures/P3_sample_size_r4_A3_rejected.png`.

## 3. "Really love the split-half consistency" — **kept as is**

§8 is untouched; its nine correlations are still pinned and pass.

## 4. "Show when CHEWS is *better*, not just different. Is the underestimation in the 0–0.05 stratum good or bad?" — **partly answered; the simulation failed its own calibration test**

Real data cannot say which estimator is closer to truth, so we ran two new
GTEx passes on the same donors, halves, covariates and estimators.

**B2: permutation null (model-free).** We permuted genotype rows within each
half, so every gene's true h² is 0 and the expression is real:

| estimator | mean ĥ² | median | at the floor* |
|---|---:|---:|---:|
| GCTA (constrained REML) | 0.013 | 0.000 | 57% |
| CHEWS | 0.013 | 0.009 | 40% |
| HEELS | 0.011 | 0.000 | 64% |

\*GCTA/HEELS at their 0 bound; CHEWS < 0.005.

**At true h² = 0 the two estimators average the same.** GCTA's
constrained floor is not higher on average; GCTA reaches the bound more often,
and CHEWS sits slightly above 0 more often. So "GCTA's floor is too high" is
not, on its own, why CHEWS reads 0.004 lower in that column. This goes on the
round-4 poster (one sentence in §7).

From the same run we estimated the share of truly null genes, π̂₀ (the fraction
of genes whose REML hits the bound, relative to the permutation null):
**0.73 in the 0-signal, h² < 0.05 cell** (9,572 genes), 0.49 in the 1-signal
cell, and 0.50 genome-wide (0.55–0.62 under other cut-offs).

**B1: semi-synthetic GTEx (needs calibration).** Each gene got two simulated
phenotypes on its own real genotypes, with known h² ∈ {0, …, 0.4} and 1–20
causal SNPs, re-weighted to the mix of null, sparse and diffuse genes that best
reproduces the real §7 panel. The rule fixed **before** looking was: if the
best mix misses the real gap means by more than 0.010 on average (or gets the
sign wrong in more than 20% of cells), the simulation does not go on the
poster.

**It missed: 0.0139.** Signs agree in 94% of cells and the stratum sizes match
(total-variation distance 0.022). The miss sits almost entirely in the
high-h² / many-signal corner, where simulated gaps are 0.03–0.05 too
negative. The likely cause is that the simulation's h² grid stops at 0.4 while
the real .4–1 stratum holds genes above that. The best-fitting mix also wants
π₀ = 0.1, against the permutation estimate of ~0.5–0.7. **Per the rule, §7
keeps its round-3 figure, and §10 gains: "Semi-synthetic GTEx simulation
missed the real strata beyond its pre-set tolerance; §7 is disagreement, not
error."**

**What the uncalibrated simulation suggests, for discussion only and not on
the poster.** In the 0–0.05 column the simulation does reproduce the real gaps
(0-signal cell: −0.0046 simulated, −0.0044 real). There:

- CHEWS is closer to truth on **60% of genes (95% CI 59–61%)** in the 0-signal
  cell and 63% (61–65%) in the 1-signal cell.
- This holds at **every** null fraction from 0 to 0.9: 54–61% in the 0-signal
  cell, ~55% at the permutation-anchored π₀ ≈ 0.7. The CI excludes 50% throughout.
- Constrained GCTA is biased upward there (+0.009 vs CHEWS +0.005). The
  constraint does push GCTA up by roughly the size of the observed gap (mean
  0.026 constrained vs 0.020 unconstrained), but unconstrained REML is much
  noisier: CHEWS beats it on 77% of genes.
- **Against CHEWS:** at 4+ signals and h² > 0.4, GCTA is closer (CHEWS 36%,
  CI 33–40%), and at 3–4+ signals with h² 0.1–0.4 the answer depends on the mix.

So the provisional answer to "good or bad?" is **"mildly good, in the column
you asked about"**: CHEWS's lower values are closer to truth a bit more often
than not, but by a modest margin, and not from evidence we have yet allowed
onto the poster.

**To make it poster-grade:** extend the simulated h² grid to 0.8 and re-run
(≈17 min on one debug node), keeping the same pre-registered rule. That is a
change to the design after seeing a result, so it is flagged as such and
needs a decision.

---

## Questions back to you

1. "Bottom-left 0–0.05 stratum": did you mean the 0-signal cell (9,572 genes)?
   The analysis covers the whole column either way.
2. Is the extended-grid re-run worth doing before printing, or should the
   simulation stay a follow-up?
3. Is the permutation sentence in §7 the right level of claim for the poster?
