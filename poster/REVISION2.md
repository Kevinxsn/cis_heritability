# Response to the second PI review

Point-by-point against `XuS_IGVF2026_poster_feedback_round2.md`.
Round-1 response is in [`REVISION.md`](REVISION.md).

---

## Priority 0

### 1. The quadratic form — **the review was right, and my first attempt was wrong**

This is the decisive item and the reviewer's instinct was correct on both counts:
Σ is *not* diagonal, and assuming it was produced a wrong answer.

**What SuSiE's posterior actually is.** The variational posterior factorizes
*across* single effects, and within one single effect exactly one SNP is
non-zero. Writing `m_l = α_l ⊙ μ_l` and `s_l = Σ_j α_lj·mu2_lj`:

```
E[b]    = Σ_l m_l                                       ( = coef(fit)[-1] )
E[bb']  = Σ_l diag(α_l ⊙ mu2_l)  +  Σ_{l≠l'} m_l m_l'ᵀ
```

so the genetic variance is

```
h²_full = tr(R·E[bb']) = μᵀRμ  −  Σ_l m_lᵀ R m_l  +  Σ_l s_l
```

The middle term is a **within-single-effect** correction, and it is exactly the
`Σ_{j≠k} R_jk Σ_kj` contribution the review asked about. Individually it is
large — a credible set concentrates mass on SNPs in near-perfect LD, so
`m_lᵀR m_l ≫ m_lᵀm_l` — and it is *negative*, offsetting the positive
between-effect term `μᵀRμ − μᵀμ`.

**My round-1 implementation omitted it**, computing `μᵀRμ + Σσ²`. That
over-counts by `Σ_l m_lᵀ(R−I)m_l`, which is precisely why it showed a spurious
*upward* bias (+0.031) at one causal SNP. The reviewer's concern found a real
bug, not a presentational gap.

**Now implemented and verified — and it reverses my round-1 answer.**
`R/susie_batch.R` computes the closed form above. Four checks, in
`tests/test_susie_estimand.R` and wired into pytest:

| check | result |
|---|---|
| `μ` equals `coef(fit)[-1]` | 1.7e−16 |
| `sd` equals `susie_get_posterior_sd(fit)` | 2.8e−17 |
| `h²_full` equals brute-force `tr(R·E[bb'])` from an explicitly built p×p matrix | 1.7e−16 |
| at `R = I`, `h²_full` collapses to `Σ(μ²+σ²)` | 0.0e+00 |

The `R = I` identity is additionally asserted on **every fit at run time**, so a
future change to the posterior parameterization cannot silently break it.

**The corrected result:** with the within-effect term included, the two forms
agree. Over the re-run baseline grid (3,000 paired fits) they correlate at
**r = 0.996**, the difference has **median 0.000**, and **82% of fits differ by
< 0.01**. The between-effect
cross terms (positive) are almost exactly cancelled by the within-credible-set
term (also positive) — which is why dropping only the second, as my round-1
version did, produced a large spurious `+0.065`.

| causal SNPs | bias, trace | bias, exact LD-aware |
|---:|---:|---:|
| 1 | −0.008 | −0.009 |
| 5 | −0.044 | −0.046 |
| 20 | −0.130 | −0.130 |
| overall | −0.058 | −0.059 |

**So the poster's central claim survives, and is now verified rather than
assumed:** the polygenicity shortfall is a property of the sparse-effect model,
not of the diagonal moment approximation. In the reviewer's terms this is the
first of their two candidate narratives, and §2 now says so as a settled result
rather than an open question. The dashed LD-aware curve is plotted on the §4
figure so a viewer can see the two estimators coincide.

### 2. Definitions box — **corrected**

The old box called MAE "precision", which is wrong: MAE mixes bias and
variability. Now:

> **Accuracy** — closeness to truth; bias and MAE, in simulation only.
> **Precision** — sampling variability; SD of ĥ² across replicate loci.
> **Agreement** — similarity of values; offset, slope, RMSE.
> **Reliability** — reproducibility; split-half correlation.

And rather than avoid the word, the poster now **shows a direct precision
analysis**, which the accepted title needs: SD of ĥ² across replicate loci is
**0.095 for CHEWS against 0.255 for GCTA at n = 100** (2.7×), narrowing to
0.044 vs 0.054 at n = 1,000. That is a genuine precision result and it is now
panel A of §5.

### 3. Literal `$h^2$` — **fixed, and guarded**

Real bug: the banner string carried TeX that PowerPoint printed verbatim. Now
`h²`. `build_poster.py` gained `assert_no_tex()`, which scans every string bound
for a PowerPoint text box and **fails the build** if `$…$` appears. `se` → `SE`
in the workflow strip.

### 4. Claims resting on unfinished work — **resolved**

The "open, and being measured now" box is gone; §2 now reports the finished
measurement. The factorial LD-source design moved to the methods page, so the
poster no longer advertises two unfinished analyses.

---

## Priority 1

### 5. Central plot title

"Only the sparse model loses signal as effects spread" is gone. The panel is now
titled **"Where the CHEWS shortfall comes from"**, and it shows all four curves —
GCTA, HEELS, CHEWS trace, CHEWS full — so the viewer sees which mechanism
carries how much.

### 6. Sample-size conclusion — **replaced with the interaction**

"CHEWS is most useful exactly where it is most biased" is gone. §5 panel B is now
the heatmap the review asked for: `MAE_CHEWS − MAE_GCTA` by *n* × causal count.
It shows the decision rule directly.

| MAE difference | n = 100 | n = 1,000 |
|---|---:|---:|
| 1 causal SNP | **−0.17** | −0.03 |
| 5 causal SNPs | −0.07 | −0.01 |
| 20 causal SNPs | **+0.04** | +0.03 |

Small *n* **and** concentrated effects — not small *n* alone. At 20 causal SNPs
CHEWS is worse at every sample size.

### 7. Simulation accounting and intervals

> 20,400 simulated replicates attempted; 20,294 complete method triplets.
> Unconstrained GCTA/GREML failed in 106; CHEWS and HEELS returned in all.

Intervals are now **region-clustered bootstraps** (2,000 resamples of whole
regions, since a region recurs across every causal count and h²), and the poster
says what is resampled. Honest note: at this level of aggregation clustering
changes the intervals very little (width ratio 0.69–1.08 versus the naive SE),
because each region contributes only six rows per causal-count level.

### 8. GTEx disclosure and agreement metrics

The box is now two lines. Agreement metrics are tied to their sample size and
gene set, with the offset added:

> At **n = 574 over 19,169 complete-case genes**: r = 0.9996, slope 0.93,
> offset −0.0074, RMSE 0.012.

### 9. Independence of the strata — **stated precisely**

Both stratifying variables — GCTA ĥ² *and* the clumped eQTL count — are computed
in half A; the difference is computed in half B. No donor is shared by the two.
The poster says this explicitly rather than leaving the eQTL axis ambiguous.

### 10. Donor-split variability — **limitation stated at the result**

> In **one** donor split CHEWS leads on all three gene sets; a gene-level paired
> bootstrap gives Δr = +0.050 (CI +0.038 to +0.063). That interval resamples
> **genes, not donor splits**, so it omits the variability that matters most here.

Repeating the split many times means re-running GCTA + CHEWS + HEELS on ~19k
genes per split; not done, and labelled as such rather than glossed.

### 11. Decision table

Header now reads "interpretation (current evidence)". The in-sample row is the
reviewer's wording: *HEELS targets the same REML equations*.

---

## Priority 2 — density and production

- Word count **1,012 → 915**, within the 800–900 band's reach; the §2 derivation
  shrank once the estimator was final. Wald coverage and the [0,1] screen moved
  into limitations. One sentence under each figure.
- Spelling standardized to **US**, matching the submitted abstract.
- `1000 Genomes` at first mention.
- A **human-readable URL** now sits under the QR code, so a failed scan is not a
  dead end.

**Still open, and flagged on the poster or here rather than hidden:**

1. Repeated donor splits.
2. The factorial model × LD-source design.
3. Validating the MAF result on discrete/haplotype genotypes.
4. Grant numbers and the real QR destination.
5. Rerunning the GTEx panels with the exact estimator — the simulation panels
   carry it, the GTEx pass predates it. Much less pressing now that the two forms
   agree at r = 0.996 on simulated data, but they are not formally on one
   estimator. ~4 h of compute.

## Final state

- 61 numbers re-derived from the result tables by `check_numbers.py`; all pass.
- 15 pipeline tests pass, including the new brute-force estimand check.
- Deck validated: 36 × 48 in, 81 shapes, nothing off-sheet, smallest type 20 pt,
  no literal TeX.
- Word count 915.
