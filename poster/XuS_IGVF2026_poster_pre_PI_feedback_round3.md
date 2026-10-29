# Pre-PI Review of the Current IGVF 2026 CHEWS Poster

## Verdict

This is now a credible and intellectually honest draft to show your PI. The poster has improved substantially in scientific rigor, visual hierarchy, and transparency. It renders cleanly at the full 36 × 48 inch size, the raw `$h^2$` problem is fixed, the statistical terms are much better defined, and the new exact quadratic-form comparison directly addresses the largest concern from the previous review.

Before sending it, I would make five targeted corrections. One sentence in Section 4 is factually inconsistent with its own figure, the exact estimator has not yet been propagated clearly through every analysis, the precision claim still needs a more defensible simulation design or wording, one symbol in the new formula is undefined, and the footer still contains visible placeholders.

## Five changes to make before sending this to your PI

### 1. Correct the erroneous sentence in Section 4

The current text says:

> GCTA, HEELS and the LD-aware estimator (dashed) all stay put, so the shortfall is the sparse-effect assumption itself.

The dashed purple LD-aware estimator does **not** stay put. It lies on the orange CHEWS curve and becomes increasingly negative as causal count increases. This sentence will immediately confuse a careful reader.

Replace it with:

> **With twenty causal SNPs, both the trace-form and exact LD-aware CHEWS estimators underestimate `h²` by approximately 0.13, while GCTA and HEELS remain near zero bias. Including the posterior LD cross terms therefore does not remove the polygenicity-associated shortfall in this analysis.**

Do not yet conclude “the sparse-effect assumption itself” unless the matched-LD factorial comparison also supports that attribution. CHEWS uses reference-panel LD while HEELS uses in-sample LD, so LD source remains a competing explanation.

### 2. Make the exact LD-aware estimator the primary CHEWS result everywhere

The exact calculation is now the scientifically preferred estimator. The poster says it was evaluated over 3,000 paired fits, which appears to correspond mainly to the `n = 300` architecture analysis. The sample-size panel and GTEx panels may still use the original trace form.

Your PI will likely ask: “Which version of CHEWS is represented in every figure?”

Resolve this by doing one of the following:

1. **Preferred:** use the exact LD-aware estimator in every main CHEWS result and show the trace form only as a sensitivity analysis.
2. If both forms are nearly identical across the full grid and GTEx, state that explicitly and give the maximum or 95th-percentile difference, not only the median.
3. If the exact estimator was tested only at `n = 300`, label that limitation and avoid generalizing its equivalence to the sample-size and GTEx analyses.

The claim “The two forms agree” is slightly too broad based on `r = 0.996`, median difference `0.000`, and 82% within `0.01`. Eighteen percent differ by at least `0.01`, which may be meaningful for local heritability. Add the 95th percentile and maximum absolute difference, stratified by causal count if possible.

A more precise sentence is:

> **At `n = 300`, the trace and exact LD-aware forms were similar on average (`r = 0.996`; median absolute difference = ...), and their mean bias curves overlapped across causal-count conditions.**

### 3. Define `s_l` in the exact estimator

The new formula defines

`m_l = alpha_l elementwise-times mu_l`

but does not define `s_l`. Every symbol in the key estimator must be interpretable from the poster.

If this matches your implementation, define it as:

`s_l = sum_j alpha_lj(mu_lj² + sigma_lj²)`

for unit-diagonal `R`. Then explain in one short phrase that the subtraction removes within-single-effect cross products that cannot occur because each SuSiE single effect selects only one variant.

I would rewrite the entire formula block as:

> **For `b = sum_l b_l`, let `m_l = E[b_l] = alpha_l elementwise-times mu_l` and `s_l = E[b_l'Rb_l]`. Then**  
> `E[b'Rb | data] = mu'Rmu - sum_l m_l'Rm_l + sum_l s_l`.

This is clearer than first assuming a diagonal posterior and then stating that the SuSiE posterior is not diagonal. The current sequence is mathematically correct as a comparison of an approximation with the exact form, but it is cognitively difficult to follow.

Also use either `b` or `beta` consistently throughout Section 2 and the brute-force check.

### 4. Reconsider what the poster calls “precision”

The definitions box is now conceptually correct. The remaining concern is the operational analysis.

Section 5 defines precision as the SD of estimated `h²` across “replicate loci,” but the benchmark says there is one replicate per region and parameter set. Different regions have different LD structure, so the observed SD combines:

- phenotype sampling variability;
- genotype or LD-region heterogeneity;
- any region-specific numerical behavior.

That quantity is useful, but it is not a clean estimate of conditional sampling precision for a fixed locus. Therefore, “CHEWS is 2.7× more precise than GCTA” is stronger than the current design supports.

Preferred analysis: generate repeated phenotypes for the same region and parameter condition, estimate within-region variance for each method, and then summarize that variance across regions.

If you cannot add repeated phenotypes before meeting with your PI, change the wording to:

> **Across simulated regions at `n = 100`, CHEWS estimates were 2.7× less dispersed than GCTA estimates (SD 0.095 versus 0.255). This dispersion includes both sampling variability and between-region heterogeneity.**

That phrasing is defensible and gives your PI a clear choice about whether a formal repeated-sampling analysis is needed for the manuscript.

### 5. Remove the visible footer placeholders

The footer still contains:

- `[repository URL here]`
- `[funding / grant numbers]`

These make the otherwise polished poster look unfinished. Before sending it to your PI:

- insert the repository URL if it is public;
- verify that the QR code resolves to that exact location;
- add a short readable URL below the QR code;
- replace the funding placeholder with the known grant information, or write “Funding information to be confirmed” in your email to the PI rather than on the poster.

If the repository is not yet public, remove the QR block temporarily or link to a stable project page that you control.

## Scientific points your PI is likely to raise

### Can the bias really be attributed to the sparse-effect assumption?

The exact-vs-trace comparison now shows that omitted posterior cross terms do not explain the mean causal-count trend at `n = 300`. That is useful evidence. However, CHEWS and HEELS still differ in both model and LD source. The causal statement becomes convincing only after a matched factorial comparison:

| Model | In-sample LD | Reference-panel LD |
|---|---:|---:|
| Sparse CHEWS | needed | current analysis |
| Infinitesimal HEELS | current analysis | needed |

If the factorial results already exist on the methods page, bring one decisive result onto the poster. If they do not exist, say “consistent with the sparse-effect model” rather than “the sparse-effect assumption itself.”

### Why does HEELS fail at `n = 100` if it targets the REML equations?

The current caveat correctly flags `p/n approximately 82`, but “breaks down” needs a measurable definition. Give the audience one number:

- failure or nonconvergence rate;
- mean bias;
- boundary rate;
- condition-number or LD-regularization behavior.

Clarify whether HEELS receives a rank-deficient empirical LD matrix, a regularized matrix, or a low-rank approximation. This is likely to be an important technical question at an IGVF meeting.

### Does the heatmap show probability of superiority or average error?

Section 4 reports the fraction of loci where CHEWS is closer to truth. Section 5 reports mean MAE differences. Both are useful, but they answer different questions. A method may win on more loci yet lose badly on a small number of loci.

Keep both metrics and state this distinction once. Consider adding median absolute error or a high quantile of absolute error to show whether rare catastrophic estimates drive the means.

### What uncertainty do the intervals represent?

The statement that intervals resample whole regions is good. Make the caption say “region-clustered bootstrap” if that is the procedure. Also state the number of bootstrap replicates on the methods page.

### What does the gene bootstrap in Section 8 support?

The new text correctly states that the interval resamples genes rather than donor splits. That is excellent. Do not visually emphasize the confidence interval more than the one-split limitation. The strongest accurate conclusion is:

> **CHEWS had higher split-half correlation in this donor split across all three gene sets; repeated donor partitions are required to quantify robustness to the split itself.**

## Visual and communication review

### What is working well

- The poster has a clear three-column reading order.
- The full-width conclusion and workflow make the topic understandable quickly.
- The method colors remain consistent.
- The revised definitions box is genuinely useful.
- “GTEx: agreement, not truth” is an excellent section title.
- The conditioning-artifact section demonstrates strong scientific self-correction.
- The decision table is cautious and concise.
- The main title, section headings, body text, and captions now have a workable size hierarchy.
- I found no visible clipping, overlap, or missing-glyph problem in the rendered poster.

### What can still be simplified

The poster contains approximately 986 words. That is an improvement, but it still feels like a dense methods paper at poster scale. A final cut of 100–150 words would help.

The easiest reductions are:

- condense the exact-estimator explanation after defining all symbols;
- shorten the updated-after-submission box to two sentences;
- remove the long HEELS caveat from the main column if the methods page contains the details;
- reduce the robustness paragraph to the single-split conclusion;
- move Wald-interval coverage and out-of-range behavior to the QR-linked supplement if space remains tight.

### Make panel B in Section 5 easier to read

The MAE heatmap communicates the interaction well, but the causal-count row values are difficult to see at normal poster distance. Make sure the y-axis shows explicit values `1, 3, 5, 10, 20` and enlarge them. Keep the colorbar label explicit:

`MAE_CHEWS − MAE_GCTA; negative values favor CHEWS`.

This will prevent readers from misreading which architectures reverse the advantage.

### Use more neutral language in two places

Replace:

> Where CHEWS beats GCTA

with:

> **Fraction of simulations with lower CHEWS absolute error**

Replace:

> Where the CHEWS shortfall comes from

with:

> **CHEWS bias by causal-effect count**

The neutral titles sound more like a mature statistical comparison and less like a competition between methods.

## Suggested final wording for the top banner

The current banner is good but can be slightly sharper:

> **Estimator performance depends jointly on regulatory architecture, sample size, and LD source. CHEWS lowers error in small samples with concentrated effects but underestimates `h²` as effects become diffuse.**

Use “performance” rather than “choice” because the current data compare statistical behavior; the poster does not yet validate every branch of a clinical or production decision rule.

## Suggested 60-second explanation to your PI

> I revised the poster around the bias–variance tradeoff rather than claiming CHEWS always outperforms GCTA. I corrected the GTEx covariate analysis, separated accuracy, precision, agreement, and reliability, and tested whether the diagonal posterior-moment approximation caused the polygenicity bias. The exact LD-aware posterior quadratic form nearly overlaps the trace form at `n = 300`, so the cross-term approximation does not explain the average trend there. The remaining limitation is that CHEWS and HEELS still use different LD sources, so I describe the sparse-model interpretation as provisional pending the factorial LD comparison. I also show that the previous GTEx stratification pattern largely disappears when strata are defined in independent donors.

That framing will help your PI see immediately what you changed, what is now defensible, and where you still want their guidance.

## Pre-send checklist

- [ ] Correct the Section 4 sentence about which estimators “stay put.”
- [ ] Define `s_l` and standardize `b` versus `beta` notation.
- [ ] State which CHEWS form appears in every figure.
- [ ] Report how the exact and trace forms differ in the tail, not only at the median.
- [ ] Replace “2.7× more precise” or add repeated phenotypes within fixed regions.
- [ ] Make the matched-LD limitation constrain the causal language about sparsity.
- [ ] Define what “HEELS breaks down” means numerically.
- [ ] Make Section 5 heatmap row labels legible.
- [ ] Replace competitive figure titles with neutral statistical titles.
- [ ] Remove repository and funding placeholders.
- [ ] Test the QR code and add a short URL.
- [ ] Cut approximately 100–150 additional words if time permits.

## Final assessment

You have done enough that sending the draft to your PI will demonstrate strong independent scientific judgment. The poster no longer looks like a preliminary visualization exercise; it now reads like a serious method comparison that recognizes its own inferential limits. Fix the Section 4 sentence, the undefined formula term, the precision wording, and the footer placeholders before sending. Then use your PI meeting to decide whether the matched-LD factorial experiment and within-region repeated simulations are essential for the conference poster or can be reserved for the manuscript.
