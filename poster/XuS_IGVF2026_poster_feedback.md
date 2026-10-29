# PI Review of the IGVF 2026 CHEWS Poster

## Bottom line

This is a strong and unusually thoughtful scientific story. The three-estimator comparison is a good idea, the paired simulations are appropriate, and the conclusion that performance depends on regulatory architecture is more credible than a simple claim that CHEWS always wins.

I would not print the poster in its current form. The main problems are scientific consistency and claim calibration, followed by excessive density. Most importantly, the GTEx numbers do not agree with the submitted abstract, the displayed CHEWS estimator needs a mathematically complete definition, and several conclusions are stronger than the analyses shown. The poster also contains about 1,237 words, with body text mostly around 19–23 pt and captions around 17 pt on a 36 × 48 inch poster. That is too much material at too small a size for normal conference viewing.

## What already works

- The central question matters to IGVF: can local heritability be estimated when individual-level genotypes cannot be shared?
- Including HEELS is conceptually valuable because it can help separate the effect of data representation from the effect-model assumption.
- The paired simulation design is stronger than comparing methods across independently generated phenotypes.
- The poster correctly recognizes a bias–variance tradeoff. CHEWS can have more bias but lower mean absolute error when REML is very noisy.
- The architecture-dependent result is the most interesting part of the study. Sparse, concentrated effects favor CHEWS, whereas diffuse effects produce downward bias.
- The poster appropriately states that GTEx has no ground truth. That distinction should remain central.
- Method colors are used consistently across most figures, and the numbered sections create a visible reading order.

## Priority 0: issues to resolve before printing

### 1. Reconcile the GTEx results with the submitted abstract

The poster and abstract report materially different real-data results.

| Quantity | Submitted abstract | Poster |
|---|---:|---:|
| Mean CHEWS − GCTA difference at `n = 100` | −0.151 | approximately −0.015 |
| Mean CHEWS − GCTA difference at the largest sample | −0.023 at `n = 574` | approximately −0.007 at `n = 574` |
| CHEWS–GCTA correlation across sample sizes | 0.61 to 0.88 | approximately 0.81 to 0.90 |

This is a tenfold discrepancy in the first mean difference, not a rounding issue. Determine whether the poster uses a revised analysis, different gene filters, complete-case genes, removal of out-of-range estimates, different covariate residualization, or a different LD panel. Then:

1. Use one analysis consistently throughout the poster.
2. State the gene inclusion rule and the number of genes contributing at every sample size.
3. If the analysis changed after abstract submission, label the poster result as an updated analysis and be prepared to explain the change in one sentence.
4. Recheck every numerical statement in the title, text, figure, and caption against a single frozen results table.

Also reconcile the counts of 19,432 analyzed genes, 19,169 genes in the HEELS comparison, and 4,171 genes in the split-half analysis. A small inclusion flow diagram would help.

### 2. Show the complete CHEWS estimand

The poster currently gives

`h_hat² = sum_j(mu_j² + sigma_j²)`.

For correlated standardized genotypes, genetic variance generally involves the LD matrix:

`E(beta' R beta | data) = tr[R E(beta beta' | data)]`.

The displayed diagonal sum is only sufficient if the effects have been transformed so that the relevant cross terms vanish, or if `mu_j` and `sigma_j²` have a special definition that already incorporates LD. This is the first equation a statistical-genetics audience will challenge because the simulation truth is explicitly defined as `beta' R beta`.

Add a compact derivation or exact estimator definition. Define:

- the scale of genotype and phenotype standardization;
- whether the posterior moments are for joint SNP effects, single-effect components, or a transformed basis;
- how LD and posterior cross-covariances enter the estimate;
- the denominator used to convert genetic variance to heritability;
- the SuSiE settings, including the maximum number of effects and prior-variance handling.

If the simple diagonal formula is correct under your construction, state the condition that makes it correct. Do not leave the audience to infer it.

### 3. The HEELS control does not currently isolate sparsity cleanly

Section 2 says HEELS shares CHEWS's inputs and GCTA's model, so the CHEWS–GCTA difference is attributable to sparsity alone. Section 12 says CHEWS used reference-panel LD while HEELS used in-sample LD. Those statements cannot both be true.

The clean experiment is a factorial comparison:

- CHEWS with in-sample LD;
- CHEWS with reference-panel LD;
- HEELS with in-sample LD;
- HEELS with the same reference-panel LD;
- GCTA/GREML on individual-level data.

This separates model assumption from LD-source mismatch. Until that experiment is shown, replace “attributable to the sparsity assumption alone” with a narrower statement.

There is a second inconsistency to explain. HEELS is described as solving the same REML equations as GCTA, yet in the simulation at `n = 100` HEELS shows large positive bias while GCTA does not. If this reflects rank deficiency, LD regularization, boundary constraints, or solver failure at high `p/n`, report it explicitly. HEELS cannot serve as a clean bridge between CHEWS and GCTA in conditions where it numerically breaks down.

### 4. Distinguish accuracy, precision, agreement, and reliability

The title says “precision,” but the simulation evidence emphasizes bias and mean absolute error, while the real-data analysis uses Pearson correlation across split samples.

- **Accuracy** means closeness to a known truth and can be assessed in simulation.
- **Precision** means sampling variability, such as variance or standard error.
- **Agreement** asks whether two methods give similar values, including calibration and scale.
- **Reliability or repeatability** asks whether estimates are stable across repeated samples.

Pearson correlation alone is not agreement and is insensitive to systematic offsets and rescaling. A high correlation does not show that HEELS and GCTA are “one method.” Show at least the mean difference, calibration slope and intercept, and RMSE or a Bland–Altman-style summary.

The statement that “a constant estimator scores a perfect ICC” should be removed or corrected. A constant estimator has no between-gene variance, so common ICC and correlation definitions are undefined or uninformative, not automatically perfect.

If the accepted title can still be changed, a more defensible title is:

> **CHEWS reduces cis-heritability error in small samples when regulatory effects are sparse**

If the title must remain identical to the accepted abstract, define “precision” operationally near the top and avoid using accuracy, precision, and reliability interchangeably elsewhere.

### 5. Do not infer sparse whole-blood architecture from method agreement alone

“Whole Blood cis-architecture is mostly sparse” and “the failure mode is real but rarely triggered in this tissue” are stronger conclusions than the GTEx comparison supports. A small CHEWS–GCTA difference can also reflect limited power, filtering, shrinkage, LD mismatch, or shared model behavior. There is no truth set in GTEx.

Use language such as:

> **The observed CHEWS–GCTA agreement in GTEx whole blood is consistent with predominantly concentrated cis effects under the CHEWS model, but does not independently establish sparsity.**

### 6. Rework the heatmap interpretation

Section 9 stratifies genes by GCTA `h²` and plots CHEWS − GCTA. Because GCTA appears in both the stratification variable and the outcome, mathematical coupling and regression to the mean can create or amplify a negative pattern. The independent cis-eQTL count also depends on power, sample size, and heritability.

Before calling this a Simpson's paradox or a biological signature of diffuse regulation:

- stratify on an independent quantity, such as GCTA from a disjoint donor subset, HEELS from independent data, or true `h²` in simulation;
- show the number of genes in every heatmap cell;
- repeat the analysis across multiple sample splits;
- model the continuous variables directly rather than relying only on bins;
- report uncertainty for the cell means.

At minimum, describe the current result as an observed conditional pattern, not proof that genes without detected cis-eQTLs have diffuse architecture.

### 7. The split-half result is vulnerable to selection bias

The robustness analysis keeps 4,171 genes with estimated `h² > 0.05` for every estimator in both halves. Conditioning on all six estimates can strongly distort correlations and favor methods differently. It also removes precisely the low-signal genes where the methods may behave most differently.

Show results for a prespecified gene set selected independently of the split-half estimates, or present both the unfiltered and filtered analyses. Repeat many random half-splits and report a distribution or confidence interval rather than one split. The statement that the method difference is no larger than donor sampling noise requires a direct paired comparison with uncertainty; comparing 0.835 with 0.848 is not sufficient by itself.

## Priority 1: strengthen the benchmark and recommendations

### Report failures, boundaries, and exclusions

The total of 20,294 fits is not transparently connected to 350 regions and 13 experimental conditions. State the number attempted, successful, nonconverged, boundary-constrained, and excluded for each method. Specify whether GCTA used constrained or unconstrained REML and whether CHEWS and HEELS were treated equivalently. At `n = 100` with a median of 8,200 SNPs, numerical behavior is part of the scientific result, not a minor implementation detail.

### Do not silently remove out-of-range estimates

The abstract recommends removing CHEWS estimates outside `[0,1]`. Simply dropping them can create selection bias and hide numerical failures. Predefine how these values are handled, retain raw results for method evaluation, report them as failures or boundary violations, and give a sensitivity analysis. If values reach 3–5 despite clean convergence, investigate LD mismatch, scaling, and posterior-variance calculation before presenting this as a routine screening step.

### Narrow the standard-error claim

“Any reported `h²` standard error: do not trust” is too broad. Replace it with something like:

> **Nominal 95% Wald intervals under-covered in these simulations (coverage 68–87%, depending on method and condition). Report the method-specific coverage and consider bootstrap or profile-likelihood intervals.**

Show coverage by estimator, sample size, true `h²`, and architecture. Boundary behavior often drives poor Wald coverage, so pooling hides the mechanism.

### Validate the simulation architecture

The multivariate-normal genotype simulation preserves pairwise LD but not discrete genotype counts, Hardy–Weinberg structure, or haplotypes. This is especially important for the MAF-stratified result and the recommendation to use an MAF-stratified GRM. Validate that result with sampled or haplotype-based genotypes before making it a headline recommendation.

Other simulation details that should be visible or accessible through the QR code:

- how the 350 regions were sampled and whether regions were reused across conditions;
- the number of phenotype replicates per condition;
- causal-variant selection and effect-size scaling;
- the reference-panel size and ancestry match;
- SuSiE and HEELS tuning or regularization;
- clustered uncertainty by region when the same region appears repeatedly;
- the practical meaning of “polygenic,” since 20 causal variants among approximately 8,200 SNPs remains sparse in a literal sense.

### Make the decision guide conditional

The current recommendations are too categorical. In particular, “small study (`n <= 300`) use CHEWS” conflicts with the finding that CHEWS performs poorly for diffuse architectures. A more defensible guide is:

| Data and architecture | Recommended interpretation |
|---|---|
| Summary statistics plus exact in-sample LD | HEELS can approximate the REML target under the same infinitesimal model; verify numerical agreement. |
| Reference-panel LD only | Assess LD match and run sensitivity analyses; the current benchmark does not isolate LD source from estimator choice. |
| Small `n` and concentrated signals | CHEWS may reduce MAE despite downward bias; report both bias-sensitive and variance-sensitive metrics. |
| Diffuse or weakly localizable cis signal | CHEWS may underestimate `h²`; prefer an infinitesimal model or present both estimates. |
| Strong MAF dependence | Use an appropriately stratified model only after validating the result with discrete or haplotype-based genotypes. |

## Priority 2: redesign the poster around one message

### Reduce the poster from 12 sections to about 5 modules

The current poster reads like a compressed manuscript. A viewer will not absorb 12 sections and more than 1,200 words while standing at the board. Aim to cut the prose by roughly 40–50% and make the three take-home messages visible within 20–30 seconds:

1. With exact in-sample LD, summary statistics can closely reproduce the REML target under the same model.
2. CHEWS trades variance for downward bias and lowers error mainly when sample size is small and effects are concentrated.
3. Diffuse effects and LD-source mismatch are the main reasons to distrust a fine-mapping-based `h²` estimate.

### Suggested three-column structure

**Full-width top band**

- Shorter title if permitted.
- One-sentence conclusion.
- A simple CHEWS workflow: summary statistics + LD → SuSiE posterior → posterior genetic variance → cis-`h²`.

**Column 1: Question and benchmark**

- Why summary-data local heritability matters.
- Condensed three-estimator comparison.
- Simulation design with a small schematic rather than paragraphs.

**Column 2: The main simulation result**

- Make the architecture-by-`h²` heatmap the dominant figure.
- Add the sample-size bias/MAE panel below it.
- Use one interpretation sentence: CHEWS lowers error in small samples only when signal is sufficiently concentrated.

**Column 3: GTEx and practical interpretation**

- GTEx agreement plot with exact sample counts and uncertainty.
- One cautious conditional analysis of disagreement.
- A revised decision guide and two limitations.

**Footer**

- Three to five essential references.
- Filled funding and acknowledgments.
- QR code to methods, full result tables, code, and additional sensitivity analyses.

Move the current MAF/LD sensitivity panels, detailed robustness analysis, full references, and secondary recommendations to the QR-linked supplement unless one becomes essential after the scientific revisions.

### Increase type size and simplify the visual hierarchy

Approximate current sizes are 52 pt for the title, 27 pt for section headers, 19–23 pt for most body text, and 17 pt for several captions. For a 36 × 48 inch poster, target approximately:

- title: 72–90 pt;
- authors: 30–36 pt;
- section headings: 32–38 pt;
- body: 24–28 pt;
- captions and references: at least 20–22 pt.

Test by printing one dense figure and its caption at 100% scale on letter paper, then view it from several feet away. If increasing the font causes overflow, remove text rather than shrinking it.

### Improve figure communication

- Give each figure a declarative takeaway title rather than labels such as “The one place GCTA breaks.” A neutral title would be “MAF restriction shifts GCTA bias upward in this simulation.”
- Direct-label lines when possible and reduce repeated legends.
- Use point shapes or line styles in addition to color.
- Keep the CHEWS orange, GCTA blue, and HEELS green mapping everywhere. Align the heatmap palette and annotations with that logic or explain the diverging scale clearly.
- Increase axis, tick, annotation, and panel-letter sizes. Several plot annotations are currently unreadable at normal poster distance.
- Put sample size, number of regions or genes, and uncertainty definition in every figure caption.
- Define “95% Monte Carlo intervals” and account for repeated use of the same cis-regions, ideally with a region-clustered bootstrap.

### Reduce argumentative phrasing

The poster is commendably opinionated, but several phrases sound more definitive than the data:

- “HEELS reproduces REML exactly” should become “HEELS closely agrees with GCTA/GREML when both use the same donors and in-sample LD,” followed by agreement metrics.
- “Nothing is lost by sharing summary statistics” should become “Under the same linear mixed model and with exact in-sample sufficient statistics, summary-data estimation can closely recapitulate REML.”
- “GCTA is not an unconditional gold standard” is rhetorically unnecessary. State the specific model misspecification and observed bias instead.
- “The worst cell is the gene...” should become “The largest mean difference occurs among genes with...” because each cell contains a group of genes.
- Replace all-caps emphasis such as `SAME`, `OPPOSITE`, `SUFFICIENT`, and `DISAGREEMENT` with restrained bold text.

## A concise oral story for the poster session

Use a 60–90 second opening rather than walking through all numbered sections:

> CHEWS estimates local heritability from fine-mapping posteriors, so it can be used when only summary statistics and LD are available. We compared it with individual-level GCTA and summary-data HEELS. The key result is a bias–variance tradeoff: in small samples with concentrated causal effects, CHEWS has lower error than GCTA even though it is downward biased. As effects become diffuse, CHEWS increasingly misses genetic variance. In GTEx whole blood the methods are broadly concordant, but those data do not provide ground truth, so we use them as an agreement analysis rather than a validation of accuracy. The practical message is that estimator choice should depend on LD availability, sample size, and how concentrated the cis signal appears to be.

Then ask the visitor whether they are most interested in the estimator, the simulation architecture, or the GTEx application and take them to the relevant figure.

## Final pre-print checklist

- [ ] Poster and abstract use reconciled GTEx values.
- [ ] Every reported number comes from one frozen results table.
- [ ] The CHEWS equation includes LD or states why cross terms vanish.
- [ ] CHEWS and HEELS are compared under matched LD inputs.
- [ ] Attempted fits, failures, boundaries, and exclusions are reported.
- [ ] Accuracy, precision, agreement, and reliability are used correctly.
- [ ] Heatmap conditioning and split-half selection bias are addressed.
- [ ] Claims about sparsity, exact equivalence, standard errors, and out-of-range values are softened or supported.
- [ ] The recommendation table is conditional on both sample size and architecture.
- [ ] Body text is at least about 24 pt and the poster has roughly half as much prose.
- [ ] Funding, grant numbers, repository URL, and QR code replace the placeholders.
- [ ] A final 100% PDF export is checked for missing fonts, clipped equations, and readable figures.

## Method references relevant to these revisions

- Li H, Mazumder R, Lin X. [Accurate and efficient estimation of local heritability using summary statistics and the linkage disequilibrium matrix](https://www.nature.com/articles/s41467-023-43565-9). *Nature Communications*. The paper emphasizes the in-sample LD requirement for HEELS and cautions that out-of-sample LD can introduce bias.
- Zou Y et al. [Fine-mapping from summary data with the Sum of Single Effects model](https://doi.org/10.1371/journal.pgen.1010299). *PLOS Genetics*. The paper distinguishes exact sufficient statistics from the approximation that uses marginal statistics and reference-panel LD.
- Yang Lab. [GCTA-GREML guidance for small samples](https://yanglab.westlake.edu.cn/software/gcta/greml_ldms.html). The guidance notes that small-sample GREML estimates have large sampling variance and can hit parameter boundaries.
