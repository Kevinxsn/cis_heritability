# Second PI Review of the IGVF 2026 CHEWS Poster

## Overall assessment

This version is substantially better. It is more transparent, more statistically literate, and easier to navigate. You addressed nearly every major concern from the first review: the revised GTEx analysis is disclosed, accuracy and agreement are discussed separately, the HEELS comparison is no longer presented as a clean sparsity control, the conditioning artifact is tested with independent donors, selection sensitivity is shown, and the recommendation table is appropriately conditional.

The poster is now close to being a strong conference poster, but I would still not print this exact version. One unresolved analysis affects the central scientific conclusion: the poster says the omitted LD-weighted posterior cross terms may explain part of the apparent loss of heritability as effects become more polygenic. If that is true, the current results do not yet distinguish a limitation of the sparse SuSiE model from a limitation of the diagonal or trace approximation used to convert the posterior into `h²`.

My recommendation is to finish the full quadratic-form rerun before finalizing the poster. Most remaining revisions are then straightforward.

## Major improvements since the previous version

- The GTEx discrepancy with the abstract is now acknowledged and explained.
- The poster distinguishes simulation truth from real-data agreement.
- The incomplete equivalence between HEELS and GCTA is now demonstrated with slope and RMSE, rather than correlation alone.
- The confounding of model choice and LD source is stated directly.
- The same-donor conditioning pattern is compared with an independent-donor analysis, and the Simpson's paradox claim has been withdrawn.
- Split-half results are shown under multiple gene-selection rules.
- Attempted analyses and GCTA failures are reported.
- The recommendations no longer claim one universally superior method.
- The poster has fallen from approximately 1,237 to 1,012 words.
- Typography has improved substantially: the title is now approximately 72 pt, section headings are mostly 32–33 pt, body text is generally 22–25 pt, and captions are approximately 20–21 pt.
- The top workflow gives the viewer an immediate conceptual map.
- The QR code and concise reference list are much better suited to a poster than the previous full bibliography.

## Priority 0: complete before printing

### 1. Resolve the quadratic-form issue and rerun the key analyses

The revised mathematical explanation is a major improvement, but it exposes a potentially central problem.

You correctly state that

`E(beta' R beta | data) = mu' R mu + tr(R Sigma)`.

The implemented estimator uses only

`sum_j(mu_j² + sigma_j²)`.

The poster then says that the difference from the full expression is “exactly `mu'Rmu − mu'mu`” under diagonal posterior covariance. This is only exact if `Sigma` is genuinely diagonal. A SuSiE posterior can induce dependence among SNP effects, especially among alternative variants within a single effect. If off-diagonal posterior covariances are present, the omitted contribution also includes LD-weighted covariance terms:

`sum_{j != k} R_jk Sigma_kj`.

Therefore, confirm rather than assume that the relevant posterior covariance is diagonal. The full computation should ideally use posterior draws or the complete second-moment matrix:

`tr[R E(beta beta' | data)]`.

After implementing the full expression, rerun at least:

1. bias against causal count;
2. the architecture-by-`h²` heatmap;
3. bias and MAE against sample size;
4. the GTEx offset and calibration analyses.

Then compare the diagonal and full estimators directly. This comparison may become the most interesting result in the poster.

Until that rerun is complete, the following claims are not supported cleanly:

- “Only the sparse model loses signal as effects spread.”
- “The shortfall is specific to the sparse-effect estimator.”
- “CHEWS underestimates `h²` for diffuse signal.”
- Any recommendation based on polygenicity.

The current data could instead show that the diagonal posterior-moment approximation loses positive LD cross terms as more effects enter the model.

**Preferred path:** finish the full-estimator rerun and update the conclusions.

**Fallback if timing makes that impossible:** rename the method shown here as the “trace-form CHEWS estimator,” replace causal claims about sparsity with descriptive claims about the current implementation, and frame the full quadratic estimator as necessary future work. This would be honest, but scientifically less satisfying.

### 2. Correct the definitions box

The current box says:

> Precision = sampling variability: MAE and split-half spread.

This is still incorrect. MAE measures overall error relative to truth and combines bias with variability. Split-half Pearson correlation measures reliability across genes, not “spread.”

Use:

> **Accuracy** = closeness to truth; assessed with bias and MAE in simulation.  
> **Precision** = sampling variability; assessed with the SD or SE across repeated samples.  
> **Agreement** = similarity of numerical values; assessed with offset, slope, and RMSE.  
> **Reliability** = reproducibility across samples; assessed here with split-half correlation.

If you do not actually show an estimator-level SD or SE comparison, avoid claiming improved “precision” in the poster body. Because the accepted title still uses “precision,” either add a direct precision analysis or state that the title uses the submitted terminology while the poster distinguishes lower MAE from higher repeatability.

### 3. Remove the raw LaTeX from the top message

The summary banner visibly contains `$h^2$`. PowerPoint is displaying the markup literally. Replace it with the formatted character `h²` or a native PowerPoint equation.

Also change `se` in the workflow to `SE` if it denotes standard error.

### 4. Do not build the central conclusion around analyses marked “in progress”

The poster currently contains both:

- “Open, and being measured now” for the full quadratic estimator;
- “The factorial run that separates these is in progress” for the model-by-LD comparison.

This degree of candor is admirable, but it tells the audience that two analyses required to interpret the main result are unfinished. Complete them if at all possible. If they cannot be completed, reduce the claims and present this explicitly as preliminary method development.

## Priority 1: statistical and interpretive refinements

### 5. Replace the central architecture headline

The title above the causal-count plot currently says:

> Only the sparse model loses signal as effects spread

That conclusion is contradicted by the caveat immediately below it. Before the quadratic rerun, use:

> **Trace-form CHEWS increasingly underestimates `h²` as effects spread**

After the rerun, the title should name whichever mechanism survives the comparison. For example:

- **Diffuse effects drive CHEWS bias even after including LD cross terms**, or
- **Omitted LD cross terms, rather than sparsity alone, explain most downward bias**.

### 6. Tighten the sample-size conclusion

“CHEWS is most useful exactly where it is most biased” is memorable but too absolute. The plotted MAE is pooled over causal count and true `h²`, so it does not establish usefulness for every architecture at small `n`.

Use:

> **Across the pooled simulation grid, lower sampling error offsets CHEWS's downward bias most strongly at `n = 100`.**

More importantly, show the interaction between sample size and architecture. A small two-dimensional panel or supplementary heatmap of `MAE_CHEWS − MAE_GCTA` by `n` and causal count would support the actual decision rule: small sample **and** concentrated effects.

### 7. Clarify the simulation accounting

“20,400 cells attempted” is confusing. A cell usually means one parameter combination, not one simulated dataset.

Use something like:

> **20,400 simulation replicates attempted; 20,294 complete method triplets. Unconstrained GCTA/GREML failed in 106 replicates; CHEWS and HEELS returned estimates in all replicates.**

If that is not the exact unit, state the correct unit explicitly. Also clarify whether the reported Monte Carlo intervals are ordinary bootstrap intervals, standard-error intervals, or region-clustered bootstrap intervals. Because regions recur across the grid, cluster resampling by region is preferable.

### 8. Make the updated GTEx analysis more concise and complete

The disclosure is scientifically appropriate, but it currently occupies a large box and reports only the old mean difference. A shorter version would read:

> **Updated after abstract submission:** covariates are now included directly in REML rather than used to pre-residualize expression. The values shown here replace the submitted GTEx estimates.

Keep the detailed degrees-of-freedom explanation for conversation or the QR-linked methods page unless you have explicitly verified the direction and magnitude of that bias.

For the HEELS statement, specify the sample size and gene set attached to `r = 0.9996`, slope `0.93`, and RMSE `0.012`. If these are values at `n = 574`, say so. Add the mean offset or intercept at that same sample size so agreement is fully characterized.

### 9. Preserve the conditioning-artifact analysis, but define independence precisely

Section 7 is one of the strongest changes in the poster. Keep it.

Clarify whether both stratification variables, GCTA `h²` and independent cis-eQTL count, come from the disjoint donor half. If only GCTA `h²` is estimated in the other half, the eQTL-count axis can still carry same-sample selection effects. State exactly which quantities are computed in each half.

The heading “A conditioning artifact” is accurate. The text can be shortened to:

> **Using GCTA in both the stratum definition and the CHEWS–GCTA difference creates mathematical coupling. Defining strata in an independent donor half removes most of the apparent pattern.**

### 10. Repeat the donor split before making a precision claim

The multi-filter table is useful and shows that the result is not solely caused by the post-hoc `h² > 0.05` filter. However, all three rows still come from one split of the same 574 donors.

The gene bootstrap confidence interval estimates uncertainty across genes conditional on that split. It does not capture variation caused by which donors happen to enter each half. Repeat the donor split many times and report the distribution of `Delta r = r_CHEWS − r_GCTA`.

Until then, use:

> **In one donor split, CHEWS had higher split-half correlation across all three gene sets. A gene-level paired bootstrap gave `Delta r = 0.050` on the prespecified eQTL set, but this interval does not include donor-split variability.**

Also replace “all genes (18.5–19.4k)” with exact method-specific or common-intersection sample sizes. A correlation comparison should make clear whether all methods use the same genes.

### 11. Keep the decision table explicitly provisional

The table is much improved, but the polygenicity-dependent rows rely on the unresolved estimator calculation. Add a small label such as “Current evidence” and update it after the quadratic rerun.

For exact in-sample LD, a more precise row is:

| Situation | Interpretation |
|---|---|
| In-sample sufficient statistics and LD | HEELS targets the same infinitesimal-model REML equations; verify numerical agreement and calibration. |

For reference-panel LD, the current “not isolated here” wording is appropriately cautious.

## Priority 2: poster communication and visual polish

### 12. The layout is now workable, but the poster remains slightly too dense

The poster is visibly easier to read than the first version, and the current font sizes are close to acceptable. The remaining problem is cognitive density rather than simple font size. Approximately 1,012 words plus six substantial figures and multiple caveat boxes still asks too much of a standing viewer.

Aim for approximately 800–900 words. The easiest cuts are:

- shorten the updated-analysis box;
- move details of the unfinished factorial design to the QR-linked methods page;
- reduce the Section 2 derivation after the estimator is finalized;
- shorten the prose below the sample-size plot;
- keep only one sentence below each heatmap or plot;
- move the Wald-coverage and out-of-range details into the limitations section.

### 13. Preserve the new hierarchy

The following design choices now work well and should remain:

- the 72 pt two-line title;
- the full-width conclusion banner;
- the input-to-estimate workflow;
- three clearly separated columns;
- consistent blue, orange, and green method colors;
- the architecture heatmap as the main left-column result;
- the conditional recommendation table;
- the compact references and QR code.

The central orange caveat boxes are visually prominent. Once the analyses are complete, replace them with one concise “What changes the interpretation?” box or remove them so the principal results regain visual priority.

### 14. Final production details

- Replace `[funding / grant numbers]` before export.
- Test the QR code from the printed poster size with at least two phones and include a short human-readable URL below it.
- Confirm that the QR destination is public and does not require institutional authentication.
- Use `1000 Genomes` rather than `1000G` at first mention.
- Standardize `cis-h²`, `h²`, `h_hat²`, and `h^2` notation.
- Use one spelling convention throughout. The current version mixes US and UK forms in places.
- Ensure every reported interval states what was resampled and whether the same regions or genes recur.
- Export the final PDF at the full 36 × 48 inch size and inspect it at 100% for font substitution, clipped superscripts, and thin plot lines.

## Recommended final narrative after the rerun

If the full quadratic estimator still shows the same qualitative behavior, the poster can make this clean claim:

> **CHEWS estimates cis-heritability from summary statistics by integrating posterior genetic variance across fine-mapped effects. In simulations, it lowers error relative to GCTA/GREML in small samples when regulatory effects are concentrated, but becomes downward biased as effects become diffuse. With exact in-sample LD, HEELS closely recapitulates the infinitesimal REML target. GTEx whole-blood results show method agreement rather than accuracy and are consistent with, but do not prove, concentrated cis architecture.**

If the full quadratic estimator substantially repairs the polygenicity bias, the more interesting narrative is:

> **Posterior cross terms are essential for estimating cis-heritability from fine-mapping output. A diagonal moment approximation loses genetic variance as the number of correlated effects increases, while the full LD-aware posterior quadratic form reduces that bias.**

The second outcome would change the paper, but it may produce a more novel and methodologically important result.

## Final decision checklist

- [ ] Full posterior quadratic form implemented and compared with the diagonal approximation.
- [ ] Off-diagonal posterior covariance assessed rather than assumed absent.
- [ ] Main simulation and GTEx panels rerun with the final estimator.
- [ ] Accuracy, precision, agreement, and reliability box corrected.
- [ ] Literal `$h^2$` replaced with formatted `h²`.
- [ ] Central plot title no longer attributes the effect to sparsity prematurely.
- [ ] Simulation unit and interval construction clarified.
- [ ] HEELS agreement metrics tied to a specific sample size and gene set.
- [ ] Independent-half construction fully defined.
- [ ] Donor split repeated or single-split limitation stated next to the result.
- [ ] Word count reduced toward 800–900.
- [ ] Funding information completed and QR destination tested.

## PI verdict

The revision shows excellent scientific judgment and responsiveness. The poster is clearer, more cautious, and more credible than the first version. The unresolved quadratic-form calculation is now the decisive issue. Resolve that first, because it may change both the biological interpretation and the main methodological contribution. Once the estimator is frozen and the definitions box is corrected, the remaining work is primarily compression and production polish.
