# Round 4: update for Tiffany

## Email draft

**Subject:** Poster update: when CHEWS improves estimates, and what we still need to check

Hi Tiffany,

Thank you for your feedback. I updated the poster to show more directly when CHEWS gives a better estimate than GCTA. I also ran two follow-up analyses to address your question about genes with very low estimated heritability.

The clearest finding is that **CHEWS performs well when a few variants account for most of a gene's genetic effect. When the effect is spread across many variants, GCTA often performs better.** The new analyses give us some evidence about the low-heritability genes, but that part of the answer is still preliminary.

**1. The revised figure shows how often CHEWS is closer to the true value.**

In the simulations, we know the true heritability, so we can compare the two methods directly. I replaced the figure showing their difference in average error with one showing the percentage of simulated cases in which CHEWS has the smaller error.

- With one causal variant, CHEWS is closer to the truth in **75–82% of cases** across sample sizes from 100 to 1,000. Its average error is 73% lower than GCTA's at 100 individuals and 58% lower at 1,000 individuals.
- With 20 causal variants, **GCTA is closer in 72–77% of cases** at sample sizes of 300 or more. This result is also shown on the poster so that the conditions favoring each method are clear.

On your question about percentages: each part of this figure combines several true heritability values, from 0.05 to 0.8. Expressing every error as a percentage of the true value would give more weight to cases with low heritability. The revised figure answers “how often does CHEWS improve the estimate?” directly. The percentage reduction in average error answers the separate question, “how much does it improve?”

I kept the standard-error comparison and the split-half consistency table as they were, as you suggested.

**2. For the low-heritability genes, lower CHEWS estimates may be helpful, but we cannot yet make that claim for the real data.**

I focused on genes with estimated cis-heritability below 0.05 and no detected independent association signal. This group contains 9,572 genes in the real-data comparison. “No detected signal” does not necessarily mean that a gene has no genetic contribution; weak effects may be difficult to detect.

Because we do not know the true heritability of these genes, I ran two checks:

- **A test with the genotype–expression relationship removed.** I shuffled the donor labels on the genotypes while keeping expression and covariates unchanged. In this no-signal setting, GCTA and CHEWS both gave an average heritability estimate of about **0.013**. This does not support the simple explanation that GCTA's restriction against negative estimates, by itself, explains why its estimates are higher than CHEWS's in the real low-heritability group.
- **A simulation using the actual GTEx genotypes.** I generated expression values with known genetic contributions, using the same donors, donor split, and adjustment variables as the real analysis. In the simulated group corresponding to those low-heritability genes, CHEWS was closer to the truth in about **60% of cases** (95% confidence interval: 59–61%). When I varied the assumed proportion of genes with no genetic contribution from 0% to 90%, that figure ranged from **54% to 61%**.

The second result is encouraging, but the simulation did not reproduce the overall real-data pattern closely enough to meet the check set before examining its results. Across the gene groups, the average difference between the simulated and observed CHEWS–GCTA gaps was **0.014**, above the allowed **0.010**. Most of the mismatch was among genes with high estimated heritability. The confidence interval around the 60% result does not account for this mismatch between the simulation and real data.

For now, the poster includes the shuffled-data result, and the simulated accuracy result remains exploratory. My answer to “is the lower estimate good or bad?” is: **it looks modestly helpful in the simulation, but we do not yet have enough evidence to say that it is more accurate for these real genes.**

One possible reason for the mismatch is that the simulation only included true heritability values up to 0.4. Extending that range to 0.8 is a useful next check, although it may not resolve the problem. The previous run took about 20 minutes of computing time; the new run may differ. We would keep the same evaluation criteria and clearly report that the simulation design was expanded after reviewing the first result.

**3. I have also drafted an outline for the paper.**

I suggest organizing the paper around when CHEWS improves estimation, when it misses genetic effects, and what the GTEx results can tell us about its practical use. The outline separates the completed findings from the analyses we still need before making stronger claims.

I would appreciate your thoughts on three points:

1. When you referred to the “bottom-left” group, did you mean genes with no detected signal and estimated heritability below 0.05? I checked the other groups below 0.05 as well.
2. Should I run the wider-range simulation before we print the poster, or keep that as follow-up work for the paper?
3. Does this framing work for the paper, and which of the remaining analyses would you prioritize?

The detailed response is in `poster/REVISION4.md`, and the updated poster is `XuS_IGVF2026_poster_round4.pptx`.

Best,
Sunan

---

## Short message alternative

Hi Tiffany — I updated the figure to show how often CHEWS is closer to the true heritability in the simulations. With one causal variant, CHEWS is closer in 75–82% of cases across sample sizes of 100–1,000. With 20 causal variants, GCTA is closer in 72–77% of cases once the sample size reaches 300. I kept the standard-error comparison and split-half table unchanged.

For your question about the low-heritability genes: after shuffling the genotype–expression relationship, GCTA and CHEWS both average about 0.013. A separate simulation using real GTEx genotypes suggests that CHEWS's lower estimates may be modestly more accurate in the low-heritability, no-detected-signal group. However, that simulation did not reproduce the overall real-data pattern well enough to pass our planned check, so I have kept its accuracy results off the poster. Extending the simulated heritability range is a reasonable next step, but we do not yet know whether it will resolve the mismatch.

Would you prefer that follow-up before printing, or as part of the paper? I also drafted a paper outline and would appreciate your thoughts on the main story and the remaining analyses.
