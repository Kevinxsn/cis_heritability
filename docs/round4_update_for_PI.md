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

---

## Proposed next steps to discuss with Tiffany

*These are proposed analyses, not additional results. They are intended to strengthen the paper's main claims before we expand its scope.*

I think we already have a useful story: CHEWS often improves estimation when genetic effects are concentrated, but misses heritability when effects are more distributed. The most useful next work would tell us **how researchers can recognize when CHEWS is useful, why its performance changes, and whether that guidance survives a fairer comparison.** I would prioritize the following.

### 1. Main priority: can observable features help researchers choose a method?

**The question:** Before knowing a gene's true heritability, can we use information available in a real analysis to decide whether CHEWS or another estimator is likely to give a more accurate result?

**Why this could be the main contribution:** Our simulations show that CHEWS often helps when genetic effects are concentrated. But researchers do not know the true number or distribution of causal effects in their own data. If we can recognize favorable conditions from observable features, we could turn the benchmark into practical guidance: **which method to use, for which genes, and when the evidence is too uncertain to recommend one.**

**A manageable first analysis:** Use the existing simulation results to test whether a few observable features predict the difference in absolute error between CHEWS and GCTA. This tells us both which method is closer and how much is gained or lost by choosing it. The true heritability is used to evaluate the decision in simulations, never as an input to the decision rule.

Start with features already saved in the result tables:

| Observable feature | Why it might help |
|---|---|
| Sample size and number of variants analyzed | Indicate how much information is available relative to the size of the region |
| Strength of the strongest association | Helps distinguish a clearly detectable signal from a weak one |
| Highest fine-mapping probability for a variant | Measures how strongly the model supports a particular variant as contributing to the signal |
| Number and size of SuSiE credible sets—groups of plausible causal variants | May distinguish signals that are easier to locate from signals whose location is uncertain |

These are hypotheses to test. In particular, confidence about a variant's location does not prove that it explains most of the genetic contribution. The association-based signal count already available in the GTEx analyses could be an additional feature once it is calculated consistently in the main simulations. “No detected signal” must not be treated as “no genetic effect.”

**Keep the first rule simple.** Begin with plots of error differences across these features, then test a small decision tree or a few clearly stated rules using two or three features. For example, we could ask whether sample size combined with strong support for a small number of signals is informative. We should learn and validate the cutoffs rather than assume that this pattern favors CHEWS. The number of independent regions is much smaller than the number of simulated fits, which argues for a simple first model.

**Test it on genuinely separate regions.** Keep every simulation from the same or overlapping genomic regions together when dividing data into development and evaluation sets, including repeats across experiments. Choose features and cutoffs using the development data, then freeze the rule before evaluating the held-out regions. Otherwise, the rule could appear useful because it has already seen almost the same region during development.

**Define success by better estimates, not just correct winner predictions.** On held-out regions, compare the rule with always using CHEWS and always using GCTA. The main question is whether it reduces average absolute error relative to the better fixed choice. Also report the size of losses when it chooses poorly, performance for small but nonzero heritability, and uncertainty from resampling whole regions. Report results by simulation setting and under different mixtures of settings, so success does not depend on an arbitrary balance of easy and difficult cases. A rule that predicts many tiny wins but makes a few costly mistakes would not be useful.

The available data must also determine the choices: GCTA requires individual-level data, while the current HEELS implementation requires in-sample LD. Evaluate a separate CHEWS-versus-HEELS rule where those inputs are available. A rule using SuSiE outputs requires running fine-mapping first; its initial purpose is better estimation, not saving computing time. Failed fits should be flagged explicitly rather than counted as ordinary choices.

**Allow an uncertain answer.** For genes where the evidence does not clearly favor a method, the guide could recommend reporting both estimates and their difference. Evaluate how often the rule can make a useful recommendation, and the error among those recommendations, alongside its overall performance. The goal is a trustworthy guide, not a forced choice for every gene.

**What we could show in the paper:** A simple decision diagram, a figure showing where each method has lower error on held-out regions, and a comparison of the rule's error with always choosing one method. Applying the frozen rule to GTEx could show how often each recommendation occurs and whether recommendations are stable across donor splits. It would not establish real-data accuracy, because true heritability remains unknown there.

**First milestone:** Determine whether observable features provide a reproducible improvement over the better fixed method in the existing simulations. If they do, validate the rule in the more controlled comparisons in proposal 2 and with simulations on real genotypes. This priority can start now; those follow-up checks are needed before presenting it as guidance for real studies. If no robust improvement emerges, we still learn that the current features are insufficient for reliable method selection.

### 2. Test what is causing the underestimation

**Why this matters:** We currently interpret the shortfall at twenty causal variants as consistent with the fine-mapping model. However, three features of the experiment could contribute: CHEWS uses reference-panel LD, SuSiE allows only ten effect components, and the two summary-statistic methods exclude some variants that remain available to GCTA. Increasing the effect limit might help, but could also increase variability or unstable fits. We should test the trade-off rather than assume that a larger limit is better.

**What I propose:** Start with a focused experiment covering the cases where CHEWS clearly helps and clearly struggles: one versus twenty causal variants, small versus moderate sample sizes, and low versus higher true heritability. Save the simulated genotypes and expression once, so every comparison uses exactly the same data, including when a fit fails.

- Give all methods the same variant set, retaining strand-ambiguous variants where allele alignment is already known. Check the effect of the old filter separately.
- Compare CHEWS and HEELS using LD calculated from the same simulated individuals. Then change only CHEWS's LD input to the reference panel to measure the cost of that substitution. The current `sim` LD option is the population correlation matrix used to generate genotypes; it is **not** the realized in-sample LD, so this requires a small pipeline change.
- Compare SuSiE effect limits of 10, 20, and 40, with the other inputs unchanged. Calculate both CHEWS formulas from each fitted model, using the fuller LD-aware calculation as the proposed primary result.

Report average bias, average absolute error, how often CHEWS is closer, and failed or extreme estimates. Change factors in controlled comparisons so we can identify which change matters.

**What this would tell us:** If the shortfall largely disappears, the paper should emphasize the relevant implementation or input limitation. If it persists after these checks, we have stronger evidence that distributed effects are intrinsically harder for this approach under the tested conditions. Either outcome makes the paper more informative.

### 3. Establish exactly what “more precise” and “more reliable” mean here

**Why this matters:** Lower variability is central to CHEWS's appeal, but lower variability is not automatically lower error or better uncertainty estimates.

**What I propose:** For the simulation check, keep the genotypes, causal variants, and effect sizes fixed, and generate repeated expression datasets by changing only the expression noise. Repeat this at several regions and genetic-effect settings. The existing repeat experiment redraws causal effects too, so it does not isolate this source of sampling variability.

Compare each method's average bias, variation across repeats, total estimation error, and how often its stated 95% interval contains the chosen true value. Define that true value consistently across repeats. This directly tests whether CHEWS's smaller reported errors reflect genuinely more stable estimates, overly narrow uncertainty intervals, or both.

For GTEx, repeat the donor partition, for example across twenty random splits. Report both correlation and the absolute differences between estimates from the two halves. Keep gene-selection rules consistent and show the all-gene analysis first. The splits reuse donors, so their results describe sensitivity to the partition; they are not twenty independent replications.

**What this would tell us:** We could distinguish three claims clearly: CHEWS varies less, CHEWS has lower total error, and CHEWS reports trustworthy uncertainty. The current evidence does not establish all three.

### 4. Strengthen the low-heritability analysis without trying to force a positive answer

**Why this matters:** This directly addresses your question about whether CHEWS's lower estimates are helpful. The current simulation suggests a modest benefit, but it failed the overall check against real data. Its preferred proportion of genes with no genetic contribution also differs substantially from the shuffled-data estimate.

**What I propose:** Extend the simulated true heritability range to 0.8 while retaining the low-heritability settings and the existing evaluation criteria. Treat this as a follow-up design change. Report whether it improves agreement across the full figure, including the high-heritability groups, and whether the disagreement about the no-effect proportion remains.

For a stronger validation, use one set of genomic regions to choose the mixture of simulated gene types and a separate set to evaluate it. Keep overlapping regions together so that nearby genes do not appear on both sides. This is separate from the donor split used to define gene groups. It would test whether the simulation reproduces patterns beyond the genes used to tune it.

In the low-heritability group, report **how much** error changes as well as how often CHEWS is closer. Also separate genes with truly zero simulated heritability from those with small but nonzero heritability. A method could benefit the former by reducing positive estimates while overlooking real effects in the latter. Show sensitivity to the assumed mixture rather than presenting one fitted mixture as known biology.

**What this would tell us:** If the benefit survives these checks, we can make a more convincing, explicitly simulation-based argument. If it does not, the useful conclusion is that we cannot yet distinguish correction of overestimation from loss of weak genetic signal in these real genes.

### Suggested order and decisions

I would make **observable-feature-based method selection in proposal 1 the main research priority**. We can begin with the existing results to see whether a simple, interpretable rule improves estimation on held-out regions. Proposal 2 then becomes an essential validation of that rule: we need to know that it captures useful patterns rather than peculiarities of the current LD inputs, variant filtering, or effect limit. We should revisit the rule if those changes alter which method performs better.

Proposal 3 supports the precision and reliability claims and can test how stable the recommendations are across donor splits. Proposal 4 addresses the unresolved low-heritability question; it becomes especially important if the proposed rule recommends CHEWS for genes with weak or undetected signals.

I would add another GTEx tissue after the rule and evaluation plan are fixed, to assess its applicability and stability beyond whole blood. Real-data replication still cannot establish accuracy without known truth. Allele-frequency validation should move earlier only if that finding becomes central to the guide. I would estimate computing time from a small pilot before committing to full reruns.

For our discussion:

1. Can we make a practical, interpretable guide for choosing a method the paper's main goal, with the benchmark explaining and validating that guide?
2. Should the first milestone be a simple CHEWS-versus-GCTA rule using the saved simulation results, followed by a CHEWS-versus-HEELS comparison for summary-statistic users?
3. What reduction in error would make a recommendation useful, and when should the guide recommend reporting both estimates instead of choosing one?
