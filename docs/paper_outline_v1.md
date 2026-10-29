# Paper outline for discussion

**Working title:** *CHEWS improves cis-heritability estimation when genetic effects are concentrated in a few variants*

**Authors:** Sunan Xu, Stephen Golzari, Tiffany Amariuta-Bartell

**Draft date:** October 28, 2026

**Purpose of this outline:** Agree on the main scientific message, the order of the results, and the analyses needed to support the final paper. “Completed” below means that an analysis has been run; it does not mean that every interpretation has been established.

## 1. The main message

We want to estimate how much of the variation in a gene's expression is explained by nearby genetic variants. This quantity is called **cis-heritability**, written as cis-h². For example, cis-h² = 0.10 means that nearby variants explain 10% of the expression variance under the model being used.

CHEWS uses the output of SuSiE, a fine-mapping method, to estimate cis-heritability from association summary statistics and information about correlations between variants. We compare it with GCTA, which uses individual-level data, and HEELS, another method that uses summary statistics.

**Our main finding is a trade-off: CHEWS can give less variable and more accurate estimates when the genetic contribution is concentrated in a few variants, but it tends to underestimate heritability when that contribution is spread across many variants.** Its advantage depends on both sample size and the distribution of genetic effects.

In GTEx whole blood, CHEWS broadly agrees with GCTA and gives estimates that are more consistent between two separate groups of donors. However, consistency and agreement do not establish accuracy: we do not know the true heritability in the real data. The additional simulation intended to address this question remains exploratory because it did not reproduce the overall real-data pattern closely enough.

The paper should explain where CHEWS is useful, where it falls short, and what additional evidence is needed to choose between methods in practice.

## 2. Terms we should keep distinct

| Term | Question it answers | How we assess it |
|---|---|---|
| Accuracy | Is the estimate close to the true value? | Simulations with known truth |
| Bias | Does a method systematically estimate too high or too low? | Average estimate minus the true value |
| Precision | How much does an estimate vary when we repeat an experiment under the same conditions? | The current repeat experiment measures variability within a region, including changes in the simulated causal effects |
| Agreement | Do two methods give similar values? | Differences and relationships between their estimates in GTEx |
| Reliability | Are estimates consistent in separate samples? | Correlation between estimates from two non-overlapping donor groups |

**Mean absolute error (MAE)** is the average distance from the truth, ignoring whether an estimate is too high or too low. **Win rate** means the percentage of paired simulated cases in which CHEWS is closer to the truth than GCTA. It does not mean the percentage in which CHEWS is best among all three methods.

**LD**, or linkage disequilibrium, describes correlations between genetic variants. **In-sample LD** is calculated from the same individuals who supplied the association statistics. **Reference-panel LD** comes from a separate dataset. This distinction matters for interpreting the current comparisons.

## 3. Draft abstract

Cis-heritability measures the contribution of nearby genetic variants to variation in gene expression. Estimating it from association summary statistics could extend its use when individual-level data are unavailable. We introduce CHEWS (Cis-Heritability Estimation With SuSiE), which uses fine-mapping estimates of genetic effects and their uncertainty to estimate cis-heritability. We compared CHEWS with GCTA and HEELS in 20,294 paired simulation fits and applied the methods to GTEx v8 whole blood from 574 donors of European ancestry. In baseline simulations with one causal variant, CHEWS had little average bias and was closer to the truth than GCTA in 76% of cases. With 20 causal variants, CHEWS underestimated heritability by 0.13 on average and was closer to the truth in only 27% of cases. In the sample-size experiments, CHEWS retained an advantage at one causal variant across sample sizes of 100–1,000, while GCTA generally performed better at 20 causal variants and sample sizes of at least 300. In GTEx, CHEWS and GCTA gave broadly similar estimates, and CHEWS estimates were more consistent between two non-overlapping donor groups (correlation 0.912 versus 0.865). These real-data results establish agreement and reliability, rather than accuracy. Our findings support CHEWS as a useful approach when genetic effects are concentrated, while highlighting its tendency to miss heritability when effects are distributed across many variants.

*This abstract reflects the current analyses. Update it after deciding which version of CHEWS to use throughout the paper and completing the key comparison described below.*

## 4. Introduction: why this question matters

1. **Start with the biological quantity.** Explain cis-heritability and its role in studying gene regulation, selecting genes for expression prediction, and informing transcriptome-wide association studies (TWAS). Keep the focus on what a useful estimate enables researchers to do.
2. **Explain the practical problem.** Analysts may have association summary statistics but lack access to individual-level genotypes and expression data. Introduce GCTA as the individual-level comparison and HEELS as a summary-statistic comparison.
3. **Introduce the opportunity.** Fine-mapping already estimates which variants may affect expression and how large their effects may be. CHEWS uses those estimates, including their uncertainty, to estimate the total genetic contribution.
4. **State the study question.** When does this approach improve estimation, and when does it miss part of the genetic contribution? Explain that simulations assess accuracy, while GTEx assesses agreement and reliability.

The earlier outline lists Gusev 2016, Gamazon 2015, GTEx 2020, Amariuta 2023, Yang 2011, Wang 2020, Zou 2022, and Li, Mazumder & Lin 2023 as starting references. Check the relevant papers before making a novelty claim such as “the first study to evaluate this approach.”

## 5. Results: the proposed scientific sequence

### Result 1. How CHEWS turns fine-mapping results into a heritability estimate

**Question:** What does CHEWS calculate, and does its calculation account for correlations between variants?

CHEWS combines estimated genetic effects with the uncertainty around those effects. The current main analyses use a simpler calculation, called the **trace form**. We also evaluated a fuller calculation that includes the contribution of correlations between variants to the final genetic-variance estimate, called the **LD-aware form**. Both use LD during fine-mapping; the difference is in the final calculation of heritability.

In 3,000 paired fits at a sample size of 300, the two versions were strongly correlated (r = 0.996). Their median absolute difference was 0.002, but the largest difference was 0.368. Including the additional terms did not remove the average underestimation at 20 causal variants.

**Takeaway:** The simpler calculation gives similar results in most tested cases, but the versions are not interchangeable for every gene. We should choose the primary version before finalizing the manuscript.

**Status:** Completed at sample size 300; not yet evaluated throughout all analyses.

**Figure 1:** CHEWS workflow and comparison of the two calculations. Put the derivation in Methods.

### Result 2. CHEWS performs best when a few variants account for the genetic contribution

**Question:** How does performance change as the same total genetic contribution is spread across more variants?

The baseline experiment used 300 individuals, five causal-variant counts, six true heritability values, and 100 regions per setting. Each method analyzed the same simulated expression data.

| Baseline setting | CHEWS average bias | Cases where CHEWS is closer than GCTA |
|---|---:|---:|
| One causal variant | −0.008 | 76% |
| Twenty causal variants | −0.130 | 27% |

With 20 causal variants, GCTA and HEELS had much smaller average biases: −0.005 and +0.003, respectively. The fuller LD-aware CHEWS calculation still had a bias of about −0.130.

**Takeaway:** CHEWS works well in the tested settings with concentrated effects and misses part of the genetic contribution when effects are more distributed. This is consistent with the assumptions of the fine-mapping model, but the current experiment does not isolate the cause: CHEWS received reference-panel LD, whereas HEELS received in-sample LD. The current SuSiE setting also allows ten effect components, while the simulations include up to twenty causal variants.

**Status:** Comparison completed; the explanation needs further testing.

**Figure 2:** Average bias and the percentage of cases favoring CHEWS across causal-variant counts and true heritability values.

### Result 3. Sample size changes the trade-off between bias and variability

**Question:** Can a method with some downward bias still give a better estimate because its estimates vary less?

In a separate experiment, we generated 25 expression datasets per setting at each of 30 regions. With 100 individuals and one causal variant, the average within-region standard deviation was **0.108 for CHEWS and 0.246 for GCTA**. The code holds the region and simulation settings fixed but redraws the causal variants and effects for each replicate, as well as expression noise. These numbers therefore measure variability across repeated simulations within a region, rather than sampling variability for a single fixed set of genetic effects.

Across the sample-size experiment:

- With **one causal variant**, CHEWS was closer to the truth in **75–82% of cases** at every tested sample size from 100 to 1,000. Its average error was **73% lower** than GCTA's at 100 individuals and **58% lower** at 1,000.
- With **twenty causal variants**, GCTA was closer in **72–77% of cases** at sample sizes of 300 or more.
- When averaging over **all causal counts and true heritability values**, CHEWS's average absolute error was 0.142 versus GCTA's 0.205 at 100 individuals, and 0.039 versus 0.044 at 1,000. These averages describe a broader set of cases than the one-causal-variant percentages above.

HEELS also showed a limitation at the smallest sample size: at 100 individuals, 46% of its estimates exceeded 0.9 despite the algorithm reporting convergence. This needs to be reported alongside its good performance in other settings.

**Takeaway:** Sample size alone does not determine which method performs better. CHEWS's lower variability is most useful when its downward bias is limited, particularly when effects are concentrated.

**Status:** Within-region repeat experiment completed. A claim about sampling precision at fixed genetic effects would require an additional experiment holding those effects fixed.

**Figure 3:** Variability under repeated simulation and the percentage of cases favoring CHEWS by sample size and causal-variant count.

### Result 4. The distribution of effect sizes matters, beyond the number of causal variants

**Question:** Do all settings with the same number of causal variants behave alike?

No. Holding total genetic variance fixed, CHEWS performed better when a few variants carried larger effects than when the effects were equally distributed. Its average bias was **−0.030 with concentrated, heavy-tailed effects**, compared with **−0.094 with equal-magnitude effects**. It was closer than GCTA in 61% versus 43% of cases, respectively.

Two additional checks give context:

- Restricting causal variants to common variants gave GCTA an average upward bias of **+0.045**. This shows why GCTA should not be treated as known truth, but this allele-frequency result needs validation with realistic discrete genotypes.
- Selecting regions with high versus low LD changed CHEWS's average bias by only **0.009** in this experiment. That does not establish that LD mismatch is unimportant; regional correlation strength and a mismatch between LD datasets are different questions.

**Status:** Completed; the allele-frequency result needs validation before a strong biological interpretation.

**Figure 4:** Effect-size concentration, with allele frequency and regional LD as additional comparisons or supplementary panels.

### Result 5. In GTEx whole blood, the methods broadly agree, and CHEWS is more consistent across donor groups

**Question:** How do the methods behave on real gene-expression data?

We analyzed **19,432 genes from 574 donors of European ancestry**, using five nested sample sizes from 100 to 574. GCTA returned estimates for 96.5–98.6% of genes across these subsets, so the number entering each comparison is smaller than the number analyzed. Both summary-statistic methods used in-sample LD.

CHEWS estimates were lower than GCTA's by **0.015 on average at 100 individuals** and **0.007 at 574 individuals**. HEELS agreed very closely with GCTA at 574 individuals (r = 0.9996), although their values were not identical.

We then divided the donors into two non-overlapping groups of 287 and estimated heritability separately in each group.

| Genes included | GCTA correlation between halves | CHEWS | HEELS |
|---|---:|---:|---:|
| All genes with usable estimates for the relevant method | 0.865 | **0.912** | 0.870 |
| Genes with at least one detected cis-eQTL in the full dataset | 0.850 | **0.898** | 0.857 |
| Genes with estimates above 0.05 for all three methods in both halves | 0.848 | **0.892** | 0.855 |

The last group was selected using the estimates themselves and should be presented as a sensitivity analysis. The cis-eQTL group was defined using the full dataset, not an independent cohort.

**Takeaway:** CHEWS gives more consistent estimates across these two donor groups. This supports reliability, but an estimator can be consistently biased. The result should therefore be interpreted alongside the simulation accuracy results. We have tested only one donor partition so far.

**Status:** Completed for one tissue and one donor partition.

**Figure 5 and Table 1:** Agreement across sample sizes and reliability across donor halves.

### Result 6. We are still testing whether lower CHEWS estimates are better for genes with low estimated heritability

**Question:** When CHEWS gives a lower value than GCTA, is it correcting an overestimate or missing real genetic effects?

First, the way genes are grouped affects the apparent answer. If we group genes by a noisy GCTA estimate and compare methods using the same donors, the grouping itself can create a pattern of differences. We therefore used one donor half to define the groups and the other to compare the methods. Switching to this approach changed the group-average differences by 0.023 on average.

We then ran two follow-up analyses:

1. **Remove the genotype–expression relationship.** After shuffling genotype donor labels, constrained GCTA and CHEWS both gave average estimates of about **0.013**. A positive estimate can therefore occur even without a genetic signal. This result does not support attributing the real-data gap solely to GCTA's restriction against negative estimates.
2. **Simulate expression using real GTEx genotypes.** With a known genetic contribution, CHEWS was closer to the truth in about **60% of simulated cases** (95% confidence interval: 59–61%) corresponding to the real group with no detected signal and estimated heritability below 0.05. That real group contains 9,572 genes. Across the tested assumptions about the overall proportion of genes with no genetic contribution, the percentage ranged from **54% to 61%**.

**The important limitation:** The second analysis did not reproduce the overall real-data pattern closely enough. Its average mismatch in the CHEWS–GCTA differences was **0.014**, above the **0.010** limit set before examining the simulation results. The largest discrepancies were in groups with high estimated heritability. The quoted confidence interval describes uncertainty within the simulation; it does not capture uncertainty about whether that simulation represents the real data.

There is also uncertainty about how many genes truly have no genetic contribution. The shuffled-data comparison gives a rough estimate of 73% within the low-heritability, no-detected-signal group and 50% across the analyzed genes overall. These estimates depend on the cutoff used and can include genes with weak, undetected effects. They are not observed ground truth. The simulation's best-fitting overall proportion was only 10%, adding to the reasons for caution; overall and group-specific proportions should not be used interchangeably.

**Takeaway:** The simulation suggests that lower CHEWS estimates may be modestly helpful in this group, but it does not yet establish greater accuracy in real GTEx data. Extending the simulated heritability range from 0.4 to 0.8 is a reasonable next check, not a guaranteed solution. Even if that check passes, the accuracy conclusion remains dependent on the simulation assumptions.

**Status:** Shuffled-data analysis completed; simulated accuracy result exploratory.

**Figure 6:** Gene groups defined in separate donors and the shuffled-data comparison. Include simulated accuracy in the main text only if follow-up validation supports it; otherwise present it as exploratory in the Supplement.

## 6. Discussion: what the findings mean

The Discussion should begin with the conditional benefit of CHEWS: it can reduce estimation error when genetic effects are concentrated, while underestimating heritability when effects are more distributed. The effect-size experiment is particularly useful here because it shows that causal-variant count alone does not explain performance.

Next, discuss the practical opportunity to estimate heritability from fine-mapping results already available in eQTL studies. Compare CHEWS with GCTA, HEELS, and other local-heritability methods after checking the relevant literature. Avoid recommending a method solely because summary statistics are available; the source and quality of LD also matter.

Then explain what the GTEx analyses add: evidence of agreement and reliability in one tissue, plus an unresolved question about accuracy for genes with low estimated heritability. More reproducible estimates may be useful, but reproducibility alone does not show that a method captures the full genetic contribution.

The main limitations to address are:

- **The current simulations do not isolate why methods differ.** LD source, the number of allowed SuSiE effects, and variant filtering may contribute alongside the underlying modeling assumptions.
- **Most simulations use continuous, normally distributed genotype values.** They reproduce the chosen correlation structure but not discrete genotypes or haplotypes. This especially limits the interpretation of the allele-frequency result.
- **The real-data evidence comes from one tissue and one donor partition.** The GTEx analysis also uses a smaller set of variants per region than the main simulations.
- **Reported uncertainty is too small in some settings.** In the reported simulation summary, nominal 95% intervals included the true value about 87% of the time for GCTA and 68% for CHEWS. Lower standard errors should not be presented as automatically better uncertainty estimates.
- **Numerical failures and extreme estimates matter.** CHEWS occasionally returns estimates above the possible heritability range despite reporting convergence. Report these cases and how they affect average errors, rather than silently excluding them.
- **GCTA is fitted differently in simulation and GTEx.** The main simulations allow unconstrained estimates; GTEx uses constrained estimates because unconstrained fitting frequently failed with covariates. This matters especially near zero.

End with a practical conclusion tied to the tested conditions, rather than a universal recommendation that one method is best.

## 7. Methods: details needed for reproducibility

Keep the main Results readable and place the following details in Methods:

- **CHEWS calculation and model settings.** Define the trace form, Σⱼ(E[βⱼ]² + Var[βⱼ]), and the LD-aware form, E[βᵀRβ | data], under the fitted SuSiE posterior. Give the derivation, scaling assumptions, and validation. Report the current SuSiE settings, including ten effect components and a maximum of 100 iterations, and explain any follow-up changes.
- **Comparison methods.** Describe GCTA's covariates and constraints, and the project's HEELS implementation and its requirement for in-sample LD. Document checks against the published algorithm.
- **Main simulations.** Describe the 1000 Genomes European reference data, the ±500 kb regions, genotype simulation, scaling to known genetic variance, and all experimental settings. Distinguish the 20,400 attempted fits from the 20,294 fits with estimates from both CHEWS and GCTA. Explain reuse of regions across settings and bootstrap confidence intervals that resample whole regions.
- **Repeated simulations within regions.** Explain what is held fixed and what is redrawn. In the current code, causal variants and effects are redrawn for each replicate. If we add a test of sampling precision alone, hold the genetic effects fixed and generate repeated expression noise.
- **GTEx analysis.** Describe donor and gene inclusion, variant selection, covariate adjustment, nested subsets, non-overlapping donor halves, and the association-based signal count. A detected signal is not the same as a known causal variant.
- **Agreement and reliability.** Define correlations, average differences, regression comparisons, and resampling procedures. State the gene set and usable sample size for each result.
- **Shuffled-data and simulated-expression analyses.** Describe how genetic signal is removed or generated, the assumed mixture of genes, the target heritability used for accuracy, and the sensitivity check using heritability after covariate adjustment. Report the evaluation criteria set before the first analysis and any subsequent changes to the design.
- **Code and data availability.** Identify the analysis code and explain the controlled-access requirements for GTEx data.

**Version note:** The current GTEx analysis passes covariates directly to GCTA as fixed effects. An earlier analysis adjusted expression before fitting GCTA and gave an inflated small-sample difference. The current difference at 100 individuals is −0.015, replacing −0.151 in the submitted IGVF abstract. Use the corrected results consistently and document the change where relevant.

## 8. Decisions for discussion with Tiffany

We can begin drafting the established findings now. The following decisions determine how strong the final claims can be.

| Decision | Why it matters | Suggested next step |
|---|---|---|
| Compare CHEWS and HEELS using the same in-sample LD | Separates differences in their models from differences in LD inputs | Prioritize this comparison. Test CHEWS with reference-panel LD separately; do not treat HEELS with reference-panel LD as a valid comparison under the current implementation. |
| Choose the primary CHEWS calculation | The fuller formula includes LD contributions that can matter for individual genes | Prefer the LD-aware calculation for the final paper, after evaluating it throughout the main analyses; retain the trace form as a comparison. |
| Check the ten-effect limit and variant filtering | These choices may contribute to CHEWS's performance at high causal counts | Test a larger effect limit and removal of the unnecessary strand-ambiguous-variant filter before attributing the shortfall to the model alone. |
| Extend the GTEx simulation to higher true heritability | May explain the mismatch in high-heritability groups | Extend to 0.8, retain the evaluation criteria, and report this as a follow-up design change. Keep the result exploratory if the mismatch persists. |
| Repeat the donor split | Tests whether the reliability result depends on one partition | Repeat across multiple random partitions, for example twenty. |
| Clarify the precision claim | The current within-region repeats also redraw causal effects | Describe the current result as within-region variability, or add repeats with fixed causal effects to isolate sampling precision. |
| Validate the allele-frequency result with discrete genotypes | The current genotype simulation limits its interpretation | Do this before making it a central claim; otherwise keep it exploratory or supplementary. |
| Add another GTEx tissue | Tests whether the findings extend beyond whole blood | Decide whether replication in another tissue is essential for the paper's intended scope. |
| Agree on the title and target journal | Determines the scope and level of evidence to prioritize | Choose after agreeing on the main claim and follow-up analyses. |

## 9. Proposed figures and tables

| Item | Main purpose | Current status |
|---|---|---|
| Figure 1 | Explain CHEWS and compare its two calculations | Comparison completed at sample size 300 |
| Figure 2 | Show how accuracy changes with causal-variant count and true heritability | Completed; interpretation needs the LD comparison |
| Figure 3 | Show how sample size changes variability and accuracy | Completed |
| Figure 4 | Show why effect-size concentration matters | Completed; allele-frequency validation pending |
| Figure 5 | Show agreement among methods in GTEx | Completed for whole blood |
| Figure 6 | Explain the low-heritability question and follow-up checks | Simulated accuracy remains exploratory |
| Table 1 | Compare consistency across donor halves | Completed for one partition |
| Supplement | Detailed error distributions, interval coverage, extreme estimates, gene inclusion, LD-aware comparison, and GTEx simulation diagnostics | Existing analyses provide much of the material |

**Evidence behind this outline:** The numerical results come from the saved simulation and GTEx results, the tables in `poster/tables/`, and the round-4 analysis record in `docs/pi_feedback_round4_implementation.md`. The project's `poster/check_numbers.py --round 4` was run during this revision and passed all 91 numerical checks. Those checks confirm the numbers they cover; they do not establish the broader interpretations or validate every proposed manuscript claim.
