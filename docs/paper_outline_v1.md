# Manuscript outline (v1, for review)

**Working title:** *CHEWS reduces cis-heritability error in small samples when regulatory effects are concentrated*

**Authors:** Sunan Xu, Stephen Golzari, Tiffany Amariuta-Bartell · **Draft:** 2026-10-28

**What I would like from this review:** (1) whether the story and section order are right; (2) whether each claim is stated at the strength the evidence supports; (3) the decisions in §8, especially which missing analyses must be done before writing.

**How to read the status tags:**
- **[done]**: analysis finished; every number is re-derived by the poster's `check_numbers.py`.
- **[partial]**: run, but only in part of the design.
- **[needed]**: not yet run.

---

## 1. The one-paragraph story

Cis-SNP heritability (cis-h²) of gene expression is used to prioritize genes and to decide which genes to model in TWAS and related analyses. The standard estimator, GCTA/GREML, needs individual-level genotypes, but consortium eQTL data are increasingly shared only as summary statistics. We introduce CHEWS (Cis-Heritability Estimation With SuSiE), which estimates cis-h² from the posterior second moments of fine-mapped effects. We benchmark it against GCTA (individual-level, infinitesimal model) and HEELS (summary-level, infinitesimal model).

Across 20,294 paired simulation fits, CHEWS trades variance for downward bias. With concentrated cis effects and small samples it is closer to the truth than GCTA on most loci, whereas it increasingly underestimates h² as effects become diffuse. In GTEx whole blood, where no truth is available, CHEWS and GCTA broadly agree. CHEWS estimates were more reproducible across disjoint donor halves, and the residual disagreement is concentrated in genes with low estimated h².

The practical conclusion is conditional: choose the estimator according to data availability, LD source, sample size and regulatory architecture.

---

## 2. Terminology used throughout

Your earlier reviews asked for these to be kept distinct. I propose fixing them in the first Results paragraph and never using them interchangeably.

| term | meaning | how we measure it | where |
|---|---|---|---|
| **Accuracy** | closeness to truth | bias, MAE; pairwise win rate | simulation only |
| **Precision** | sampling variability at a fixed locus | SD over repeated phenotypes at the same region | simulation |
| **Agreement** | similarity of values between methods | mean offset, Deming slope, RMSE | GTEx |
| **Reliability** | reproducibility across samples | split-half correlation on disjoint donors | GTEx |

Other conventions:
- **Win rate** always means *pairwise*: P(\|CHEWS − truth\| < \|GCTA − truth\|). The pipeline also writes a three-way "closest of three" share (70% at one causal SNP vs 76% pairwise), which the paper does not use.
- "CHEWS" means the trace form unless stated; "exact CHEWS" is the LD-aware form (§4, Results 1).
- MAE is reported next to RMSE, because single runaway fits can dominate squared error.

---

## 3. Abstract (skeleton in the target register)

> Cis-SNP heritability (cis-h²) quantifies the genetic regulation of gene expression and informs gene prioritization. However, standard estimators such as GCTA require individual-level genotypes, which are often unavailable for consortium-scale eQTL data. In the present study, we introduce CHEWS, which estimates cis-h² from the posterior second moments of effects fine-mapped with SuSiE using summary statistics and LD. We benchmarked CHEWS against GCTA and HEELS in 20,294 paired simulations across 350 1000 Genomes EUR cis-regions, and applied all three estimators to 19,432 genes in GTEx v8 whole blood (574 individuals of European ancestry). With a single causal variant, CHEWS was nearly unbiased (bias = −0.008) and closer to the truth than GCTA at 76% of loci; with 20 causal variants, CHEWS underestimated cis-h² (bias = −0.13), whereas GCTA and HEELS remained nearly unbiased. CHEWS estimates were 2–4× less variable than GCTA estimates at a fixed locus, and at one causal variant CHEWS was closer to the truth at 75–82% of loci at every sample size from 100 to 1,000. In GTEx, CHEWS and GCTA differed by −0.015 at n = 100 and −0.007 at n = 574, and CHEWS estimates were more reproducible across disjoint donor halves (split-half r = 0.91 versus 0.86). In conclusion, CHEWS is a precise estimator of cis-h² for concentrated regulatory effects when only summary statistics are available, whereas infinitesimal-model estimators better capture regulation distributed across many weak effects.

Note: the GTEx numbers supersede the submitted IGVF abstract (−0.151 → −0.015 at n = 100). The analysis changed after submission: covariates now enter REML as fixed effects instead of pre-residualizing expression. That change should be explained once, in Methods.

---

## 4. Introduction (four paragraphs)

1. **Why cis-h² of expression matters.** Its uses: gene prioritization, selecting cis-heritable genes for TWAS and expression prediction (refs: Gusev 2016, Gamazon 2015, GTEx 2020), tissue co-regulation analyses (Amariuta 2023 used GCTA cis-h² to define significantly cis-heritable genes). The quantity is local and per gene, typically estimated with GCTA/GREML (Yang 2011).
2. **The data problem.** Individual-level genotypes are increasingly unavailable to the analyst (consent, dbGaP, consortium sharing); summary statistics plus LD are what is shared. Summary-level local-heritability estimators exist. HEELS (Li, Mazumder & Lin 2023) recovers the REML target but requires in-sample LD. Fine-mapping (SuSiE, Wang 2020; susie_rss, Zou 2022) is already run routinely on these data and yields posterior effect moments.
3. **The gap.** To our knowledge, no study has evaluated fine-mapping posteriors as a cis-h² estimator against individual-level REML. The same applies to *when* a sparse-effect estimator should be preferred over an infinitesimal one: across sample size, effect concentration, allele frequency and LD.
4. **This study.** We introduce CHEWS, derive its exact LD-aware estimand, and benchmark it against GCTA and HEELS. Including HEELS allows a partial separation of data representation (summary vs individual-level) from the effect model (sparse vs infinitesimal). We then test the simulated predictions in GTEx whole blood, measuring agreement and reliability rather than accuracy.

---

## 5. Results

### R1. Overview of CHEWS [done]
- What CHEWS computes. SuSiE writes β = Σₗ βₗ, with each single effect selecting one variant. With mₗ = E[βₗ] = αₗ⊙μₗ and sₗ = E[βₗᵀRβₗ]:
  **E[βᵀRβ | data] = μᵀRμ − Σₗ mₗᵀRmₗ + Σₗ sₗ**.
  The pipeline's trace form Σⱼ(μⱼ² + σⱼ²) is what this reduces to at R = I.
- The exact form was checked against brute-force tr(R·E[ββᵀ]) to 2×10⁻¹⁶.
- On 3,000 paired fits at n = 300 the two forms track each other on average: r = 0.996; median \|difference\| 0.002, 95th percentile 0.034, maximum 0.368. Their mean bias curves overlap at every causal count.
- **Fig. 1**: workflow (summary statistics + LD → SuSiE posterior → posterior genetic variance → cis-h²) and the trace-vs-exact comparison.
- [partial] The exact form has been evaluated at n = 300 only (see §8, decision 2).

### R2. CHEWS is nearly unbiased for concentrated effects and underestimates h² as effects become diffuse [done]
- Baseline grid: n = 300, 5 causal counts × 6 true h² values, 50 regions per cell.
- One causal SNP: CHEWS bias −0.008 and closer than GCTA at 76% of loci.
- Twenty causal SNPs: bias −0.130 for both trace and exact CHEWS, and CHEWS closer at only 27% of loci. GCTA (−0.005) and HEELS (+0.003) stay near zero.
- Win landscape over h² × causal count: from 85% (1 causal, h² = 0.4) down to 3% (20 causal, h² = 0.8).
- **Interpretation, as calibrated on the poster:** this is *consistent with* the sparse-effect model. Including the LD cross terms does not remove the shortfall. Attribution to sparsity remains **provisional**, because CHEWS (reference-panel LD) and HEELS (in-sample LD) differ in LD source as well as model. The matched-LD factorial ([needed], §7) resolves this.
- **Fig. 2**: (a) bias vs causal count, three estimators plus exact CHEWS, region-clustered bootstrap CIs; (b) pairwise win-rate heatmap (h² × causal count).

### R3. Lower sampling variability offsets the downward bias mainly in small samples with concentrated effects [done]
- **Precision, in a clean design:** 25 phenotypes per region, 30 regions. Within-locus SD is 0.108 for CHEWS vs 0.246 for GCTA (n = 100, one causal SNP), and CHEWS is 2–4× lower across conditions. Within- and across-locus SDs agree to ≤ 0.007.
- **Bias falls with n:** CHEWS −0.122 (n = 100) → −0.021 (n = 1,000). MAE 0.142 vs 0.205 for GCTA at n = 100, and 0.039 vs 0.044 at n = 1,000.
- **How often vs by how much (round 4):** at one causal SNP, CHEWS is closer at 75–82% of loci at every n; what shrinks with n is the margin (MAE −73% at n = 100, −58% at n = 1,000). At 20 causal SNPs GCTA is closer at 72–77% of loci for n ≥ 300. So the advantage is an **interaction** of small n *and* concentrated effects, not a rule about n alone.
- HEELS at n = 100 (p/n ≈ 70): it converges on 100% of fits, but 46% of estimates exceed ĥ² = 0.9. That is saturation, not non-convergence, and should be reported as a result.
- **Fig. 3**: (a) dispersion / within-locus SD vs n; (b) win rate by n × causal count (the round-4 panel).

### R4. Effect concentration and causal allele frequency modify performance; regional LD has little effect [done]
- **Effect concentration** (total βᵀRβ fixed; only its distribution varies): CHEWS bias −0.030 for heavy-tailed effects vs −0.094 for equal-magnitude effects. Win rate 61% vs 43%.
- **Causal MAF:** GCTA becomes upward-biased (+0.045 for common causal variants), the only sign reversal in the study. This is consistent with the known GREML MAF/LD-stratification behaviour, and is a reason not to treat GCTA as an unconditional reference.
- **Regional LD:** top vs bottom LD decile changes CHEWS bias by only 0.009.
- **Fig. 4**; and **Supplementary**: failure profiles. Unconstrained GCTA leaves [0,1] on 22% of fits at n = 100; CHEWS leaves it on 0.1–0.3% of fits but can reach 3–5 while reporting convergence. HEELS is bounded by construction.

### R5. In GTEx whole blood, HEELS closely agrees with REML and CHEWS–GCTA disagreement is small [done]
- Data: 19,432 genes; 574 EUR donors; five nested subsets (100–574); in-sample LD; covariates as REML fixed effects. GCTA returned estimates for 96.5–98.6% of genes per subset.
- **Agreement with GCTA:**
  - CHEWS mean offset: −0.015 (n = 100) → −0.007 (n = 574).
  - CHEWS Deming slope: 0.67 → 1.06.
  - HEELS at n = 574: r = 0.9996, slope 0.93, offset −0.0074, RMSE 0.012. That is close, but **not exact**.
- **Fig. 5**: offset and calibration slope vs n.
- **Wording, as in your round-1 review:** agreement in whole blood is *consistent with* predominantly concentrated cis effects, but does not independently establish sparsity.

### R6. Disagreement in low-h² genes: a conditioning artefact, a shared null floor, and (pending) simulated truth [partial]
- **Conditioning:** stratifying on GCTA from the *same* donors manufactures much of the apparent pattern. Strata defined on a disjoint half change the cell means by 0.023 on average. No Simpson's paradox is claimed.
- **Permutation null, round 4 [done]:** with genotypes permuted (true h² = 0), constrained GCTA and CHEWS both average 0.013. GCTA sits at its bound for 57% of genes; CHEWS is below 0.005 for 40%. Estimated share of truly null genes (π̂₀): 0.73 in the 0-signal, h² < 0.05 stratum; 0.50 genome-wide.
- **Semi-synthetic GTEx, round 4 [partial]:** phenotypes simulated on the real genotypes give a known truth per stratum. In the 0-signal, h² < 0.05 stratum, CHEWS was closer for 60% of genes (59–61%), and for 54–61% across any null fraction. **However, it failed the pre-specified calibration check** (mean miss 0.014 vs a 0.010 limit). The miss is in the high-h² corner, probably because the simulated h² stopped at 0.4.
  - If a re-run with h² up to 0.8 passes, this becomes the "better, not just different" figure.
  - If it does not pass, it goes to the Supplement as exploratory.
- **Fig. 6**: (a) disjoint-donor strata; (b) permutation null distributions; (c) simulated win rate by stratum, conditional on calibration.

### R7. CHEWS estimates are more reproducible across disjoint donor halves [done; one split]
- Two disjoint halves of 287 donors. Split-half Pearson r for GCTA / CHEWS / HEELS:
  - all genes: 0.865 / **0.912** / 0.870
  - prespecified genes with ≥ 1 cis-eQTL: 0.850 / **0.898** / 0.857
  - post-hoc ĥ² > 0.05 in all six estimates: 0.848 / **0.892** / 0.855
- Paired gene bootstrap: Δr = +0.050 (95% CI +0.038 to +0.063).
- **Interpretation:** reliability is not accuracy. A more shrunken estimator is expected to be more reproducible, so R7 must be read together with R2–R3. One donor split only; [needed] repeated partitions (§7).
- **Table 1** (or Fig. 7).

### R8. A conditional guide to choosing an estimator [done]
**Table 2**, conditional on both data and architecture:

| situation | interpretation (current evidence) |
|---|---|
| summary statistics + exact in-sample LD | HEELS targets the REML equations; verify agreement |
| reference-panel LD only | assess LD match; not isolated in this benchmark |
| small n, concentrated signal | CHEWS may lower MAE despite downward bias |
| diffuse or weakly localizable signal | CHEWS underestimates h²; prefer an infinitesimal model or report both |
| strong MAF dependence | stratified GRM, after validation with discrete genotypes |

---

## 6. Discussion (planned paragraphs)

1. **Summary.** "We developed CHEWS…; we benchmarked…; we found…", in the register of the TCSC Discussion.
2. **Relation to existing estimators.** GCTA/GREML; GREML-LDMS for MAF stratification; HEELS and its in-sample LD requirement; summary-level local h² methods (e.g. HESS); why fine-mapping posteriors are attractive (already computed, interpretable per variant).
3. **Why sparsity helps and hurts:** variance reduction by localization vs missing unlocalized variance. Link to the effect-concentration result: causal count alone does not predict error.
4. **Implications for eQTL resources and TWAS:** gene selection by cis-h²; consortium settings with summary statistics only.
5. **Limitations**, enumerated ("We note several limitations of our work. First, …"):
   1. Attribution to the sparse-effect model is provisional until the matched-LD factorial is run.
   2. Simulated genotypes are multivariate normal: pairwise LD is preserved, but discrete genotypes and haplotypes are not. This matters most for the MAF result.
   3. One tissue, HapMap3 density (~400 SNPs per window vs ~8,200 simulated), one donor split.
   4. Nominal 95% Wald intervals under-cover (87% GCTA, 68% CHEWS). The CHEWS standard error reflects only the effects it localizes.
   5. Out-of-range CHEWS estimates are rare but extreme. We report them as failures rather than silently removing them.
   6. GTEx GCTA uses constrained REML (unconstrained failed for 60–75% of genes with covariates as fixed effects), while the simulations use unconstrained REML.
6. **Closing:** "Despite these limitations, …", ending with the conditional recommendation.

---

## 7. Methods (sections)

- **CHEWS estimator:** SuSiE/susie_rss settings (L = 10, coverage 0.95, max 100 iterations, `check_prior = TRUE`), the trace and exact estimands with derivation, standardization, and handling of out-of-range estimates.
- **GCTA/GREML:** binary GRM, unconstrained in simulation, constrained in GTEx with covariates as `--qcovar`.
- **HEELS:** our implementation (direct profile-likelihood maximization; same fixed point as the published iteration); in-sample LD requirement and why (the S′R⁺S ≤ y′y bound).
- **Simulation framework:** 350 regions of ±500 kb from 1000 Genomes Phase 3 EUR (503 individuals), MVN genotypes, effect rescaling so that βᵀRβ = h², grid axes, paired design, the precision design with 25 replicates, and region-clustered bootstrap (2,000 resamples).
- **GTEx analysis:** v8 whole blood, EUR, HapMap3; covariates (5 genotype PCs, inferred factors scaled with n per the GTEx rule, sex, platform, PCR); nested and disjoint subsets; LD clumping (p < 10⁻⁴, r² < 0.1) as a SuSiE-independent signal count.
- **Agreement and reliability statistics:** Deming regression, Bland–Altman, ICC(A,1), paired bootstrap.
- **Permutation null and semi-synthetic GTEx (round 4):** design prior, π̂₀ estimator, importance reweighting, and decision rules fixed before analysis.
- **Code and data availability:** `cish2` repository; GTEx is dbGaP controlled access.

---

## 8. Decisions needed before writing

1. **Matched-LD factorial** (CHEWS with in-sample LD; HEELS with panel LD). This is the single most important missing run: it decides whether R2 can say "sparse-effect model" or must stay "consistent with". About 25 min of compute. *Recommend: run before drafting Results.*
2. **Which CHEWS is the primary estimator?** Exact (LD-aware) everywhere, with trace as sensitivity (needs about 4 h to propagate through the sample-size and GTEx passes), or trace with the exact form as a check at n = 300. *Recommend: exact everywhere, since it is the correct estimand.*
3. **Semi-synthetic GTEx:** re-run with h² up to 0.8 under the same pre-registered rule, or move it to the Supplement as exploratory?
4. **Repeated donor splits** for R7 (e.g. 20 random partitions). About 2–3 h of compute. *Recommend: yes; a single split is the main reviewer target in R7.*
5. **Discrete-genotype validation** of the MAF result: resample real 1000G haplotypes for a subset of regions. Needed if R4's MAF finding stays in the main text.
6. **Scope of real data:** whole blood only, or add 2–3 GTEx tissues with larger and smaller n? A second tissue would turn "one tissue" from a limitation into a replication.
7. **`legacy_ambiguous_filter`:** current simulations drop strand-ambiguous SNPs before fine-mapping (about 15% of SNPs), which handicaps CHEWS unnecessarily. Re-run without it, or report as is?
8. **Title and target journal.** The working title follows your round-1 suggestion. Journal options: *AJHG*, *Genome Biology*, *Nature Communications*, *PLoS Genetics*.

---

## 9. Figure and table plan

| item | content | status |
|---|---|---|
| Fig. 1 | CHEWS workflow; exact vs trace estimand | done (exact at n = 300 only) |
| Fig. 2 | bias by causal count; win landscape | done |
| Fig. 3 | precision (within-locus SD); win rate by n × causal count | done |
| Fig. 4 | effect concentration, causal MAF, regional LD | done (MAF: validation needed) |
| Fig. 5 | GTEx agreement: offset and slope vs n | done |
| Fig. 6 | low-h² strata: disjoint donors, permutation null, simulated truth | partial |
| Table 1 | split-half reliability by gene set | done (one split) |
| Table 2 | conditional estimator guide | done |
| Supp. | failure profiles and out-of-range estimates; interval coverage by condition; trace-vs-exact distribution; HEELS saturation at n = 100; gene inclusion flow; agreement tables; calibration diagnostics for the semi-synthetic run | mostly done |

All main-text numbers above are pinned by `poster/check_numbers.py` (91 checks). The manuscript will reuse that mechanism so text and tables cannot drift.
