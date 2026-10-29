# Round-4 update for Tiffany: paste-ready

## Email version

**Subject:** Poster round 4: how often CHEWS wins, and what happens in the 0–0.05 stratum

Hi Tiffany,

Thank you for the comments on the poster. Here is where each one landed.

**1. Standard-error comparison (§5A) and split-half table (§8).** Both are unchanged, as you suggested.

**2. "Show how frequently CHEWS improves" (§5B).** The MAE-difference panel is now a win-rate heatmap: the % of loci where CHEWS is closer to the truth than GCTA, by sample size and causal-SNP count, on the same colour scale as the §4 heatmap.
- With one causal SNP, CHEWS is closer on 75–82% of loci at every sample size. What shrinks as n grows is the size of the advantage (MAE 73% lower at n = 100, 58% lower at n = 1,000), not how often it happens.
- With 20 causal SNPs the ordering reverses: GCTA is closer on 72–77% of loci once n ≥ 300. That row stays on the poster.
- On "% difference": true h² is not constant in that panel; it pools h² from 0.05 to 0.8. Dividing each error by h² lets the h² = 0.05 fits dominate, and in one cell it even flips the sign. The win rate answers "how often" with no rescaling, and % change in MAE answers "by how much".

**3. "Is the underestimation in the 0–0.05 stratum good or bad?" (§7).** Real data has no ground truth, so I ran two new GTEx analyses on the same donors, halves and covariates.
- **Permutation null.** I shuffled genotypes so that every gene's true h² is 0. Constrained GCTA and CHEWS both average 0.013. So GCTA's zero floor is not, on its own, why CHEWS reads lower there. This is now one sentence in §7.
- **Simulation on the real GTEx genotypes, with known truth.** In the 0-signal, h² < 0.05 cell (9,572 genes), CHEWS was closer to the truth on about 60% of genes (95% CI 59–61%). That held whatever share of genes I assumed truly null: 54–61%. The permutation run suggests about 70% of these genes are truly null, and at that level the figure is about 55%. So the provisional answer is "mildly good".
- **The caveat.** Before running it, I fixed a check for whether the simulation reproduces the real §7 pattern. It did not pass: the average miss was 0.014 against my limit of 0.010. The miss is almost entirely in the high-h² corner, most likely because the simulated h² only went up to 0.4. So the simulation is not on the poster yet. Re-running with h² up to 0.8 takes about 20 minutes of compute.

**Questions for you**
1. By "bottom-left," did you mean the 0-signal, h² < 0.05 cell? I analysed the whole column either way.
2. Should I re-run the simulation with the wider h² range before we print, or keep it as follow-up?
3. I think the results are ready to start the paper. I've attached a draft outline. Could you look at the structure and the "decisions" section at the end?

The point-by-point response is in `poster/REVISION4.md`, and the round-4 deck is `XuS_IGVF2026_poster_round4.pptx`.

Best,
Sunan

---

## Slack version

> Round-4 poster update 🧵
> • §5: replaced the MAE-difference panel with "% of loci where CHEWS is closer to truth". With 1 causal SNP, CHEWS wins 75–82% at *every* n; what shrinks with n is the margin, not the frequency. With 20 causal SNPs, GCTA wins 72–77% for n ≥ 300. (A % of h² isn't meaningful there because true h² is pooled over 0.05–0.8.)
> • §7, your 0–0.05 question: with permuted genotypes (true h² = 0), GCTA and CHEWS both average 0.013. In a simulation on the real GTEx genotypes, CHEWS is closer to truth on ~60% of genes in that cell, at every assumed null fraction. But the simulation failed the calibration check I set beforehand (it misses the high-h² corner), so it's not on the poster. A ~20-min re-run with h² up to 0.8 should fix that. Worth doing before print?
> • §5A and §8 unchanged.
> • Also: I drafted a paper outline and would love your eyes on it.
