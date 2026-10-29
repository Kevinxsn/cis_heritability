# PI feedback, round 4: implementation record and results

**Done:** 2026-10-28 on SDSC Expanse (account `sdp202`) · **Design:**
[`pi_feedback_round4_design.md`](pi_feedback_round4_design.md) · **Response to the PI:**
[`../poster/REVISION4.md`](../poster/REVISION4.md)

## TL;DR

- **Item A (§5 panel B) is done as designed.** It is now a win-rate heatmap on
  the §4 colour scale. CHEWS is closer on 75–82% of loci at one causal SNP for
  every n; GCTA is closer on 72–77% at 20 causal SNPs for n ≥ 300. Panel A is
  pixel-identical to round 3.
- **Item B (§7) ran in full, and the simulation failed its pre-registered
  calibration rule.** The best mix missed the real gap means by 0.0139 against
  a 0.010 limit. Per design §3.3 row 1, §7 keeps its round-3 figure. The
  model-free permutation null did go on the poster: at true h² = 0, constrained
  GCTA and CHEWS both average 0.013.
- **Exploratory answer to the PI (not on the poster):** in the 0-signal,
  h² < 0.05 cell, CHEWS is closer to truth on 60% (59–61%) of simulated genes,
  and above 50% at every null fraction from 0 to 0.9.
- **Round 3 still reproduces.** Every round-4 change is opt-in, and the default
  commands give the round-3 tables byte for byte. Figures are byte-identical
  between pre-edit and post-edit code on this machine, all 76 round-3 number
  checks pass, and real-mode GTEx reproduces the archived split-half rows to 1e-8.

---

## 1. How to reproduce

### Round 3 (unchanged commands, unchanged results)

```bash
cd poster
python revision_analyses.py        # parts A–D only, as before
python make_figures.py             # round-3 figures only
python build_poster.py             # XuS_IGVF2026_poster.pptx + poster_proof.png
python check_numbers.py            # 76 checks
```

### Round 4 (on top)

```bash
# GTEx passes (Expanse; resumable)
sbatch cish2/scripts/expanse/round4_tests.sbatch   # pytest (20) + 200-gene smoke of each mode
sbatch cish2/scripts/expanse/gtex_perm.sbatch      # -> results/gtex_perm_split_half  (~7 min, 32 cores)
sbatch cish2/scripts/expanse/gtex_sim.sbatch       # -> results/gtex_sim_split_half   (~17 min, 40 cores)

# poster: round-3 rebuild + checksum check, then round 4
sbatch cish2/scripts/expanse/poster_round4.sbatch
#   = revision_analyses.py e f g ; make_figures.py --round 4 ;
#     build_poster.py --round 4 ; check_numbers.py --round 4
```

On Expanse, everything goes through `sbatch`: `import pandas` and
`import matplotlib` hang on the login node. The `gtex_sim.sbatch` header asks
for `shared`; I submitted it with `-p debug -c 40 --mem=60G` because the debug
node was free sooner.

### Safety net

- `/home/sux002/work/cis_heritability_round3_snapshot_2026-10-28.tar.gz`:
  `poster/`, `cish2/`, `docs/`, `PROJECT.md`, `MIGRATION.md` as they were
  before any round-4 edit. It includes the mediator-rendered round-3 figures
  and deck.
- `logs/round4/round3_output_md5.txt`: checksums of those original outputs.
  `logs/round4/round3_output_md5_expanse.txt` holds the same files as rendered
  on Expanse (see §3, fonts).

---

## 2. What changed

Everything is additive. Real mode and the round-3 build take exactly the paths
they always did.

| file | change |
|---|---|
| `cish2/cish2/simulate.py` | `sim_beta_empirical(Z, k, h2, rng)`: effects on real genotypes, rescaled so `var(Zβ) = h²` exactly. The effect draw moved to a shared `_draw_effects` making the identical RNG calls, so `sim_beta` is bit-identical (tested). |
| `cish2/cish2/gtex_campaign.py` | New keys, all default-off: `phenotype_mode` (`real`/`simulated`/`permuted`), `sim_design`, `gcta_also_unconstrained`; validated. |
| `cish2/cish2/gtex_runner.py` | Real mode is the original loop, untouched. Simulated: K draws per gene from the design prior, one β shared by both halves, extra columns `phenotype_mode, draw, h2_target, ncausal_target, n_causal_snps, sim_beta_hash, h2_sample, h2_resid`. Permuted: genotype rows permuted within each half via a new `_prepare_cell(geno_idx=…)`. Optional second, unconstrained REML (`gcta_h2_uc`, …). RNG = `(seed, crc32(gene_id), draw)`. |
| `cish2/configs/gtex_{sim,perm}_split_half.yaml` | Copies of `gtex_split_half.yaml` (same seed, so same halves) plus the mode. |
| `cish2/tests/test_pipeline.py` | +5 tests (below). |
| `cish2/scripts/expanse/{round4_tests,gtex_sim,gtex_perm,poster_round4}.sbatch` | New job scripts. |
| `cish2/README.md`, `poster/README.md` | Additive "round 4" sections. |
| `poster/revision_analyses.py` | `part_e` (§5 win rates), `part_f` (permutation null, π̂₀), `part_g` (simulation + calibration, decision rules in its docstring). No argument runs A–D as before; `e f g` are opt-in. |
| `poster/make_figures.py` | `WIN_CMAP/WIN_VMIN/WIN_VMAX` shared by P2 and the new panels. §5 panels factored into helpers used by both rounds. New `fig_sample_size_winrate` → `P3_sample_size_r4.png` (+ rejected `…_A3_rejected.png`) and `fig_better_or_different` → `P6_better_or_different.png` (stamped EXPLORATORY). New figures only with `--round 4`. |
| `poster/build_poster.py` | `--round 4` takes the round-3 block list and swaps only §5 (figure, caption, text), §7 (text) and adds one §10 limitation; writes `*_round4.pptx/png`. Prints a word count. |
| `poster/round4_text.py` | Hand-typed round-4 §7 text, chosen by the calibration verdict, plus its number pins. |
| `poster/check_numbers.py` | `--round 4` runs every round-3 check plus 15 new pins. (I briefly broke the entry point while editing; fixed and re-verified: 76/76.) |
| `poster/REVISION4.md` | Response to the PI. |
| `~/.local/share/fonts/LiberationSans-*.ttf` | Font the figures require; it was missing on Expanse (§3). |

New result directories: `results/gtex_perm_split_half` (38,864 cells),
`results/gtex_sim_split_half` (77,728 cells), and the smoke runs
`results/gtex_{sim,perm}_smoke`. New tables are in `poster/tables/`:
`win_by_n_causal, perm_null, pi0_by_stratum, calibration_fit,
calibration_summary.json, sim_strata_win, sim_strata_bias, sim_win_vs_pi0`.

---

## 3. Verification

| check | result |
|---|---|
| `pytest tests/test_pipeline.py` | **20 passed** (15 old + 5 new) |
| new: `sim_beta` refactor bit-identical to the old code | pass (3 effect laws × 5 seeds) |
| new: `sim_beta_empirical` gives `var(Zβ) = 0.3` to 1e-10 with exactly k non-zeros | pass |
| new: simulated mode uses the same β in both halves; deterministic across workers | pass |
| new: permuted mode gives bit-identical phenotypes, only genotype rows reordered | pass |
| new: **real mode on chunk 0 vs archived `results/gtex_split_half`, 11 columns, 1e-8** | **pass** |
| round-3 tables (6 CSVs) after rebuild vs pre-round-4 checksums | **identical** |
| round-3 figures: snapshot code vs current code, same machine | **13/13 byte-identical** |
| round-3 proof PNG: snapshot code vs current code | byte-identical; `.pptx` members all identical (zip timestamps differ) |
| §5 panel A, round 4 vs round 3 | left 40% of the figure: 0 differing pixels |
| `check_numbers.py` / `--round 4` | 76/76 · 76 + 15 = 91/91 |
| layout | both decks fit; round 4 tightest slack +0.05 in (after the 0.7 in safety margin) |

**Why rebuilt round-3 figures don't match the mediator checksums.** The PNGs
rendered on Expanse differ from the archived ones even with the untouched
snapshot code: same pixel size and content, different anti-aliasing. That
comes from the matplotlib/FreeType build, not from round 4. One real
environment gap turned up: **Liberation Sans was not installed**, so matplotlib
silently used DejaVu Sans, which is wider and pushed cell labels past their
cells. I installed Liberation Sans 2.1.5 (SIL OFL) into `~/.local/share/fonts`
and cleared the matplotlib font cache. Note that the pass-1 rebuild overwrote
`poster/figures/*.png`, the round-3 deck and the proof with the Expanse
renders; the mediator originals are in the snapshot tarball.

---

## 4. Results

### A. §5: how often CHEWS wins (`tables/win_by_n_causal.csv`)

| causal ↓ / n → | 100 | 300 | 500 | 800 | 1000 |
|---|---:|---:|---:|---:|---:|
| **win rate** | | | | | |
| 1 | 82% | 81% | 78% | 78% | 75% |
| 3 | 65% | 64% | 62% | 62% | 64% |
| 5 | 64% | 57% | 54% | 61% | 57% |
| 10 | 53% | 41% | 37% | 40% | 43% |
| 20 | 40% | 28% | 24% | 23% | 24% |
| **MAE_C / MAE_G − 1** | | | | | |
| 1 | −73% | −65% | −61% | −57% | −58% |
| 20 | +24% | +91% | +112% | +94% | +91% |

These match the design doc's quick check exactly. 290–300 paired fits per
cell; region-clustered CIs exclude 50% in 23 of 25 cells.

### B2. Permutation null (`tables/perm_null.csv`, `pi0_by_stratum.csv`)

38,864 cells, 0 failed chunks. At true h² = 0 on real expression:

| | mean | median | p90 | at floor |
|---|---:|---:|---:|---:|
| GCTA (constrained) | 0.0132 | 0.000 | 0.045 | 57% at bound |
| CHEWS | 0.0129 | 0.009 | 0.032 | 40% < 0.005 |
| HEELS | 0.0110 | 0.000 | 0.037 | 64% at bound |

Within the permuted 0-signal, h² < 0.05 stratum the picture is the same
(0.013 vs 0.013).

Null fraction π̂₀ = P_real(REML at bound) / P_perm(REML at bound), with real
strata from half A and p-values from half B:

| signals \ h² | 0–.05 | .05–.1 | .1–.2 | .2–.4 |
|---|---:|---:|---:|---:|
| 0 | **0.73** | 0.50 | 0.49 | (6 genes) |
| 1 | 0.49 | 0.29 | 0.16 | 0.03 |
| 2 | 0.43 | 0.30 | 0.10 | 0.02 |

Genome-wide 0.50; 0.55 (p > 0.25) and 0.62 (p > 0.1) as sensitivity.

### B1 + B3. Semi-synthetic GTEx, calibrated (`tables/sim_strata_*`, `calibration_*`)

77,728 cells (19,432 genes × 2 draws × 2 halves), 0 failed. Constrained GCTA
returned on 100%, unconstrained on 91.9% (vs 25–40% on real expression).
`h2_resid` matches `h2_target` to < 0.002 at every h², so the covariate PCs
absorb almost no simulated signal.

**Calibration (480 priors).** Best fit π₀ = 0.1, h2_scale = 0.5, q = 1
(uniform over causal counts). S = 203 over 36 terms; 2 priors acceptable,
both with π₀ = 0.1.

| rule-1 criterion (fixed before the run) | value | limit | |
|---|---:|---:|---|
| sign agreement of gap means | 94% | ≥ 80% | ok |
| mean \|sim − real\| gap | **0.0139** | ≤ 0.010 | **FAIL** |
| TVD of stratum shares | 0.022 | ≤ 0.15 | ok |

→ **Rule 1 fails: B1 is not on the poster.** The misfit sits in the
high-h² / many-signal cells (e.g. 1 signal, .4–1: simulated −0.051 vs real
+0.003; 4+ signals, .4–1: −0.084 vs −0.031). In the 0–0.05 and .05–.1 columns
the simulated gaps are within 0.004 of the real ones. Likely cause: the design
grid's maximum h² is 0.4, so the simulation can't produce the genes that
populate the real .4–1 stratum. Calibration also prefers π₀ = 0.1, which
contradicts B2's ~0.5–0.73. Both point to the simulated architecture, not the
estimators.

**Exploratory per-stratum results (best prior).** % of genes where CHEWS is
closer to truth; † = acceptable priors or CIs straddle 50%.

| signals \ h² | 0–.05 | .05–.1 | .1–.2 | .2–.4 | .4–1 |
|---|---:|---:|---:|---:|---:|
| 0 | **60%** | 59% | 56% | — | — |
| 1 | **63%** | 61% | 57% | 60% | 68% |
| 2 | 67% | 61% | 53%† | 54% | 59% |
| 3 | — | — | 47%† | 47%† | 57% |
| 4+ | — | — | — | 42%† | **36%** (GCTA better) |

In the 0-signal, 0–0.05 cell (truth mean 0.016): GCTA bias +0.009, CHEWS
+0.005, unconstrained GCTA +0.004. MAE is 0.016 for CHEWS vs 0.022 for GCTA.
CHEWS beats unconstrained GCTA on 77%. Win rate against π₀ (other knobs at
best fit): 61% at π₀ = 0 down to 54% at π₀ = 0.9, CI above 50% throughout;
~55% at the B2-anchored π₀ ≈ 0.7. Using `h2_resid` as truth moves the 0–0.05
column by ≤ 1.3 points, and any cell by ≤ 4.4 (the high-h² cells); no cell
changes side of 50%.

---

## 5. Deviations from the design, and why

1. **π̂₀ event.** The design's λ = 0.5 is degenerate: GCTA reports p = 0.5
   exactly at the bound and never above it, so P(p > 0.5) = 0 in both runs. I
   used "REML at the bound" as the primary event and report λ = 0.25 and 0.1.
   This is documented in `part_f`'s docstring.
2. **§5 shipped A2, not A3.** The two-line labels collided at the required
   ≥ 16 pt. The relative-MAE numbers stayed in the body text (−73%/−58%); the
   caption has no room for the full range after layout trimming.
3. **§7 not redesigned.** That follows from the pre-registered rule 1, not a
   choice made after seeing results. `P6_better_or_different.png` is built but
   stamped EXPLORATORY and is not in the deck.
4. **Rule 1 had to be made concrete.** The design said "gap means badly
   missed". I fixed MAD ≤ 0.010, TVD ≤ 0.15, the S-score and the
   acceptable-set tolerance (χ²₃,₉₅ × reduced χ²) in `part_g`'s docstring
   before the simulation finished. The 0.010 limit is the one that failed.
5. **§7 coupling paragraph trimmed to one sentence** (as the design
   suggested), plus one permutation-null sentence. Word count 910 vs 864
   (my counter: body, bullets, callouts, captions, tables; the design's "964"
   counted differently).
6. **P2's colour-bar label not changed.** Same map and limits as the new panels
   (read from shared constants); only the label wording differs, so the
   round-3 P2 stays byte-identical.
7. **Panel A geometry frozen.** The round-4 P3 lays out the round-3 figure
   first, then swaps panel B into the same axes. Without that, `tight_layout`
   moved panel A by a few pixels.
8. **Job sizing.** The design suggested a full node for 1 h. The debug
   partition was faster to start: perm 32 cores, 7 min; sim 40 cores, 17 min.
   Total round-4 compute ≈ 20 core-hours (~20 SU).

---

## 6. Open items: need your (or the PI's) call

1. **Re-run B1 with the h² grid extended to 0.8?** About 17 min. It targets
   the identified misfit, but changing the design after a failed check must be
   reported as such; keep the same rule-1 thresholds.
2. **Is the §7 permutation sentence the right claim for the poster?** Removing
   it is a one-line edit in `poster/round4_text.py`.
3. The design's open questions 1–5 still stand (which "bottom-left" cell; poster
   vs slides; dropping the coupled panel; acceptable prior dependence;
   unconstrained GCTA as a series).
4. Before printing (unchanged from round 3): QR URL and grant numbers
   (`TODO(Sunan)` in `build_poster.py`).
