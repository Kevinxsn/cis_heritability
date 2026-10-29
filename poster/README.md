# IGVF 2026 poster

**`XuS_IGVF2026_poster.pptx`** — 36 in wide × 48 in tall (portrait), one slide.
Open in PowerPoint / Google Slides / Keynote and print at 100%.

> **Revised after PI review** — see [`REVISION.md`](REVISION.md) for the
> point-by-point response. Two claims were withdrawn on the evidence
> (the Simpson's-paradox reading of the heatmap, and "HEELS and GCTA are one
> method"), and the GTEx numbers now supersede the submitted abstract.

Rebuild everything, in this order:

```bash
conda activate dsc180
cd poster
python revision_analyses.py   # tables/ : agreement, disjoint strata, reliability
python make_figures.py        # figures/P0..P8 at true column width
python build_poster.py        # the .pptx + poster_proof.png
python check_numbers.py       # all 76 quoted numbers, re-derived
```

## Round 4 (PI feedback, Oct 2026) — built on top, opt-in

Round 3 above is still the default: the same four commands rebuild the same
deck. Round 4 adds new tables, new figures and a separate deck, selected with
a flag, and never overwrites a round-3 file:

```bash
python revision_analyses.py e f g   # tables/win_by_n_causal, perm_null, pi0_by_stratum,
                                    #   sim_strata_win/_bias, calibration_fit (+ summary json)
python make_figures.py --round 4    # round-3 figures + P3_sample_size_r4, P6_better_or_different
python build_poster.py --round 4    # XuS_IGVF2026_poster_round4.pptx + poster_proof_round4.png
python check_numbers.py --round 4   # every round-3 check + the round-4 pins
```

`f` and `g` need the two round-4 GTEx runs (`results/gtex_perm_split_half`,
`results/gtex_sim_split_half`; see `cish2/scripts/expanse/gtex_{perm,sim}.sbatch`).
On Expanse the whole thing, round-3 rebuild and checksum check included, is
`sbatch cish2/scripts/expanse/poster_round4.sbatch`. Design:
[`../docs/pi_feedback_round4_design.md`](../docs/pi_feedback_round4_design.md);
what was done and the results:
[`../docs/pi_feedback_round4_implementation.md`](../docs/pi_feedback_round4_implementation.md);
response to the PI: [`REVISION4.md`](REVISION4.md).

**Fonts on Expanse.** The figures use Liberation Sans. It was not installed on
Expanse, and matplotlib silently fell back to DejaVu Sans (wider glyphs, so
the cell labels no longer fit). It is now in `~/.local/share/fonts`. If figure
labels look too wide on a new machine, check the build log for `findfont`.

## Files

| file | what it is |
|---|---|
| `XuS_IGVF2026_poster.pptx` | the deliverable — fully editable, native text boxes and tables |
| `poster_proof.png` | flat proof, rendered from the *same* geometry as the deck |
| `make_figures.py` | builds `figures/P0..P8` at true column width |
| `build_poster.py` | content, layout engine, both renderers |
| `check_numbers.py` | re-derives all 76 quoted numbers from the result tables |
| `revision_analyses.py` | agreement stats, disjoint-donor strata, reliability by gene set |
| `REVISION.md` | response to the first PI review |
| `REVISION2.md` | response to the second PI review |
| `REVISION3.md` | response to the pre-PI review (round 3) |
| `figures/*.png,.pdf` | poster figures (PDF versions carry vector text) |

## Design decisions worth knowing

**The figures are not the manuscript figures.** A poster is read from ~1.5 m, so
`make_figures.py` redraws each panel with fewer facets, thicker strokes and
roughly 2× the type size of the versions in `../results/*/figures/`. The
disjoint-strata comparison is stacked rather than side by side, because two 5×5
grids across one column push the cell labels below the legibility floor. Colours
are the validated palette from `cish2/analysis/style.py`, unchanged, so nothing
needs re-validating for colour-vision deficiency.

**Naming.** The fine-mapping estimator is `SuSiE` everywhere in the pipeline and
**CHEWS** on the poster, matching the submitted abstract. The rename happens only
at display time — no result file was touched.

**The h² estimand has two forms, and only one is exact.** The pipeline reports
the trace `Σ(μ²+σ²)`; the exact LD-aware form is
`μᵀRμ − Σₗ mₗᵀR mₗ + Σₗ sₗ`, and the middle (within-credible-set) term is not
optional — dropping it while keeping `μᵀRμ` over-counts badly. Both ship as
`susie_h2` / `susie_h2_quad`; they agree at r = 0.996. See trap 13 in
`../PROJECT.md` and `cish2/tests/test_susie_estimand.R`.

**Win rates are head-to-head against GCTA**, not the pipeline's `susie_win_rate`.
That column is a *three-way* "closest of GCTA / CHEWS / HEELS" share, which is a
different quantity: at one causal SNP it reads 70% where the pairwise number is
76%. The abstract quotes 76%, so the poster uses pairwise throughout and
`check_numbers.py` pins it.

**The layout is computed, not hand-placed.** `build_poster.py` measures every
block, then searches all valid two-cut splits of the ordered block list for the
one that minimises the tallest column. Cuts are forbidden immediately after a
section header (that would orphan the header) and penalised inside a section, so
in practice every column starts with a header. A 0.7 in safety margin is held
back at the foot of each column because the text metric is an estimate of how
Arial wraps, not a measurement of it. If you add text and a column overflows, the
script fails loudly rather than writing a broken deck.

**Fonts.** Body text is Arial (present on any machine that will print this).
Figures use Liberation Sans, which is metrically identical to Arial, so labels
and prose match.

**Figures are drawn at the column width, and this matters.** A figure authored
11 in wide and placed in an 11.07 in column reproduces every rcParam point size
at its true printed size. Draw the same figure 22 in wide and every annotation
halves — which is exactly how the first version ended up with 8 pt cell labels
on a poster meant to be read from 1.5 m. `COLW_IN` / `FULLW_IN` in
`make_figures.py` exist to keep this honest; do not raise a figsize width
without also raising the font sizes.

**Precision means within-locus sampling SD.** The main grid runs one phenotype
per (region, architecture), so its dispersion mixes sampling variability with
between-region LD heterogeneity. `configs/precision.yaml` (via the new
`grid.replicates` axis) draws 25 phenotypes at a fixed region and architecture;
the two agree to <= 0.007 here, but they are different quantities. See trap 14
in `../PROJECT.md`.

## Before printing

- Point the QR code at the real repository and regenerate `figures/qr.png`.
- Fill in the acknowledgements placeholder in §10 — grant numbers. It is marked
  `TODO(Sunan)` in `build_poster.py`.
- Check the affiliation superscripts against the final author list.
- If the venue wants landscape, set `W, H = 48.0, 36.0` at the top of
  `build_poster.py` and re-run; the packer re-flows automatically, though at
  36 in tall you will likely want `NCOL = 4`.

## What is on the poster

| § | content | figure |
|---|---|---|
| header | take-home sentence + what CHEWS computes | `P0_workflow` |
| 1 | the question; accuracy / precision / agreement defined | — |
| 2 | the estimand, and the cross terms the trace form omits | — |
| 3 | benchmark design and fit accounting | — |
| 4 | architecture decides where CHEWS wins | `P2`, `P1_bias_only` |
| 5 | sample size trades bias against error | `P3` |
| 6 | GTEx: agreement (offset, slope, RMSE), not validation | `P5_agreement` |
| 7 | the conditioning artefact, corrected on disjoint donors | `P6_strata_compare` |
| 8 | split-half reliability under three gene-selection rules | — (table) |
| 9–10 | conditional decision guide, limitations, references, QR | — |

Built but held back for the QR supplement: `P4_architecture` (MAF/LD
sensitivity), `P7_reliability_sets`, `P8_estimand`, and the original
`P1_polygenicity` / `P5_gtex_sample_size` / `P6_gtex_architecture`.

Sources: [`../results/FINDINGS.md`](../results/FINDINGS.md),
[`../results/gtex_whole_blood/GTEX_FINDINGS.md`](../results/gtex_whole_blood/GTEX_FINDINGS.md),
[`../results/gtex_split_half/RELIABILITY.md`](../results/gtex_split_half/RELIABILITY.md).
