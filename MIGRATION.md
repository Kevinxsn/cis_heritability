# MIGRATION.md — moving `cis_heritability` to a new server

**Audience:** the coding agent that sets this project up on the new machine.
**Packed:** 2026-10-27, from `mediator.sdsc.edu:/projects/ps-renlab2/sux002/cis_heritability`.
**State at packing:** all analyses finished, IGVF 2026 poster revised three
times (latest on 2026-09-09). Nothing was running.

Read this file first, then [`PROJECT.md`](PROJECT.md) (the science, the
findings, and **15 traps that silently corrupt results**), then
[`cish2/README.md`](cish2/README.md) (the pipeline) and
[`poster/README.md`](poster/README.md) (the poster). Do not edit result
numbers or change estimator code before reading the trap list in `PROJECT.md`.

---

## 1. What the project is, in one paragraph

It compares three estimators of *cis*-heritability (h²) of gene expression:
**GCTA/REML** (individual genotypes, infinitesimal model), **SuSiE** fine-mapping
(summary stats + LD, sparse model; called **CHEWS** on the poster — same thing,
renamed at display time only) and **HEELS** (summary stats + LD, infinitesimal
model). Simulations on 1000 Genomes EUR genotypes measure *bias* (truth known);
GTEx v8 Whole Blood measures *disagreement* and *test-retest reliability* (truth
unknown). All of it runs through one config-driven Python package, `cish2`,
that shells out to GCTA and to R (`susieR`).

---

## 2. Folder map

```
cis_heritability/
├── MIGRATION.md            this file
├── PROJECT.md              START HERE for the science: question, findings, design, traps, open work
├── XuS_IGVF2026 (2).docx   submitted IGVF 2026 abstract
│
├── cish2/                  ALL ACTIVE CODE (~6,400 lines)
│   ├── README.md           pipeline docs, CLI, efficiency notes
│   ├── environment.yml     conda env to rebuild (added for this move)
│   ├── cish2/              the engine  (python -m cish2 plan|run|analyze|gtex)
│   │   ├── cli.py, __main__.py      entry points
│   │   ├── config.py               CampaignConfig dataclasses + defaults (has absolute paths)
│   │   ├── campaign.py, runner.py  simulation driver: regions → cells, multiprocess, resumable
│   │   ├── regions.py              pick ±500 kb windows, LD-architecture deciles
│   │   ├── genotypes.py            direct .bed byte-range reader + .bim index cache
│   │   ├── ld.py, simulate.py, sumstats.py   LD, phenotype simulation, marginal OLS
│   │   ├── estimators/             gcta.py (GRM + REML), susie.py (→ R), heels.py (own impl.)
│   │   ├── gtex.py, gtex_campaign.py, gtex_runner.py   real-data pass (GTEx)
│   │   ├── results.py              manifest.jsonl + sharded cells.tsv.gz writer
│   │   └── scratch.py              picks fast local scratch for LD/GRM files
│   ├── analysis/           metrics (with Monte-Carlo SEs), figures, REPORT/FINDINGS/OVERVIEW generators,
│   │                       gtex_report.py, reliability.py, compare_legacy.py
│   ├── configs/*.yaml      one file per campaign (see §5)
│   ├── R/                  susie_batch.R (simulation), susie_genes.R (GTEx)
│   ├── scripts/            env.sh (source it), run_campaign.sh, analyze_all.sh,
│   │                       relocate_paths.sh (added for this move)
│   └── tests/              test_pipeline.py (pytest, 15 checks), test_susie_estimand.R
│
├── results/                ALL OUTPUTS — finished, keep as-is
│   ├── FINDINGS.md / .pdf  simulation write-up
│   ├── OVERVIEW.md, overview.csv, conditions.csv   13 simulated conditions in one table
│   ├── figures/            cross-condition figures
│   ├── <campaign>/         baseline_n300, sample_size, effect_size, maf_arch, high_ld, low_ld,
│   │                       precision, estimand_check, smoke — each has REPORT.md, tables/,
│   │                       figures/, cells.tsv.gz, shards/, manifest.jsonl, config.yaml, logs/
│   ├── gtex_whole_blood/   GTEX_FINDINGS.md + real-data results (97,160 fits)
│   └── gtex_split_half/    RELIABILITY.md — test-retest on disjoint donor halves
│
├── poster/                 IGVF 2026 poster, 36×48 in portrait — LATEST RESULTS LIVE HERE
│   ├── XuS_IGVF2026_poster.pptx   the deliverable
│   ├── poster_proof.png           flat proof
│   ├── README.md                  rebuild order + design decisions
│   ├── revision_analyses.py → tables/   make_figures.py → figures/P0..P8
│   ├── build_poster.py            layout engine, writes .pptx + proof
│   ├── check_numbers.py           re-derives all 76 quoted numbers from results/
│   └── REVISION*.md, *_feedback*.md   three rounds of PI review and responses
│
├── expand_proposal/        ACCESS compute-allocation proposal + benchmarks (standalone, not needed to run)
├── simulation_gtex/        SUPERSEDED old scripts — ignore
└── simulations/            SUPERSEDED old scripts — ignore, BUT keep
                            result_eur_num_people_100_300_regular_LD.txt
                            (read by analysis/compare_legacy.py)
```

### Not in the zip (on purpose)

| left behind | size | why / how to get it |
|---|---|---|
| input data `/projects/ps-renlab2/sux002/data/` | 14 GB | see §3 — copy separately |
| `.git/` | 404 MB | stale (3 commits, most current files untracked). `git init` fresh on the new server |
| `simulations/temp_ld.ld` | 2.3 GB | scratch LD matrix from the old pipeline |
| `results/.cache/` | 169 MB | `.bim` position index, rebuilt automatically on first run |
| `__pycache__/` | — | regenerated |

---

## 3. Input data (not in the zip — must be staged separately)

The code expects one data root (old: `/projects/ps-renlab2/sux002/data`) with:

```
<DATA>/1000G.EUR.{1..22}.{bed,bim,fam}       ~13.7 GB  1000 Genomes Phase 3 EUR (503 people), simulations
<DATA>/gtex/GTEx_v8_genotype_EUR_HM3_exclude_dups.allchr.reorder.{bed,bim,fam}   GTEx v8 EUR, HapMap3 SNPs
<DATA>/gtex/Whole_Blood.v8.normalized_expression.bed.gz
<DATA>/gtex/Covar_all_Whole_Blood.txt
```

- **1000G EUR** is public. Either `rsync` it from the old server or re-download
  the 1000 Genomes Phase 3 EUR PLINK files (the same set distributed with LDSC
  / the Alkes Price lab). Exact file identity matters only if you want to
  reproduce the existing simulation numbers bit-for-bit (region sampling is
  seeded against the `.bim`).
- **GTEx v8 genotypes are dbGaP controlled-access.** Only move them to a
  machine covered by the project's data-use approval, server-to-server
  (`rsync -av old:/projects/ps-renlab2/sux002/data/gtex/ new:<DATA>/gtex/`),
  not via a laptop. The expression and covariate files are public GTEx portal
  downloads.
- **None of the data is needed** to read results, rebuild the poster, or run
  `poster/check_numbers.py` — those work from `results/` alone.

---

## 4. Deploy on the new server

### 4.1 Unpack

```bash
mkdir -p <NEW_ROOT_PARENT> && cd <NEW_ROOT_PARENT>
unzip cis_heritability_2026-10-27.zip          # creates ./cis_heritability/
export PROJ=$PWD/cis_heritability
```

### 4.2 Environment

The old env was `/projects/ps-renlab2/sunan/anaconda3/envs/dsc180` (Python 3.8.16,
pandas 2.0.3, numpy 1.24.3, R 4.2.3, susieR 0.12.35, data.table 1.15.4,
Rfast 2.1.0, GCTA 1.94.1, PLINK 1.9, python-pptx 1.0.2). Rebuild it:

```bash
cd $PROJ/cish2
mamba env create -f environment.yml     # or conda; env name: cish2
export CISH2_ENV="$(conda info --base)/envs/cish2"
```

`scripts/env.sh` prepends `$CISH2_ENV/bin` to PATH, pins BLAS to 1 thread and
sets PYTHONPATH. It defaults `CISH2_ENV` to the **old** path — either export
`CISH2_ENV` in your shell profile or edit the default on line 11 of
`scripts/env.sh`. Required on PATH after sourcing: `python`, `Rscript`,
`gcta64`, `plink`.

Constraints inherited from the old env (keep them unless you also re-test):
no `pyarrow` (tables are gzipped TSV), no `jsonlite` in R (matrices cross into R
as raw float64 binary).

### 4.3 Rewrite absolute paths

Code and configs hard-code the old project and data roots (`configs/*.yaml`,
`cish2/config.py`, `cish2/gtex_campaign.py`, `analysis/{overview,findings,
reliability,gtex_report}.py`, `scripts/analyze_all.sh`, `tests/test_pipeline.py`).
One script fixes all of them:

```bash
cd $PROJ/cish2
scripts/relocate_paths.sh $PROJ <NEW_DATA_ROOT> --dry-run   # list files
scripts/relocate_paths.sh $PROJ <NEW_DATA_ROOT>
```

It leaves `results/*/config.yaml` and `results/*/logs/` untouched — those are
records of past runs. Note: the per-row `fingerprint` column hashes the config
including the data path, so rows produced after the move will carry a different
fingerprint string. Nothing in `analysis/` filters on it; resume is keyed on
the region slot in `manifest.jsonl`, so resuming old campaigns still works.

### 4.4 Scratch and concurrency — tune to the new machine

- **Scratch.** Each worker writes a p×p float64 LD matrix (~620 MB at p≈8.8k)
  and a binary GRM. `cish2/scratch.py` picks the first writable dir with
  ≥ 50 GB free from: `$SLURM_TMPDIR`, `/scratch/$USER/$SLURM_JOB_ID`,
  `/scratch/$USER`, `/data/tmp`, `/tmp`. On the old box shared NFS was 13×
  slower than local NVMe and moving scratch to NVMe halved wall time. Make sure
  a **node-local** disk wins, or force it: `-o compute.scratch_root=/fast/disk`
  (simulations) / `scratch_root:` key (GTEx configs).
- **Workers.** Configs say `workers: 48` (a polite share of a 256-core box with
  no scheduler). Override per run: `-o compute.workers=N` (sim) or
  `-o workers=N` (GTEx). Memory: budget ~2–3 GB per worker.
- **If the new server has Slurm**, wrap `python -m cish2 run ...` in an
  `sbatch` with `--cpus-per-task=N` and set workers to N; `scratch.py` already
  prefers `$SLURM_TMPDIR`. Without a scheduler, launch with `nohup` as before.

### 4.5 Verify, in this order (stop at the first failure)

```bash
cd $PROJ/cish2 && source scripts/env.sh
which python Rscript gcta64 plink

# a) no data needed — results and poster are intact
cd $PROJ/poster && python check_numbers.py       # expect: 76 checks, all passing

# b) needs 1000G data — code + tools work
cd $PROJ/cish2 && python -m pytest tests/test_pipeline.py -q   # expect: 15 passed
Rscript tests/test_susie_estimand.R

# c) tiny end-to-end run (writes results/smoke/)
python -m cish2 plan configs/smoke.yaml
python -m cish2 run  configs/smoke.yaml -o compute.workers=4
python -m cish2 analyze ../results/smoke
```

If (a) passes and (b) fails, it is environment or data paths, not science.

---

## 5. Campaigns (all finished — do not rerun unless asked)

| config | result dir | what it varies | approx. cost (48 workers) |
|---|---|---|---|
| `baseline.yaml` | `baseline_n300` | 5 causal counts × 6 h², n=300 | ~23 min |
| `sample_size.yaml` | `sample_size` | n = 100…1000 | |
| `effect_size.yaml` | `effect_size` | effect concentration (normal vs heavy-tailed) | |
| `maf_arch.yaml` | `maf_arch` | causal MAF bin | |
| `high_ld.yaml`, `low_ld.yaml` | `high_ld`, `low_ld` | top/bottom LD-decile regions | |
| `precision.yaml` | `precision` | 25 replicate phenotypes per locus (within-locus SD) | |
| `estimand_check.yaml` | `estimand_check` | exact vs trace SuSiE h² form | |
| `smoke.yaml` | `smoke` | tiny sanity run | ~1 min |
| `gtex_whole_blood.yaml` | `gtex_whole_blood` | nested n = 100⊂200⊂300⊂450⊂574, ~19.2k genes | ~17 min |
| `gtex_split_half.yaml` | `gtex_split_half` | disjoint donor halves (reliability) | |

Whole simulation study ≈ 101 core-hours; GTEx ≈ 13 core-hours. Every campaign
is resumable and idempotent. Regenerate all reports with `scripts/analyze_all.sh`.

---

## 6. Where things stand — latest results (from the poster, 2026-09-09)

Headline numbers, all pinned by `poster/check_numbers.py`:

- **Polygenicity is the switch.** At 1 causal SNP SuSiE/CHEWS is nearly
  unbiased (−0.008) and beats GCTA on **76%** of loci (pairwise win rate). At 20
  causal SNPs both CHEWS forms underestimate h² by ≈0.13; GCTA and HEELS stay
  near zero. Attribution to the sparse-effect model is **provisional** until the
  matched-LD factorial runs.
- **Sample size trades bias for error.** CHEWS bias −0.122 (n=100) → −0.021
  (n=1000), yet at n=100 its MAE is ~31% lower than GCTA's.
- **Precision.** Within-locus sampling SD is 2–4× lower for CHEWS than GCTA;
  within- and across-locus SDs agree to ≤ 0.007.
- **GTEx agreement** (`poster/tables/gtex_agreement.csv`): CHEWS−GCTA mean
  offset −0.015 (n=100) → −0.007 (n=574); HEELS tracks GCTA at r = 0.997–1.000.
- **Split-half reliability** (`poster/tables/reliability_by_gene_set.csv`):
  CHEWS had the highest split-half correlation (0.91 vs GCTA 0.86, all genes)
  — one donor split only.
- **Withdrawn after PI review:** the Simpson's-paradox reading of the GTEx
  heatmap (now shown on disjoint donors), and "HEELS and GCTA are one method".

### Open work (from `poster/REVISION3.md` and `PROJECT.md`)

1. **QR URL + grant numbers** on the poster — user must supply
   (`TODO(Sunan)` in `poster/build_poster.py`; set `QR_URL`, regenerate
   `figures/qr.png`).
2. **Matched-LD factorial** — CHEWS with in-sample LD
   (`-o method.ld_for_susie=sim`), HEELS with panel LD. The single most
   valuable missing run (~25 min).
3. **Exact SuSiE estimand** (`susie_h2_quad`) through the sample-size and GTEx
   passes (~4 h).
4. **Repeated donor splits** for the reliability result.
5. `-o method.legacy_ambiguous_filter=false` rerun (removes a ~15% SNP handicap).
6. Larger: more GTEx tissues, denser-than-HapMap3 genotypes, a calibrated SuSiE SE.

---

## 7. Rules for the agent

- Read the trap list in `PROJECT.md` before touching covariates, LD,
  matrix I/O to R, or `susie_rss` arguments. Several traps produce plausible
  but wrong numbers with no error.
- Never pass `check_prior = FALSE` to `susie_rss`. Never feed HEELS
  reference-panel LD. Give GCTA covariates via `--qcovar`, never a
  pre-residualised phenotype.
- Always report MAE beside RMSE; "SuSiE beats GCTA" must use
  `pairwise_win_rate`, not the three-way `susie_win_rate`.
- Do not delete `simulations/result_eur_num_people_100_300_regular_LD.txt`.
- Keep large transient files on node-local scratch, not on network storage.
- After any change that could move a number, rerun `poster/check_numbers.py`.
