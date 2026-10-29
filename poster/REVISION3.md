# Response to the pre-PI review (round 3)

Point-by-point against `XuS_IGVF2026_poster_pre_PI_feedback_round3.md`.
Earlier rounds: [`REVISION.md`](REVISION.md), [`REVISION2.md`](REVISION2.md).
All numbers below are re-derived by `check_numbers.py` (76 checks, all passing).

---

## The five pre-send changes

### 1. The Section 4 sentence — **it was wrong, and the review is right**

I had written that GCTA, HEELS *and the LD-aware estimator* "all stay put." The
LD-aware estimator does not stay put: it runs from −0.009 at one causal SNP to
−0.130 at twenty, lying on the CHEWS curve. I meant "tracks the trace form" and
wrote something that contradicted my own figure.

Replaced with the reviewer's construction:

> One causal SNP: CHEWS nearly unbiased (−0.008), closer on 76% of loci. At
> twenty, **both** CHEWS forms underestimate h² by ≈0.13 while GCTA and HEELS
> stay near zero, so the LD cross terms do not explain it. **Consistent with the
> sparse-effect model** — but CHEWS and HEELS also differ in LD source (§5), so
> that attribution is provisional.

Both endpoints are now pinned in `check_numbers.py`, so that sentence cannot
drift from the figure again.

The causal language is downgraded throughout: "the sparse-effect assumption
itself" → "consistent with the sparse-effect model … provisional". The matched
factorial (CHEWS with in-sample LD; HEELS with panel LD) has not been run, and
the poster says the attribution waits on it.

### 2. Which CHEWS form is in which figure — **scoped, with the tail reported**

"The two forms agree" was too broad. The distribution of |exact − trace| over
3,000 fits at n = 300:

| median | p90 | p95 | p99 | max |
|---:|---:|---:|---:|---:|
| 0.002 | 0.019 | 0.034 | 0.089 | **0.368** |

§2 now reads: at n = 300 the two forms track each other on average (r = 0.996,
median 0.002, mean bias curves overlap at every causal count), **per locus they
can differ by up to 0.368**, and **the exact form was evaluated at n = 300 only;
§5–6 use the trace form.** That is the review's option 3, stated rather than
generalized. Propagating the exact form through the sample-size and GTEx passes
is listed as outstanding work below.

### 3. Every symbol defined; one notation

§2 was rewritten in the reviewer's order — define the decomposition first, then
give the estimator — rather than starting from an approximation and retracting
it:

> Truth is the quadratic form **βᵀRβ**. SuSiE writes β = Σₗ βₗ, one single-effect
> vector per component, each selecting exactly one variant. Let
> **mₗ = E[βₗ] = αₗ⊙μₗ** and **sₗ = E[βₗᵀRβₗ]**. Then
> **E[βᵀRβ | data] = μᵀRμ − Σₗ mₗᵀR mₗ + Σₗ sₗ**.

`sₗ = E[βₗᵀRβₗ]` is both the definition the reviewer proposed and exactly what
the code computes: one single effect has a single non-zero entry, so
`βₗᵀRβₗ = R_jj βₗⱼ² = βₗⱼ²` under a unit diagonal, whose expectation is
`Σⱼ αₗⱼ(μₗⱼ² + σₗⱼ²)`. `β` is used throughout, including the brute-force check.

### 4. Precision — **ran the analysis instead of rewording it**

The review was right that SD across regions conflates sampling variability with
between-region LD heterogeneity, because the main grid has one replicate per
region × parameter set. The reviewer offered a wording fallback; the pipeline
turned out to be one small change from doing it properly.

`GridConfig` gained a `replicates` axis (default 1, so nothing else changes) and
`configs/precision.yaml` runs **30 regions × {n = 100, 300} × {1, 20 causal} ×
h² = 0.2 × 25 independent phenotype draws** — 3,000 fits, region and architecture
held fixed. That makes a conditional sampling SD estimable per region.

**Within-locus sampling SD**, averaged over the 30 regions:

| n | causal | GCTA | CHEWS | HEELS | GCTA/CHEWS |
|---:|---:|---:|---:|---:|---:|
| 100 | 1 | 0.246 | **0.108** | 0.398 | 2.3× |
| 100 | 20 | 0.231 | **0.054** | 0.382 | 4.3× |
| 300 | 1 | 0.126 | **0.045** | 0.123 | 2.8× |
| 300 | 20 | 0.089 | **0.043** | 0.090 | 2.1× |

Two things follow. First, the precision claim is now supported by a clean
conditional design, and the poster states the measured range (2–4×) rather than a
single pooled ratio. Second — and this is the useful part — **the within-locus
and across-locus SDs agree to ≤ 0.007 in every condition tested**, so the
confounding the review worried about is empirically negligible here. Panel A can
keep plotting across-locus dispersion over the full *n* range, and the caption
says the within-locus design agrees to ≤ 0.007. Logged as trap 14 in
`PROJECT.md`, with the caveat that this is an empirical finding, not a licence to
conflate the two quantities elsewhere.

### 5. Footer placeholders — **removed**

`[repository URL here]` and `[funding / grant numbers]` are gone. The footer now
reads *Amariuta Lab, UC San Diego · Contact: sux002@ucsd.edu*, and the URL line
is suppressed while there is no destination. `build_poster.py` prints a loud
warning on every build until `QR_URL` is set:

```
  !! QR code still points at a placeholder (https://github.com/).
     Set QR_URL and QR_URL_TEXT, regenerate figures/qr.png, and
     re-run before printing. The URL line is hidden until then.
```

**This is the one thing left that only you can supply.** If there is no public
repository by print time, delete the QR block rather than shipping a dead code.

---

## Scientific points the review anticipated

**"Why does HEELS fail at n = 100?"** — the phrase "breaks down" was doing no
work, and the answer is more interesting than the phrase. At n = 100
(p/n ≈ 70, p measured after the allele filter) HEELS **converges on 100% of
fits** with no boundary flags — but **46% of estimates exceed ĥ² = 0.9**, mean
bias +0.17. It saturates rather than failing. The poster says exactly that.

**"Win rate or average error?"** — both are kept, and the distinction is stated
once in the §4 caption. Checked whether rare catastrophic fits drive the means:
they do not. At every causal count the median absolute error tracks the mean
(e.g. at 20 causal, CHEWS mean 0.135 / median 0.117; GCTA 0.069 / 0.059).

**"What do the intervals represent?"** — captions now say *region-clustered
bootstrap, 2,000 resamples*.

**"What does the gene bootstrap support?"** — §8 now leads with the limitation
rather than the interval, in the reviewer's words: *in this one donor split
CHEWS had the higher split-half correlation on all three gene sets; repeated
donor partitions are needed to quantify robustness to the split itself.*

---

## Communication and polish

- Neutral titles: "Where CHEWS beats GCTA" → **"Fraction of simulations with
  lower CHEWS absolute error"**; "Where the CHEWS shortfall comes from" →
  **"CHEWS bias by causal-effect count"**; the §5 panels are "Dispersion of ĥ²"
  and "Mean absolute error difference".
- Banner uses **performance**, not choice: *Estimator performance depends jointly
  on regulatory architecture, sample size and LD source.*
- §5 panel B row labels enlarged to 22 pt bold; colorbar reads
  **MAE_CHEWS − MAE_GCTA (negative favors CHEWS)**.
- Word count **1,012 → 964**.

## Still outstanding

1. **QR destination and grant numbers** — yours to supply.
2. **Matched-LD factorial** (CHEWS in-sample × HEELS panel). Until it runs, the
   sparsity attribution stays provisional on the poster.
3. **Exact estimator through the sample-size and GTEx passes** (~4 h).
4. **Repeated donor splits** for the reliability result.
5. **MAF result on discrete/haplotype genotypes.**

## Verification

- `check_numbers.py`: 76 checks, all passing.
- `pytest tests/test_pipeline.py`: 15 passed, including the brute-force estimand
  identity.
- Deck: 36 × 48 in, nothing off-sheet, smallest type 20 pt, no literal TeX.
