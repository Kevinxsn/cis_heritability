"""Self-checks for the parts of the pipeline that are easy to get silently wrong.

Run with::

    python -m pytest tests -q          # or:  python tests/test_pipeline.py
"""

from __future__ import annotations

import os
import subprocess
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from cish2.estimators.gcta import parse_hsq, write_grm_binary  # noqa: E402
from cish2.estimators.heels import HeelsSolver, heels_reference  # noqa: E402
from cish2.estimators.susie import _write_cells  # noqa: E402
from cish2.genotypes import read_bed_slice, read_bim_rows, build_chrom_index  # noqa: E402
from cish2.ld import build_ld, compute_s2g  # noqa: E402
from cish2.simulate import sim_beta, sim_geno, sim_trait  # noqa: E402
from cish2.sumstats import marginal_ols  # noqa: E402

GENO = "/expanse/lustre/projects/ddp412/sux002/cis_heritability_data/1000G.EUR."
RSCRIPT = "Rscript"


def test_bed_reader_matches_pandas_plink():
    """Our direct .bed slice must equal what pandas_plink returns."""
    from pandas_plink import read_plink

    i0, i1 = 1000, 1500
    idx = build_chrom_index(GENO, 22, None)
    mine = read_bed_slice(GENO, 22, i0, i1, idx.n_samples)
    bim, _, G = read_plink(f"{GENO}22", verbose=False)
    ref = np.asarray(G[i0:i1, :].compute()).T
    assert np.allclose(np.nan_to_num(mine, nan=-9), np.nan_to_num(ref, nan=-9))

    rows = read_bim_rows(GENO, 22, i0, i1)
    assert (rows["snp"].to_numpy() == bim["snp"].to_numpy()[i0:i1]).all()
    assert (rows["a0"].to_numpy() == bim["a0"].to_numpy()[i0:i1]).all()


def test_marginal_ols_matches_statsmodels():
    """The vectorised sumstats must equal a per-SNP statsmodels OLS fit."""
    import statsmodels.api as sm

    rng = np.random.default_rng(0)
    n, p = 200, 40
    Z = rng.standard_normal((n, p))
    Z = (Z - Z.mean(0)) / Z.std(0)
    y = rng.standard_normal(n)
    y = (y - y.mean()) / y.std()

    beta, se, pval = marginal_ols(Z, y)
    for j in range(p):
        fit = sm.OLS(y, sm.add_constant(Z[:, j : j + 1])).fit()
        assert np.isclose(beta[j], fit.params[1])
        assert np.isclose(se[j], fit.bse[1])
        assert np.isclose(pval[j], fit.pvalues[1])


def test_beta_hits_target_genetic_variance():
    """After rescaling, beta' R beta must equal the target h2 exactly."""
    rng = np.random.default_rng(1)
    p = 200
    A = rng.standard_normal((400, p))
    A = (A - A.mean(0)) / A.std(0)
    ld = build_ld(A, ridge=0.1)
    for dist in ("normal", "fixed", "pareto"):
        for h2 in (0.05, 0.4, 0.8):
            beta, causal = sim_beta(ld.L, 5, h2, rng, effect_dist=dist)
            assert np.isclose(compute_s2g(ld.L, beta), h2, rtol=1e-8)
            assert causal.size == 5


def test_sim_trait_realises_target_h2():
    rng = np.random.default_rng(2)
    g = rng.standard_normal(200_000)
    y, h2_sample = sim_trait(g, 0.4, rng)
    assert abs(h2_sample - 0.4) < 0.01
    assert abs(y.mean()) < 1e-9 and abs(y.std() - 1) < 1e-9


def test_sim_geno_reproduces_ld():
    """Genotypes sampled through L must have the LD they were sampled from."""
    rng = np.random.default_rng(3)
    p = 60
    A = rng.standard_normal((300, p))
    A = (A - A.mean(0)) / A.std(0)
    ld = build_ld(A, ridge=0.1)
    Z = sim_geno(ld.L, 40_000, rng)
    R_hat = (Z.T @ Z) / Z.shape[0]
    off = ~np.eye(p, dtype=bool)
    assert np.abs(R_hat[off] - ld.R_sim[off]).max() < 0.05


def test_matrix_round_trip_to_r():
    """Cell vectors must land as R *columns* - the bug class this guards."""
    vectors = [np.arange(5.0) + 100 * k for k in range(3)]
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "m.f64")
        _write_cells(vectors, path)
        script = (
            f'v <- readBin("{path}", "double", 15, 8, endian="little"); '
            'm <- matrix(v, nrow=5, ncol=3); cat(m[1,], m[5,], sep=" ")'
        )
        out = subprocess.run([RSCRIPT, "--vanilla", "-e", script],
                             capture_output=True, text=True, check=True).stdout
    assert [float(x) for x in out.split()] == [0.0, 100.0, 200.0, 4.0, 104.0, 204.0]


def test_grm_binary_round_trip():
    """GCTA binary GRM layout: lower triangle, row-major, float32."""
    grm = np.array([[1.0, 0.2, 0.3], [0.2, 1.1, 0.4], [0.3, 0.4, 0.9]])
    with tempfile.TemporaryDirectory() as d:
        prefix = os.path.join(d, "g")
        write_grm_binary(grm, 500, prefix)
        vals = np.fromfile(f"{prefix}.grm.bin", dtype="<f4")
        assert np.allclose(vals, [1.0, 0.2, 1.1, 0.3, 0.4, 0.9], atol=1e-6)
        assert np.allclose(np.fromfile(f"{prefix}.grm.N.bin", dtype="<f4"), 500)
        assert sum(1 for _ in open(f"{prefix}.grm.id")) == 3


def _ld_and_sumstats(h2, ncausal, rng, n=400, p=250):
    """One simulated locus with banded LD; returns (R, marginal beta, true h2)."""
    from cish2.sumstats import marginal_ols

    A = rng.standard_normal((n, p))
    for j in range(1, p):
        A[:, j] = 0.6 * A[:, j - 1] + 0.8 * A[:, j]
    Z = (A - A.mean(0)) / A.std(0)
    b = np.zeros(p)
    b[rng.choice(p, ncausal, replace=False)] = rng.normal(0, 1, ncausal)
    g = Z @ b
    g *= np.sqrt(h2) / g.std()
    y = g + rng.normal(0, np.sqrt(1 - h2), n)
    y = (y - y.mean()) / y.std()
    R = (Z.T @ Z) / n
    np.fill_diagonal(R, 1.0)
    beta, _, _ = marginal_ols(Z, y)
    return R, beta, float(np.var(g, ddof=1) / np.var(y, ddof=1)), n


def test_heels_eigenbasis_matches_explicit_matrix_solve():
    """The eigenbasis iteration must equal an explicit O(p^3) solve per step."""
    rng = np.random.default_rng(7)
    R, beta, _, n = _ld_and_sumstats(0.4, 20, rng)
    em = HeelsSolver.from_ld(R).solve_em(beta, n, tol=1e-10, max_iter=20000)
    h2_ref, s2g_ref, s2e_ref = heels_reference(R, beta, n, tol=1e-10, max_iter=20000)
    assert np.isclose(em.h2, h2_ref, rtol=1e-6)
    assert np.isclose(em.s2g, s2g_ref, rtol=1e-6)
    assert np.isclose(em.s2e, s2e_ref, rtol=1e-6)


def test_heels_direct_solver_equals_converged_iteration():
    """Maximising the profile likelihood must land on the iteration's fixed
    point - the whole justification for not iterating."""
    rng = np.random.default_rng(17)
    for h2 in (0.05, 0.3, 0.7):
        R, beta, _, n = _ld_and_sumstats(h2, 10, rng)
        solver = HeelsSolver.from_ld(R)
        direct = solver.solve(beta, n)
        em = solver.solve_em(beta, n, tol=0.0, max_iter=200_000)
        assert abs(direct.h2 - em.h2) < 1e-3, f"h2={h2}: {direct.h2} vs {em.h2}"


def test_heels_from_genotypes_matches_from_ld():
    """The SVD path and the explicit-LD path are the same decomposition."""
    rng = np.random.default_rng(23)
    n, p = 400, 150          # n > p so the LD is full rank and both agree exactly
    A = rng.standard_normal((n, p))
    for j in range(1, p):
        A[:, j] = 0.6 * A[:, j - 1] + 0.8 * A[:, j]
    Z = (A - A.mean(0)) / A.std(0)
    y = rng.standard_normal(n)
    y = (y - y.mean()) / y.std()
    beta, _, _ = marginal_ols(Z, y)
    R = (Z.T @ Z) / n
    np.fill_diagonal(R, 1.0)
    a = HeelsSolver.from_genotypes(Z).solve(beta, n)
    b = HeelsSolver.from_ld(R).solve(beta, n)
    assert np.isclose(a.h2, b.h2, atol=1e-8)


def test_heels_recovers_known_heritability():
    """HEELS assumes the infinitesimal model, so it should be unbiased even
    when the truth is fully polygenic - the regime SuSiE cannot handle."""
    rng = np.random.default_rng(11)
    for h2 in (0.2, 0.5):
        for ncausal in (5, 250):
            est, truth = [], []
            for _ in range(30):
                R, beta, h2_sample, n = _ld_and_sumstats(h2, ncausal, rng)
                res = HeelsSolver.from_ld(R).solve(beta, n)
                if np.isfinite(res.h2):
                    est.append(res.h2)
                    truth.append(h2_sample)
            bias = float(np.mean(est) - np.mean(truth))
            assert abs(bias) < 0.05, f"h2={h2} ncausal={ncausal} bias={bias:+.3f}"


def test_heels_effective_dof_changes_the_answer():
    """Residualising spends degrees of freedom; ignoring that inflates h2.

    Guards the GTEx bug where passing n instead of n - n_covar pinned every
    gene at the h2 = 1 bound.
    """
    rng = np.random.default_rng(29)
    n, p, k = 120, 300, 30
    A = rng.standard_normal((n, p))
    for j in range(1, p):
        A[:, j] = 0.6 * A[:, j - 1] + 0.8 * A[:, j]
    Z = (A - A.mean(0)) / A.std(0)
    C = np.column_stack([np.ones(n), rng.standard_normal((n, k))])
    Q, _ = np.linalg.qr(C)
    y = rng.standard_normal(n)                       # no genetic signal at all
    y = y - Q @ (Q.T @ y)
    y = (y - y.mean()) / y.std()
    Zr = Z - Q @ (Q.T @ Z)
    Zr = (Zr - Zr.mean(0)) / Zr.std(0)
    beta, _, _ = marginal_ols(Zr, y)
    solver = HeelsSolver.from_genotypes(Zr)
    naive = solver.solve(beta, n).h2
    corrected = solver.solve(beta, n, n_eff=n - Q.shape[1]).h2
    assert corrected < naive
    assert corrected < 0.5, f"null phenotype should not read h2={corrected:.3f}"


def test_heels_se_is_positive_and_finite():
    rng = np.random.default_rng(13)
    R, beta, _, n = _ld_and_sumstats(0.5, 30, rng)
    res = HeelsSolver.from_ld(R).solve(beta, n)
    assert np.isfinite(res.h2_se) and res.h2_se > 0
    assert 0.0 <= res.h2 <= 1.0


def test_parse_hsq_reads_by_label():
    body = (
        "Source\tVariance\tSE\n"
        "V(G)\t0.31\t0.09\nV(e)\t0.68\t0.09\nVp\t0.99\t0.08\n"
        "V(G)/Vp\t0.313\t0.087\nlogL\t-120.1\nlogL0\t-125.0\n"
        "LRT\t9.8\ndf\t1\nPval\t8.7e-4\nn\t300\n"
    )
    with tempfile.NamedTemporaryFile("w", suffix=".hsq", delete=False) as fh:
        fh.write(body)
        path = fh.name
    res = parse_hsq(path)
    os.unlink(path)
    assert np.isclose(res.h2, 0.313) and np.isclose(res.h2_se, 0.087)
    assert np.isclose(res.vg, 0.31) and np.isclose(res.pval, 8.7e-4)
    assert res.converged


if __name__ == "__main__":
    failures = 0
    for name, fn in sorted(globals().items()):
        if not name.startswith("test_") or not callable(fn):
            continue
        try:
            fn()
            print(f"PASS  {name}")
        except Exception as exc:
            failures += 1
            print(f"FAIL  {name}: {type(exc).__name__}: {exc}")
    print(f"\n{failures} failure(s)")
    sys.exit(1 if failures else 0)


def test_susie_estimand_matches_brute_force_second_moment():
    """h2_full in R/susie_batch.R must equal tr(R E[bb']) exactly.

    The trace form and the full quadratic form differ by the within-single-effect
    term sum_l m_l'(R - I) m_l, which is not negligible when a credible set
    spreads over SNPs in LD.  Checked against an explicitly built second-moment
    matrix, not against the algebra it was derived from.
    """
    script = os.path.join(os.path.dirname(os.path.abspath(__file__)), "test_susie_estimand.R")
    proc = subprocess.run(["Rscript", "--vanilla", script],
                                                    capture_output=True, text=True, timeout=900)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert "ok:" in proc.stdout, proc.stdout


# ---------------------------------------------------------------------------
# Round 4 (PI feedback): semi-synthetic / permutation modes for the GTEx runner
# ---------------------------------------------------------------------------
GTEX_SPLIT_HALF = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                               "configs", "gtex_split_half.yaml")


def _sim_beta_round3(L, ncausal, h2, rng, effect_dist="normal"):
    """Verbatim copy of sim_beta's draw before the round-4 refactor."""
    p = L.shape[0]
    beta = np.zeros(p)
    causal = rng.choice(np.arange(p), size=min(ncausal, p), replace=False)
    k = causal.size
    if effect_dist == "normal":
        effects = rng.normal(loc=0.0, scale=np.sqrt(h2 / k), size=k)
    elif effect_dist == "fixed":
        effects = rng.choice([-1.0, 1.0], size=k) * np.sqrt(h2 / k)
    else:
        effects = rng.choice([-1.0, 1.0], size=k) * rng.pareto(a=2.0, size=k)
    np.add.at(beta, causal, effects)
    beta *= np.sqrt(h2 / compute_s2g(L, beta))
    return beta, np.unique(causal)


def test_sim_beta_refactor_is_bit_identical():
    """Sharing _draw_effects must not move a single simulated effect."""
    rng0 = np.random.default_rng(11)
    A = rng0.standard_normal((300, 80))
    ld = build_ld((A - A.mean(0)) / A.std(0), ridge=0.1)
    for dist in ("normal", "fixed", "pareto"):
        for seed in range(5):
            new = sim_beta(ld.L, 5, 0.3, np.random.default_rng(seed), effect_dist=dist)
            old = _sim_beta_round3(ld.L, 5, 0.3, np.random.default_rng(seed), effect_dist=dist)
            assert np.array_equal(new[0], old[0]) and np.array_equal(new[1], old[1])


def test_sim_beta_empirical_hits_target_variance():
    from cish2.simulate import sim_beta_empirical

    rng = np.random.default_rng(12)
    G = rng.binomial(2, 0.3, size=(574, 300)).astype(float)
    Z = (G - G.mean(0)) / G.std(0)
    for k in (1, 5):
        beta, causal = sim_beta_empirical(Z, k, 0.3, rng)
        assert abs(np.var(Z @ beta) - 0.3) < 1e-10
        assert np.count_nonzero(beta) == k and causal.size == k
    beta, causal = sim_beta_empirical(Z, 5, 0.0, rng)
    assert not beta.any() and causal.size == 0


class _FakeGtex:
    """Just enough of GtexDataset for GeneWorker._run_chunk, held in memory."""

    def __init__(self, n=160, p=70, n_genes=3, seed=0):
        import pandas as pd

        rng = np.random.default_rng(seed)
        maf = rng.uniform(0.05, 0.5, p)
        self.G = rng.binomial(2, maf, size=(n, p)).astype(float)
        self.G[rng.integers(0, n, 5), rng.integers(0, p, 5)] = np.nan
        self.expr = rng.standard_normal((n_genes, n))
        self.genes = pd.DataFrame({"gene_id": [f"ENSG{i:011d}.1" for i in range(n_genes)],
                                   "chrom": 1, "tss": 1_000_000, "expr_row": range(n_genes)})
        self.covar = rng.standard_normal((n, 4))
        self.covar_names = ["PC1", "PC2", "sex", "InferredCov1"]
        self.n_samples = n
        self.p = p

    def window(self, row, window_bp):
        from cish2.gtex import GeneWindow

        return GeneWindow(gene_id=str(row["gene_id"]), chrom=1, tss=1_000_000,
                          start_bp=0, end_bp=1, i0=0, i1=self.p, expr_row=int(row["expr_row"]))

    def genotypes(self, win):
        return self.G.copy()


def _fake_worker(mode, tmpdir, **cfg_extra):
    """A GeneWorker whose GCTA/SuSiE calls are stubbed; the prepared cells are
    captured so tests can inspect exactly what each estimator would receive."""
    from cish2.estimators.gcta import GctaResult
    from cish2.gtex import disjoint_subsamples
    from cish2.gtex_campaign import DEFAULTS
    from cish2.gtex_runner import GeneWorker

    data = _FakeGtex()
    cfg = dict(DEFAULTS, min_snps=10, phenotype_mode=mode, **cfg_extra)
    cfg["subsets"] = disjoint_subsamples(data.n_samples, 2, cfg["seed"])
    w = GeneWorker(cfg, data, tmpdir)
    captured = []
    w._run_gcta = lambda cell, workdir, runner=None: GctaResult.failed("stub")
    w._run_susie = lambda cells, workdir: captured.extend(cells) or {}
    return w, data, captured


def test_simulated_mode_uses_one_beta_for_both_halves():
    import pandas as pd

    with tempfile.TemporaryDirectory() as d:
        w, data, cells = _fake_worker("simulated", d)
        rows = pd.DataFrame(w._run_chunk(0, list(range(len(data.genes))), d))
    assert len(rows) == len(data.genes) * 2 * 2          # genes x draws x halves
    per = rows.groupby(["gene_id", "draw"])
    assert (per["subset"].nunique() == 2).all()
    for col in ("sim_beta_hash", "h2_target", "ncausal_target"):
        assert (per[col].nunique() == 1).all(), col       # same truth in both halves
    assert len({c["cell_id"] for c in cells}) == len(cells)   # draws stay distinct
    # deterministic: a second worker draws the same phenotypes
    with tempfile.TemporaryDirectory() as d:
        w2, _, _ = _fake_worker("simulated", d)
        rows2 = pd.DataFrame(w2._run_chunk(7, list(range(len(data.genes))), d))
    assert (rows["sim_beta_hash"].to_numpy() == rows2["sim_beta_hash"].to_numpy()).all()


def test_permuted_mode_changes_only_genotype_rows():
    with tempfile.TemporaryDirectory() as d:
        w, data, real = _fake_worker("real", d)
        w._run_chunk(0, [0], d)
        wp, _, perm = _fake_worker("permuted", d)
        wp._run_chunk(0, [0], d)
    assert len(real) == len(perm) == 2
    for a, b in zip(real, perm):
        assert a["subset"] == b["subset"]
        # the phenotype REML and SuSiE see is bit-identical to real mode
        assert np.array_equal(a["y_gcta"], b["y_gcta"])
        assert np.array_equal(a["y"], b["y"])
        # the genotypes are the same donors' rows, reordered
        assert not np.allclose(a["Z_gcta"], b["Z_gcta"])
        assert np.allclose(np.sort(a["Z_gcta"], axis=0), np.sort(b["Z_gcta"], axis=0), atol=1e-10)


def test_real_mode_reproduces_existing_split_half_rows():
    """Regression: phenotype_mode=real on one chunk == results/gtex_split_half."""
    import pandas as pd
    import pytest

    from cish2.gtex import disjoint_subsamples, load_dataset
    from cish2.gtex_campaign import load_config
    from cish2.gtex_runner import GeneWorker

    cfg = load_config(GTEX_SPLIT_HALF)
    old_path = os.path.join(cfg["results_root"], "gtex_split_half", "cells.tsv.gz")
    if not (os.path.exists(cfg["geno_prefix"] + ".bed") and os.path.exists(old_path)):
        pytest.skip("GTEx data or archived split-half results not available")
    data = load_dataset(cfg["geno_prefix"], cfg["expression"], cfg["covariates"],
                        cfg["chromosomes"])
    cfg["subsets"] = disjoint_subsamples(data.n_samples, int(cfg["n_splits"]), cfg["seed"])
    gene_rows = list(range(int(cfg["genes_per_chunk"])))           # chunk 0
    with tempfile.TemporaryDirectory() as d:
        out = GeneWorker(cfg, data, d).run_chunk(0, gene_rows)
    assert out.status == "ok", out.message
    new = pd.DataFrame(out.rows).set_index(["gene_id", "subset"]).sort_index()
    old = pd.read_csv(old_path, sep="\t").set_index(["gene_id", "subset"])
    old = old.loc[old.index.intersection(new.index)].sort_index()
    assert len(old) == len(new) > 0
    for col in ("p_snps", "n", "n_clumps", "gcta_h2", "gcta_h2_se", "susie_h2",
                "susie_h2_se", "heels_h2", "max_chisq", "min_pval", "maf_mean"):
        a, b = new[col].to_numpy(float), old[col].to_numpy(float)
        both = np.isfinite(a) & np.isfinite(b)
        assert (np.isfinite(a) == np.isfinite(b)).all(), col
        assert np.allclose(a[both], b[both], rtol=0, atol=1e-8), \
            (col, float(np.max(np.abs(a[both] - b[both]))))
