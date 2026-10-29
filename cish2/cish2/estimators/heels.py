"""HEELS — heritability from summary statistics and the LD matrix.

Li, Mazumder & Lin, *Nature Communications* 14:7954 (2023),
"Accurate and efficient estimation of local heritability using summary
statistics and the linkage disequilibrium matrix".

Unlike SuSiE, HEELS makes no sparsity assumption: it fits the same
infinitesimal linear mixed model REML fits, but from summary statistics plus
LD rather than individual-level genotypes.  It is therefore the natural third
arm here — a summary-statistic method that should behave like GCTA rather than
like fine-mapping.

Model (following the paper's scaling, which the reference implementation
`huilisabrina/HEELS` confirms)::

    y = X~ b + e,   b ~ N(0, s2g I_p),   e ~ N(0, s2e I_n)

with the genotype matrix column-standardised and then divided by sqrt(p),
``X~ = X_std / sqrt(p)``, and ``y`` standardised so ``y'y / n = 1``.  That
scaling is what makes ``s2g`` the *total* local genetic variance rather than a
per-SNP variance, so ``h2 = s2g / (s2g + s2e)`` directly.

The sufficient statistics are then

    S  = X~' y      = n * beta_marginal / sqrt(p)
    R~ = X~' X~     = (n / p) * R          (R = the LD correlation matrix)

and the algorithm alternates

    W       = (s2e / s2g) I + R~
    b_hat   = W^-1 S                              (BLUP)
    s2g     <- (b_hat' b_hat + s2e * tr(W^-1)) / p
    s2e     <- y'y / n  -  S' b_hat / n

until ``s2g`` stops moving.

**In-sample LD is a hard requirement, not a preference.**  Non-negativity of
the residual variance rests on ``S' R^+ S <= y'y``, which holds because
``S' R^+ S = y' P_X y`` is a projection of y when R is the LD of the *same*
individuals the summary statistics came from.  Hand the estimator a
reference-panel LD instead and that bound is simply false: on a 1000G region
with p = 893, n = 300 and a 503-individual panel we measured
``S' R^+ S = 1383`` against ``y'y = 300``, driving ``s2e`` negative and the
iteration to divergence within a few steps.  So HEELS is always given the LD of
its own sample here, and :meth:`HeelsSolver.from_genotypes` is the entry point.

**Why this implementation is fast.**  Every iteration needs ``W^-1``, and W
changes each iteration, so a naive version pays O(p^3) per iteration per grid
cell.  But W is a shift of R~ in a fixed eigenbasis, and that basis comes free
from the SVD of the n x p genotype matrix at O(n^2 p) rather than O(p^3) --
a few hundred milliseconds at n = 300, p = 8,000 instead of ~84 s.  Every
iteration is then O(p) and every new phenotype one O(np) projection.

When p > n the LD is singular, with ``p - n`` zero eigenvalues.  Those
directions are handled exactly rather than dropped: S has no component in that
null space (it lives in the row space of X), so they contribute nothing to the
BLUP, but they *do* contribute ``(p - rank) / a`` to ``tr(W^-1)`` and therefore
to the genetic-variance update.

The paper's low-rank-plus-banded LD approximation exists to make genome-wide
scale tractable; at one cis-region we can afford the exact decomposition, so
this is the exact estimator rather than its scalable approximation.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np

DEFAULT_TOL = 1e-5
DEFAULT_MAX_ITER = 200
_MIN_VAR = 1e-8


def _robust_svd(Z: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Right singular basis of ``Z``, with fallbacks.

    NumPy's default driver (LAPACK ``gesdd``) is fast but occasionally fails to
    converge on real genotype matrices -- it did so on roughly 12% of GTEx gene
    chunks.  ``gesvd`` is slower and more robust; eigendecomposing ``Z'Z`` is
    the last resort and always succeeds on a PSD matrix.
    """
    try:
        _, sv, Vt = np.linalg.svd(Z, full_matrices=False)
        return sv, Vt
    except np.linalg.LinAlgError:
        pass
    try:
        from scipy.linalg import svd as scipy_svd

        _, sv, Vt = scipy_svd(Z, full_matrices=False, lapack_driver="gesvd")
        return sv, Vt
    except Exception:
        pass
    # Gram-matrix route: eigenvectors of Z'Z are the right singular vectors
    n, p = Z.shape
    gram = Z.T @ Z if p <= n else None
    if gram is not None:
        w, V = np.linalg.eigh(gram)
        order = np.argsort(w)[::-1]
        return np.sqrt(np.maximum(w[order], 0.0)), V[:, order].T
    w, U = np.linalg.eigh(Z @ Z.T)
    order = np.argsort(w)[::-1]
    w = np.maximum(w[order], 0.0)
    sv = np.sqrt(w)
    keep = sv > 0
    Vt = (Z.T @ U[:, order][:, keep] / sv[keep]).T
    return sv[keep], Vt


@dataclass
class HeelsResult:
    """One HEELS fit."""

    h2: float = np.nan
    h2_se: float = np.nan
    s2g: float = np.nan
    s2e: float = np.nan
    niter: float = np.nan
    loglik: float = np.nan
    converged: bool = False
    status: str = "not_run"

    @classmethod
    def failed(cls, status: str) -> "HeelsResult":
        return cls(status=status)


class HeelsSolver:
    """Fits HEELS for many phenotypes against one region's LD.

    :param eigvals: eigenvalues of the p x p LD **correlation** matrix.
    :param eigvecs: the matching eigenvectors (columns).

    Build once per region with :meth:`from_ld`, then call :meth:`solve` per
    phenotype.
    """

    def __init__(self, eigvals: np.ndarray, eigvecs: np.ndarray, p: Optional[int] = None):
        self.d = np.asarray(eigvals, dtype=np.float64)   # non-zero LD eigenvalues
        self.U = np.asarray(eigvecs, dtype=np.float64)   # their eigenvectors, p x r
        self.p = int(p if p is not None else self.d.size)
        self.rank = int(self.d.size)

    @classmethod
    def from_genotypes(cls, Z: np.ndarray, tol: float = 1e-10) -> "HeelsSolver":
        """Build from the standardised genotypes themselves - the fast path.

        The in-sample LD is ``Z'Z / n``, so its eigenbasis is the right singular
        basis of ``Z``.  Taking it from an SVD of the n x p matrix costs
        O(n^2 p) instead of the O(p^3) a full eigendecomposition of the p x p
        LD would need, and it is the *same* decomposition.
        """
        Z = np.asarray(Z, dtype=np.float64)
        n, p = Z.shape
        # economy SVD: Z = U_z diag(s) V',  so Z'Z/n = V diag(s^2/n) V'
        sv, Vt = _robust_svd(Z)
        d = (sv**2) / n
        keep = d > tol * max(float(d[0]) if d.size else 0.0, 1.0)
        return cls(d[keep], Vt[keep].T, p=p)

    @classmethod
    def from_ld(cls, R: np.ndarray) -> "HeelsSolver":
        """Build from an explicit LD correlation matrix (used by the tests)."""
        d, U = np.linalg.eigh(np.asarray(R, dtype=np.float64))
        d = np.maximum(d, 0.0)  # tiny negatives are float noise on a PSD matrix
        return cls(d, U, p=R.shape[0])

    # ------------------------------------------------------------------ fit --
    def _prepare(self, beta: np.ndarray, n: int):
        """Project the summary statistics into the LD eigenbasis."""
        p = self.p
        S = n * np.asarray(beta, dtype=np.float64) / np.sqrt(p)
        lam = (n / p) * self.d
        c2 = (self.U.T @ S) ** 2
        return c2, lam, p - self.rank

    def _profile_ll(self, a: float, c2: np.ndarray, lam: np.ndarray,
                    n_eff: float, yty: float) -> float:
        """Log-likelihood with the variance scale profiled out."""
        if not np.isfinite(a) or a <= 0:
            return -np.inf
        Q = yty - float(np.sum(c2 / (a + lam)))
        if Q <= 0:
            return -np.inf
        return -0.5 * (n_eff * np.log(Q / n_eff) + float(np.sum(np.log1p(lam / a))))

    def solve(
        self,
        beta: np.ndarray,
        n: int,
        yty_over_n: float = 1.0,
        tol: float = DEFAULT_TOL,
        max_iter: int = DEFAULT_MAX_ITER,
        init_h2: float = 0.25,
        h2_bounds: Tuple[float, float] = (1e-6, 0.9999),
        grid_points: int = 256,
        n_eff: Optional[float] = None,
    ) -> HeelsResult:
        """Estimate local heritability from marginal effects.

        :param beta: marginal (per-SNP simple regression) effect estimates, on
            the standardised-genotype / standardised-phenotype scale.
        :param n: number of individuals behind those estimates.  Sets the
            scale of the sufficient statistics.
        :param n_eff: **degrees of freedom** left in the phenotype, which is
            ``n`` only when no covariates were removed.  Residualising on k
            covariates leaves ``n - k``, and using ``n`` instead inflates the
            estimate badly -- on GTEx at n = 100 with 24 covariates it pins
            every gene at the h2 = 1 bound.  Defaults to ``n``.
        :param yty_over_n: sample variance of the phenotype; 1.0 when ``y`` is
            standardised, which is how this pipeline produces it.

        ``tol``, ``max_iter`` and ``init_h2`` are accepted for interface
        compatibility with the iterative form and are unused here.
        """
        beta = np.asarray(beta, dtype=np.float64)
        if beta.size != self.p:
            return HeelsResult.failed(f"beta length {beta.size} != p {self.p}")
        if not np.all(np.isfinite(beta)):
            return HeelsResult.failed("non-finite marginal effects")

        c2, lam, n_null = self._prepare(beta, n)
        yty = float(n) * yty_over_n
        dof = float(n_eff) if n_eff is not None else float(n)
        if dof <= 1:
            return HeelsResult.failed(f"degrees of freedom exhausted (n_eff={dof:g})")

        # a = s2e/s2g, so the h2 bounds map straight onto a search interval
        lo_h2, hi_h2 = h2_bounds
        log_lo, log_hi = np.log((1 - hi_h2) / hi_h2), np.log((1 - lo_h2) / lo_h2)

        # coarse sweep first: the profile flattens near the boundary, and a
        # bracket taken from a grid is far more robust than one guessed from a
        # starting value
        grid = np.linspace(log_lo, log_hi, grid_points)
        ll = np.array([self._profile_ll(float(np.exp(t)), c2, lam, dof, yty)
                       for t in grid])
        if not np.any(np.isfinite(ll)):
            return HeelsResult.failed("profile likelihood undefined everywhere")
        k = int(np.argmax(np.where(np.isfinite(ll), ll, -np.inf)))
        left = grid[max(k - 1, 0)]
        right = grid[min(k + 1, grid_points - 1)]

        from scipy.optimize import minimize_scalar

        opt = minimize_scalar(
            lambda t: -self._profile_ll(float(np.exp(t)), c2, lam, dof, yty),
            bounds=(left, right), method="bounded", options={"xatol": 1e-10},
        )
        t_star = float(opt.x) if opt.success else float(grid[k])
        a = float(np.exp(t_star))

        Q = yty - float(np.sum(c2 / (a + lam)))
        s2e = max(Q / dof, _MIN_VAR)
        s2g = max(s2e / a, _MIN_VAR)
        h2 = 1.0 / (1.0 + a)

        at_bound = bool(t_star <= log_lo + 1e-9 or t_star >= log_hi - 1e-9)
        loglik = self._loglik(s2g, s2e, c2, lam, dof, yty / dof, n_null)
        se = self._h2_se(s2g, s2e, c2, lam, dof, yty / dof, n_null)
        return HeelsResult(
            h2=float(h2), h2_se=se, s2g=float(s2g), s2e=float(s2e),
            niter=float(opt.nfev if opt.success else grid_points),
            loglik=loglik, converged=not at_bound,
            status="ok" if not at_bound else "at_h2_bound",
        )

    def solve_em(
        self,
        beta: np.ndarray,
        n: int,
        yty_over_n: float = 1.0,
        tol: float = DEFAULT_TOL,
        max_iter: int = DEFAULT_MAX_ITER,
        init_h2: float = 0.25,
        n_eff: Optional[float] = None,
    ) -> HeelsResult:
        """The paper's alternating iteration, kept for validation.

        Converges to the same point as :meth:`solve`, but can need tens of
        thousands of steps near the h2 -> 0 boundary.
        """
        beta = np.asarray(beta, dtype=np.float64)
        if beta.size != self.p:
            return HeelsResult.failed(f"beta length {beta.size} != p {self.p}")
        if not np.all(np.isfinite(beta)):
            return HeelsResult.failed("non-finite marginal effects")

        c2, lam, n_null = self._prepare(beta, n)
        p = self.p
        yty = float(n) * yty_over_n
        dof = float(n_eff) if n_eff is not None else float(n)
        s2g = max(init_h2 * yty_over_n, _MIN_VAR)
        s2e = max((1.0 - init_h2) * yty_over_n, _MIN_VAR)

        converged = False
        it = 0
        for it in range(1, max_iter + 1):
            a = s2e / s2g
            denom = a + lam
            if not np.all(np.isfinite(denom)) or np.any(denom <= 0):
                return HeelsResult.failed("degenerate W in iteration")
            inv = 1.0 / denom
            bhat_sq = float(np.sum(c2 * inv**2))
            trace_w = float(np.sum(inv)) + n_null / a
            s_bhat = float(np.sum(c2 * inv))

            s2g_new = max((bhat_sq + s2e * trace_w) / p, _MIN_VAR)
            s2e_new = max((yty - s_bhat) / dof, _MIN_VAR)
            delta = abs(s2g_new - s2g)
            s2g, s2e = s2g_new, s2e_new
            if delta < tol:
                converged = True
                break

        h2 = s2g / (s2g + s2e)
        return HeelsResult(
            h2=float(h2), h2_se=np.nan, s2g=float(s2g), s2e=float(s2e),
            niter=float(it),
            loglik=self._loglik(s2g, s2e, c2, lam, dof, yty / dof, n_null),
            converged=converged, status="ok" if converged else "not_converged",
        )

    # ------------------------------------------------------- log-likelihood --
    @staticmethod
    def _loglik(s2g: float, s2e: float, c2: np.ndarray, lam: np.ndarray,
                n: float, yty_over_n: float, n_null: int = 0) -> float:
        """Marginal log-likelihood, up to a constant.

        ``y ~ N(0, s2g X~X~' + s2e I)``; Woodbury and the eigenbasis reduce both
        the determinant and the quadratic form to sums over eigenvalues.
        """
        if s2g <= 0 or s2e <= 0:
            return float("nan")
        a = s2e / s2g
        # null-space directions contribute log(1 + 0) = 0 to the determinant
        logdet = n * np.log(s2e) + float(np.sum(np.log1p((s2g / s2e) * lam)))
        quad = (n * yty_over_n - float(np.sum(c2 / (a + lam)))) / s2e
        return float(-0.5 * (logdet + quad))

    def _h2_se(self, s2g: float, s2e: float, c2: np.ndarray, lam: np.ndarray,
               n: int, yty_over_n: float, n_null: int = 0) -> float:
        """Delta-method SE from the numerically differentiated Hessian.

        The paper derives the Fisher information analytically; here the
        log-likelihood is O(p) to evaluate, so a central-difference Hessian is
        both cheap and free of transcription risk.
        """
        theta = np.array([s2g, s2e], dtype=np.float64)
        step = np.maximum(np.abs(theta) * 1e-4, 1e-7)

        def ll(t: np.ndarray) -> float:
            if np.any(t <= 0):
                return float("nan")
            return self._loglik(t[0], t[1], c2, lam, n, yty_over_n, n_null)

        H = np.zeros((2, 2))
        for i in range(2):
            for j in range(i, 2):
                tpp, tpm, tmp, tmm = (theta.copy() for _ in range(4))
                tpp[i] += step[i]; tpp[j] += step[j]
                tpm[i] += step[i]; tpm[j] -= step[j]
                tmp[i] -= step[i]; tmp[j] += step[j]
                tmm[i] -= step[i]; tmm[j] -= step[j]
                val = (ll(tpp) - ll(tpm) - ll(tmp) + ll(tmm)) / (4 * step[i] * step[j])
                H[i, j] = H[j, i] = val

        try:
            cov = np.linalg.inv(-H)
        except np.linalg.LinAlgError:
            return float("nan")
        if not np.all(np.isfinite(cov)):
            return float("nan")

        total = s2g + s2e
        grad = np.array([s2e / total**2, -s2g / total**2])
        var = float(grad @ cov @ grad)
        return float(np.sqrt(var)) if var > 0 else float("nan")


def heels_reference(R: np.ndarray, beta: np.ndarray, n: int,
                    yty_over_n: float = 1.0, tol: float = DEFAULT_TOL,
                    max_iter: int = DEFAULT_MAX_ITER,
                    init_h2: float = 0.25) -> Tuple[float, float, float]:
    """Direct O(p^3)-per-iteration HEELS, for validating the fast solver.

    Solves ``W b = S`` with an explicit factorisation each iteration instead of
    reusing an eigenbasis.  Only used by the tests.
    """
    p = R.shape[0]
    S = n * np.asarray(beta, dtype=np.float64) / np.sqrt(p)
    R_t = (n / p) * np.asarray(R, dtype=np.float64)

    s2g = max(init_h2 * yty_over_n, _MIN_VAR)
    s2e = max((1.0 - init_h2) * yty_over_n, _MIN_VAR)
    for _ in range(max_iter):
        W = R_t + np.eye(p) * (s2e / s2g)
        W_inv = np.linalg.inv(W)
        bhat = W_inv @ S
        s2g_new = max((bhat @ bhat + s2e * np.trace(W_inv)) / p, _MIN_VAR)
        s2e_new = max(yty_over_n - (S @ bhat) / n, _MIN_VAR)
        delta = abs(s2g_new - s2g)
        s2g, s2e = s2g_new, s2e_new
        if delta < tol:
            break
    return s2g / (s2g + s2e), s2g, s2e
