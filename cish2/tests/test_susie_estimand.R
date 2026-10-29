# Does the closed-form h2 in R/susie_batch.R equal tr(R E[bb']) exactly?
#
# The trace form sum(mu^2 + sd^2) and the full quadratic form differ by the
# within-single-effect term sum_l m_l' (R - I) m_l, which is large whenever a
# credible set spreads over SNPs in LD.  This checks the formula against a
# brute-force second-moment matrix rather than against the algebra it came from.
suppressMessages(library(susieR))
set.seed(11)
fails <- 0
for (rep in 1:5) {
  n <- 400; p <- 60
  L0 <- matrix(rnorm(p * p), p, p)
  S <- cov2cor(crossprod(L0) + 3 * diag(p))
  X <- scale(matrix(rnorm(n * p), n, p) %*% chol(S))
  b <- rep(0, p); b[sample(p, 3)] <- c(0.4, 0.3, 0.5)
  y <- scale(X %*% b + rnorm(n))
  R <- cor(X)
  fit <- susie(X, y, L = 10)

  idx <- if (is.numeric(fit$V)) which(fit$V > 1e-9) else seq_len(nrow(fit$alpha))
  A <- fit$alpha[idx, , drop = FALSE]
  M <- fit$mu[idx, , drop = FALSE]
  M2 <- fit$mu2[idx, , drop = FALSE]
  Wl <- A * M
  mu <- colSums(Wl)
  within <- sum(A * M2)
  self_R <- sum(t(Wl) * (R %*% t(Wl)))
  self_I <- sum(Wl * Wl)
  h2_full <- as.numeric(crossprod(mu, R %*% mu)) - self_R + within
  sd <- sqrt(pmax(colSums(A * M2 - Wl^2), 0))

  Ebb <- matrix(0, p, p)
  for (l in seq_along(idx)) {
    Ebb <- Ebb + diag(A[l, ] * M2[l, ])
    for (m in seq_along(idx)) if (m != l) Ebb <- Ebb + outer(Wl[l, ], Wl[m, ])
  }
  d_full  <- abs(sum(R * t(Ebb)) - h2_full)
  d_trace <- abs((sum(mu^2) - self_I + within) - (sum(mu^2) + sum(sd^2)))
  d_mu    <- max(abs(mu - coef(fit)[-1]))
  d_sd    <- max(abs(sd - susie_get_posterior_sd(fit)))
  if (d_full > 1e-9 || d_trace > 1e-9 || d_mu > 1e-9 || d_sd > 1e-9) {
    fails <- fails + 1
    cat(sprintf("rep %d FAIL full=%.2e trace=%.2e mu=%.2e sd=%.2e\n",
                rep, d_full, d_trace, d_mu, d_sd))
  }
}
if (fails > 0) stop(sprintf("%d/5 reps failed", fails))
cat("ok: h2_full == tr(R E[bb']) and collapses to the trace form at R = I\n")
