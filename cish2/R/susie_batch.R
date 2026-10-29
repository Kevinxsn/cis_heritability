#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Batched summary-statistic fine-mapping for one cis-region.
#
# The legacy run_susie.R was invoked once per grid cell and each invocation
# paid for: R startup, library loading, parsing a multi-GB PLINK text LD file,
# and rebuilding a p x p matrix from long-format pairs.  For a 30-cell grid on
# one region that is 30x the same work on the same matrix.
#
# This script is invoked once per (region, sample size).  It reads the LD
# matrix once, as raw float64 written by Python, then loops over every cell of
# the grid.  Per-SNP posteriors stay in memory; only the h2 summary each cell
# contributes is written out.
#
# Inputs (all inside --dir):
#   ld.f64        p*p float64, column-major-symmetric (symmetry makes the
#                 row/column-major distinction irrelevant)
#   bhat.f64      p*K float64, column k = marginal betas of cell k
#   shat.f64      p*K float64, column k = their standard errors
#   cells.tsv     cell_id, n  (one row per cell, in column order)
#   params.tsv    key, value
# Output:
#   susie.tsv     one row per cell
# ---------------------------------------------------------------------------

suppressMessages(library("susieR"))
suppressMessages(library("data.table"))

args <- commandArgs(trailingOnly = TRUE)
dir <- NULL
for (i in seq_along(args)) if (args[i] == "--dir") dir <- args[i + 1L]
if (is.null(dir)) stop("--dir <region scratch dir> is required")

read_f64 <- function(path, nrow, ncol) {
  v <- readBin(path, what = "double", n = as.double(nrow) * ncol, size = 8, endian = "little")
  if (length(v) != as.double(nrow) * ncol) stop(sprintf("short read: %s", path))
  matrix(v, nrow = nrow, ncol = ncol)
}

params <- fread(file.path(dir, "params.tsv"), header = TRUE, sep = "\t")
pget <- function(key, default = NA) {
  hit <- params$value[params$key == key]
  if (length(hit) == 0) default else hit[1]
}

p <- as.integer(pget("p"))
L <- as.integer(pget("susie_L", 10))
max_iter <- as.integer(pget("susie_max_iter", 100))
coverage <- as.numeric(pget("susie_coverage", 0.95))

cells <- fread(file.path(dir, "cells.tsv"), header = TRUE, sep = "\t")
K <- nrow(cells)

R <- read_f64(file.path(dir, "ld.f64"), p, p)
bhat <- read_f64(file.path(dir, "bhat.f64"), p, K)
shat <- read_f64(file.path(dir, "shat.f64"), p, K)

out <- data.table(
  cell_id = cells$cell_id,
  susie_h2 = NA_real_, susie_h2_se = NA_real_,
  susie_h2_quad = NA_real_,
  susie_sum_mu2 = NA_real_, susie_sum_sd2 = NA_real_,
  susie_n_cs = NA_integer_, susie_cs_size = NA_integer_,
  susie_max_pip = NA_real_, susie_sum_pip = NA_real_,
  susie_niter = NA_integer_, susie_sigma2 = NA_real_,
  susie_converged = FALSE, susie_status = "not_run"
)

for (k in seq_len(K)) {
  res <- tryCatch(
    withCallingHandlers(
      # check_prior stays at its default TRUE: with a rank-deficient panel LD
      # (p > n_reference) the IBSS fit can run away to an absurd prior
      # variance, and this guard is what turns that into a clean error the
      # driver can retry -- exactly how the legacy run_susie.R behaved.
      susie_rss(
        bhat = bhat[, k], shat = shat[, k], R = R, n = cells$n[k],
        L = L, max_iter = max_iter, coverage = coverage
      ),
      warning = function(w) invokeRestart("muffleWarning")
    ),
    error = function(e) e
  )

  if (inherits(res, "error")) {
    set(out, k, "susie_status", paste0("error: ", conditionMessage(res)))
    next
  }

  # --- two estimands, both from the same posterior ------------------------
  # SuSiE's variational posterior is independent ACROSS single effects, and
  # within one single effect exactly one SNP is non-zero.  Writing
  # m_l = alpha_l * mu_l and s_l = sum_j alpha_lj * mu2_lj,
  #
  #   E[b]      = sum_l m_l                                    = mu
  #   E[b b']   = sum_l diag(alpha_l * mu2_l)  +  sum_{l != l'} m_l m_l''
  #
  # so the genetic variance the simulation defines as truth is
  #
  #   h2_full = tr(R E[b b'])
  #           = mu' R mu  -  sum_l m_l' R m_l  +  sum_l s_l
  #
  # The middle term is the WITHIN-single-effect correction.  It is not
  # optional: a credible set puts its mass on SNPs in near-perfect LD, so
  # m_l' R m_l is large and positive, and dropping it inflates h2.  At R = I
  # the whole expression collapses to sum(mu^2 + sd^2), the trace form the
  # pipeline has always reported -- which is asserted below on every fit.
  idx <- if (is.numeric(res$V)) which(res$V > 1e-9) else seq_len(nrow(res$alpha))
  if (length(idx) == 0L) {
    mu <- rep(0, p); sum_mu2 <- 0; sum_sd2 <- 0; h2_full <- 0
    se_h2 <- NA_real_
  } else {
    A  <- res$alpha[idx, , drop = FALSE]
    M  <- res$mu[idx, , drop = FALSE]
    M2 <- res$mu2[idx, , drop = FALSE]
    Wl <- A * M                       # L x p, per-component mean contribution
    mu <- colSums(Wl)                 # posterior mean per SNP  (= coef(res)[-1])
    within <- sum(A * M2)             # sum_l s_l
    Rmu <- R %*% mu
    RWt <- R %*% t(Wl)                # p x L
    self_R <- sum(t(Wl) * RWt)        # sum_l m_l' R m_l
    self_I <- sum(Wl * Wl)            # sum_l m_l' m_l
    h2_full <- as.numeric(crossprod(mu, Rmu)) - self_R + within

    sum_mu2 <- sum(mu^2)
    sum_sd2 <- within - self_I        # == sum_j Var(b_j), the susieR quantity
    # R = I reduction, checked every fit rather than argued for once
    trace_from_full <- sum(mu^2) - self_I + within
    if (abs(trace_from_full - (sum_mu2 + sum_sd2)) > 1e-8)
      stop("trace/full identity violated")

    sd <- sqrt(pmax(colSums(A * M2 - Wl^2), 0))
    e_b2 <- mu^2 + sd^2
    e_b4 <- mu^4 + 6 * mu^2 * sd^2 + 3 * sd^4
    se_h2 <- sqrt(sum(pmax(e_b4 - e_b2^2, 0)))
  }

  cs_sizes <- if (length(res$sets$cs)) vapply(res$sets$cs, length, integer(1)) else integer(0)

  set(out, k, "susie_h2", sum_mu2 + sum_sd2)
  set(out, k, "susie_h2_se", se_h2)
  set(out, k, "susie_h2_quad", h2_full)
  set(out, k, "susie_sum_mu2", sum_mu2)
  set(out, k, "susie_sum_sd2", sum_sd2)
  set(out, k, "susie_n_cs", length(cs_sizes))
  set(out, k, "susie_cs_size", as.integer(sum(cs_sizes)))
  set(out, k, "susie_max_pip", max(res$pip))
  set(out, k, "susie_sum_pip", sum(res$pip))
  set(out, k, "susie_niter", as.integer(if (is.null(res$niter)) NA else res$niter))
  set(out, k, "susie_sigma2", as.numeric(if (is.null(res$sigma2)) NA else res$sigma2))
  set(out, k, "susie_converged", isTRUE(res$converged))
  set(out, k, "susie_status", if (isTRUE(res$converged)) "ok" else "not_converged")
}

fwrite(out, file.path(dir, "susie.tsv"), sep = "\t", quote = FALSE, na = "NA")
