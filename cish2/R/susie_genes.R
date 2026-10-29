#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Summary-statistic fine-mapping for a batch of *genes*.
#
# The simulation script (susie_batch.R) shares one LD matrix across every cell
# of a region.  Real genes each have their own cis-window, so each carries its
# own LD -- but paying R startup and library loading once per gene would cost
# more than the fitting does.  This script is therefore invoked once per chunk
# of genes and loops over them.
#
# Each gene arrives as a single little-endian float64 file laid out as
#     [ p*p LD ][ p bhat ][ p shat ]
# which keeps the handoff to one file per gene instead of three.
#
# Inputs (inside --dir):
#   cells.tsv    cell_id, n, p, file   (one row per gene)
#   params.tsv   key, value
# Output:
#   susie.tsv    one row per gene
# ---------------------------------------------------------------------------

suppressMessages(library("susieR"))
suppressMessages(library("data.table"))

args <- commandArgs(trailingOnly = TRUE)
dir <- NULL
for (i in seq_along(args)) if (args[i] == "--dir") dir <- args[i + 1L]
if (is.null(dir)) stop("--dir <chunk scratch dir> is required")

params <- fread(file.path(dir, "params.tsv"), header = TRUE, sep = "\t")
pget <- function(key, default = NA) {
  hit <- params$value[params$key == key]
  if (length(hit) == 0) default else hit[1]
}
L <- as.integer(pget("susie_L", 10))
max_iter <- as.integer(pget("susie_max_iter", 100))
coverage <- as.numeric(pget("susie_coverage", 0.95))

cells <- fread(file.path(dir, "cells.tsv"), header = TRUE, sep = "\t")
K <- nrow(cells)

out <- data.table(
  cell_id = cells$cell_id,
  susie_h2 = NA_real_, susie_h2_se = NA_real_,
  susie_sum_mu2 = NA_real_, susie_sum_sd2 = NA_real_,
  susie_n_cs = NA_integer_, susie_cs_size = NA_integer_,
  susie_max_pip = NA_real_, susie_sum_pip = NA_real_,
  susie_niter = NA_integer_, susie_sigma2 = NA_real_,
  susie_converged = FALSE, susie_status = "not_run"
)

for (k in seq_len(K)) {
  p <- as.integer(cells$p[k])
  path <- file.path(dir, cells$file[k])
  total <- as.double(p) * p + 2 * p

  blob <- tryCatch(
    readBin(path, what = "double", n = total, size = 8, endian = "little"),
    error = function(e) e
  )
  if (inherits(blob, "error") || length(blob) != total) {
    set(out, k, "susie_status", "input_read_failed")
    next
  }

  R <- matrix(blob[seq_len(as.double(p) * p)], nrow = p, ncol = p)
  bhat <- blob[(as.double(p) * p + 1):(as.double(p) * p + p)]
  shat <- blob[(as.double(p) * p + p + 1):total]
  rm(blob)

  # check_prior left at its default TRUE: it is the guard that turns a runaway
  # IBSS fit into a clean error rather than an absurd heritability
  res <- tryCatch(
    withCallingHandlers(
      susie_rss(bhat = bhat, shat = shat, R = R, n = cells$n[k],
                L = L, max_iter = max_iter, coverage = coverage),
      warning = function(w) invokeRestart("muffleWarning")
    ),
    error = function(e) e
  )
  rm(R)

  if (inherits(res, "error")) {
    set(out, k, "susie_status", paste0("error: ", conditionMessage(res)))
    next
  }

  mu <- coef(res)[-1]
  sdv <- susie_get_posterior_sd(res)
  sum_mu2 <- sum(mu^2)
  sum_sd2 <- sum(sdv^2)
  e_b2 <- mu^2 + sdv^2
  e_b4 <- mu^4 + 6 * mu^2 * sdv^2 + 3 * sdv^4
  se_h2 <- sqrt(sum(pmax(e_b4 - e_b2^2, 0)))
  cs_sizes <- if (length(res$sets$cs)) vapply(res$sets$cs, length, integer(1)) else integer(0)

  set(out, k, "susie_h2", sum_mu2 + sum_sd2)
  set(out, k, "susie_h2_se", se_h2)
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
