#!/usr/bin/env bash
# Environment for every cish2 entry point.  Source, don't execute.
#
#   source scripts/env.sh
#
# The dsc180 conda env carries GCTA 1.94.1, PLINK 1.90, R + susieR and the
# Python scientific stack.  BLAS threads are pinned to 1 because the driver
# gets its parallelism from many single-threaded worker processes, and nested
# threading would oversubscribe the node.

export CISH2_ENV="${CISH2_ENV:-/home/sux002/miniforge3/envs/cish2}"
export PATH="${CISH2_ENV}/bin:${PATH}"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

CISH2_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export CISH2_ROOT
export PYTHONPATH="${CISH2_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
