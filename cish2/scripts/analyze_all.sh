#!/usr/bin/env bash
# Rebuild every campaign's tables/figures/report, then the cross-campaign
# overview and the narrative write-up.
#
#   scripts/analyze_all.sh [results_root]
set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/env.sh"
cd "${CISH2_ROOT}"

RESULTS="${1:-/home/sux002/work/cis_heritability/results}"

for run in "${RESULTS}"/*/; do
  name="$(basename "${run}")"
  [[ "${name}" == smoke || "${name}" == figures || "${name}" == .* ]] && continue
  [[ -d "${run}/shards" ]] || continue
  echo "--- ${name}"
  python -m cish2 analyze "${run}"   # routes real-data runs to the GTEx report
done

python -m analysis.overview "${RESULTS}"
python -m analysis.findings "${RESULTS}"
