#!/usr/bin/env bash
# Run one campaign, then build its tables and figures.
#
#   scripts/run_campaign.sh configs/baseline.yaml [extra cish2 args...]
#
# Safe to re-run: the campaign resumes from its manifest and only the missing
# region slots are computed.
set -euo pipefail

CONFIG="${1:?usage: run_campaign.sh <config.yaml> [extra args]}"
shift || true

source "$(dirname "${BASH_SOURCE[0]}")/env.sh"
cd "${CISH2_ROOT}"

python -m cish2 plan "${CONFIG}" "$@"
python -m cish2 run  "${CONFIG}" "$@"

RUN_DIR="$(python - "${CONFIG}" <<'PY'
import sys
from cish2.config import CampaignConfig
print(CampaignConfig.from_yaml(sys.argv[1]).run_dir)
PY
)"
python -m cish2 analyze "${RUN_DIR}"
echo "results: ${RUN_DIR}"
