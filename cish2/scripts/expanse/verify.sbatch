#!/usr/bin/env bash
# MIGRATION.md §4.5 checks, in order, on a compute node.  Writes results/smoke_expanse/.
#   sbatch scripts/expanse/verify.sbatch
#SBATCH --job-name=cish2_verify
#SBATCH --account=sdp202
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=/home/sux002/work/cis_heritability/logs/%x_%j.out
set -euo pipefail
source /home/sux002/work/cis_heritability/cish2/scripts/expanse/config.sh
source "$PROJ/cish2/scripts/env.sh"
which python Rscript gcta64 plink

echo "== a) poster numbers (no data needed)"
cd "$PROJ/poster" && python check_numbers.py

echo "== panel) rebuilt 1000G vs regions recorded by the old runs"
python -I "$PROJ/cish2/scripts/expanse/check_panel.py" "$DATA/1000G.EUR." "$PROJ/results"

echo "== b) unit tests + SuSiE estimand"
cd "$PROJ/cish2"
python -m pytest tests/test_pipeline.py -q -p no:cacheprovider
Rscript tests/test_susie_estimand.R

echo "== c) smoke campaign (new run name so the archived results/smoke stays untouched)"
python -m cish2 plan configs/smoke.yaml -o name=smoke_expanse
python -m cish2 run  configs/smoke.yaml -o name=smoke_expanse -o compute.workers=4
python -m cish2 analyze "$PROJ/results/smoke_expanse"
