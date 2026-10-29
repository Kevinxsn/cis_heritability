#!/usr/bin/env bash
# Round 4 (docs/pi_feedback_round4_design.md §3.5): unit + regression tests,
# then a 200-gene smoke run of each new GTEx mode.  Writes
# results/gtex_{sim,perm}_smoke/ only - nothing existing is touched.
#   sbatch scripts/expanse/round4_tests.sbatch
#SBATCH --job-name=cish2_r4_tests
#SBATCH --account=sdp202
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=60G
#SBATCH --time=00:30:00
#SBATCH --output=/home/sux002/work/cis_heritability/logs/round4/%x_%j.out
set -euo pipefail
source /home/sux002/work/cis_heritability/cish2/scripts/expanse/config.sh
source "$PROJ/cish2/scripts/env.sh"
cd "$PROJ/cish2"

echo "== unit + regression tests"
python -m pytest tests/test_pipeline.py -q -p no:cacheprovider -rs

echo "== smoke: simulated mode, 200 genes"
python -m cish2 gtex configs/gtex_sim_split_half.yaml -o name=gtex_sim_smoke \
    -o max_genes=200 -o workers=$SLURM_CPUS_PER_TASK --no-resume
echo "== smoke: permuted mode, 200 genes"
python -m cish2 gtex configs/gtex_perm_split_half.yaml -o name=gtex_perm_smoke \
    -o max_genes=200 -o workers=$SLURM_CPUS_PER_TASK --no-resume
