#!/usr/bin/env bash
# Round 4: full GTEx split-half pass in phenotype_mode=simulated.
# Resumable: re-submitting computes only chunks missing from manifest.jsonl.
#   sbatch scripts/expanse/gtex_sim.sbatch
#SBATCH --job-name=cish2_gtex_sim
#SBATCH --account=sdp202
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=120G
#SBATCH --time=00:45:00
#SBATCH --output=/home/sux002/work/cis_heritability/logs/round4/%x_%j.out
set -euo pipefail
source /home/sux002/work/cis_heritability/cish2/scripts/expanse/config.sh
source "$PROJ/cish2/scripts/env.sh"
cd "$PROJ/cish2"
python -m cish2 gtex configs/gtex_sim_split_half.yaml -o workers=$SLURM_CPUS_PER_TASK
