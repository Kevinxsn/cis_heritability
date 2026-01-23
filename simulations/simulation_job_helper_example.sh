#!/bin/bash
#SBATCH --job-name=simulation_job
#SBATCH --output="working_err/simulation_job.%j.%N.out"
#SBATCH --error="working_err/simulation_job.%j.%N.err"
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2gb
#SBATCH --account=csd832
#SBATCH --export=ALL
#SBATCH -t 10:00:00

source ~/.mbashrc
mamba activate henv

PLINK_PATH=$1
GENO_PREFIX=$2
TMPDIR=$3
CURRENT_DATE=$4

python simulations.py --plink_path $PLINK_PATH \
  --geno_prefix $GENO_PREFIX \
  --tmpdir $TMPDIR \
  --current_date $CURRENT_DATE \

