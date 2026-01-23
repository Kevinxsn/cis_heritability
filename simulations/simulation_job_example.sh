#!/bin/bash
#SBATCH --job-name="simulation_job"
#SBATCH --output="sim_02_22_25.%j.%N.out"
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2gb
#SBATCH --account=csd832
#SBATCH --export=ALL
#SBATCH -t 00:15:00

source ~/.mbashrc
mamba activate henv

PLINK_PATH=/expanse/lustre/projects/ddp412/sgolzari/plink/plink
GENO_PREFIX=/expanse/lustre/projects/ddp412/kakamatsu/ldref/eur/1000G_eur_chr
TMPDIR=tmpdir_1000G_eur_regular/
CURRENT_DATE=2025_02_22

echo $GENO_PREFIX
echo $TMPDIR
echo $CURRENT_DATE

# Create & submit 10 sbatch jobs
for i in {1..40}; do
    sbatch simulation_job_helper.sh $PLINK_PATH $GENO_PREFIX $TMPDIR $CURRENT_DATE
done
