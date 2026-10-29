#!/usr/bin/env bash
# Convert chr22 EUR with no frequency filter and compare candidate filters
# against the surviving extract of the old panel.  Needs the env + downloads.
#   sbatch --dependency=afterok:<env>:<dl> scripts/expanse/calibrate_1000g.sbatch
#SBATCH --job-name=cish2_cal1kg
#SBATCH --account=sdp202
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/home/sux002/work/cis_heritability/logs/%x_%j.out
set -euo pipefail
source /home/sux002/work/cis_heritability/cish2/scripts/expanse/config.sh
source "$PROJ/cish2/scripts/env.sh"

W="$DATA/calib"; mkdir -p "$W"   # small; kept for inspection
plink --vcf "$VCF_DIR/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz" \
      --double-id --keep "$VCF_DIR/EUR.keep" --set-missing-var-ids @:#:\$1:\$2 \
      --snps-only just-acgt --biallelic-only strict \
      --make-bed --freq --threads 2 --out "$W/raw22" > /dev/null
grep -E "variants|people" "$W/raw22.log"
python -I "$PROJ/cish2/scripts/expanse/calibrate_1000g.py" "$W/raw22" \
    "$PROJ/expand_proposal/bench_regions/reg_22_kb1000.bim"
