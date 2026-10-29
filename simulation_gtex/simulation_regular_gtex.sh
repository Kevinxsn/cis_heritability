#!/bin/bash
#PBS -q hotel
#PBS -l nodes=1:ppn=1
#PBS -l walltime=24:00:00
#PBS -N gtex_h2_analysis
#PBS -o gtex_h2_analysis.out
#PBS -e gtex_h2_analysis.err

# ── Configuration ────────────────────────────────────────────────────────────
# plink and gcta64 are found via conda environment (no hardcoded paths needed)

# GTEx data paths
GTEX_DIR=/projects/ps-renlab2/sux002/data/gtex
GENO_PREFIX=${GTEX_DIR}/GTEx_v8_genotype_EUR_HM3_exclude_dups.allchr.reorder
EXPRESSION_BED=${GTEX_DIR}/Whole_Blood.v8.normalized_expression.bed.gz
# Optional: uncomment to include covariates (not yet implemented in pipeline)
# COVARIATE_FILE=${GTEX_DIR}/Covar_all_Whole_Blood.txt

# Output configuration
CURRENT_DATE=$(date +%Y_%m_%d)
TMPDIR=/projects/ps-renlab2/sux002/output/tmpdir_gtex_whole_blood
N_GENES=5000

# ── Setup ────────────────────────────────────────────────────────────────────
mkdir -p ${TMPDIR}

echo "============================================"
echo "GTEx Heritability Analysis: SuSiE vs GCTA"
echo "============================================"
echo "Date:           ${CURRENT_DATE}"
echo "Tissue:         Whole_Blood"
echo "Genotype:       ${GENO_PREFIX}"
echo "Expression:     ${EXPRESSION_BED}"
echo "Temp dir:       ${TMPDIR}"
echo "N genes:        ${N_GENES}"
echo "============================================"

# ── Run ──────────────────────────────────────────────────────────────────────
# No --plink or --gcta flags needed; the Python script defaults to
# "plink" and "gcta64" and finds them on PATH via conda
python simulation_regular_gtex.py \
    --geno_prefix ${GENO_PREFIX} \
    --expression_bed ${EXPRESSION_BED} \
    --tmpdir ${TMPDIR} \
    --current_date ${CURRENT_DATE} \
    --n_genes ${N_GENES}

echo "============================================"
echo "Analysis complete."
echo "============================================"

# ── Cleanup (optional — uncomment to remove temp files) ──────────────────────
# echo "Cleaning up temp files..."
# rm -rf ${TMPDIR}