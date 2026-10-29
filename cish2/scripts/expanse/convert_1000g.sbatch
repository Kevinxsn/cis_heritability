#!/usr/bin/env bash
# VCF -> PLINK 1000G.EUR.<chr> in $DATA, several chromosomes in parallel on one node.
# Filter chosen by calibrate_1000g (2026-10-28): every SNP polymorphic in the 503
# EUR samples, minus the one chr22 site the old panel lacked.  That reproduces the
# old chr22 (308,143 SNPs; all 17,287 in the reference window match on pos+alleles).
#   sbatch scripts/expanse/convert_1000g.sbatch            # chr1-22
#   sbatch --export=ALL,CHROMS="21 22" scripts/expanse/convert_1000g.sbatch
#SBATCH --job-name=cish2_conv1kg
#SBATCH --account=sdp202
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --output=/home/sux002/work/cis_heritability/logs/%x_%j.out
set -euo pipefail
source /home/sux002/work/cis_heritability/cish2/scripts/expanse/config.sh
source "$PROJ/cish2/scripts/env.sh"
CHROMS="${CHROMS:-$(seq -s ' ' 1 22)}"
export SELECT="--freq mac>=1 --exclude-ids $DATA/calib/exclude_ids.txt"
export W="${SLURM_TMPDIR:-/scratch/$USER/job_$SLURM_JOB_ID}"

convert() {
    local c=$1 d="$W/chr$1"; mkdir -p "$d"
    [ -s "$DATA/1000G.EUR.$c.bed" ] && { echo "skip chr$c (exists)"; return 0; }
    plink --vcf "$VCF_DIR/ALL.chr$c.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz" \
          --double-id --keep "$VCF_DIR/EUR.keep" \
          --snps-only just-acgt --biallelic-only strict \
          --make-bed --freq --threads 1 --memory 3000 --out "$d/raw" > /dev/null
    # shellcheck disable=SC2086
    python -I "$PROJ/cish2/scripts/expanse/select_variants.py" "$d/raw" "$d/keep.txt" $SELECT
    plink --bfile "$d/raw" --extract "$d/keep.txt" --make-bed --threads 1 --memory 3000 \
          --out "$d/1000G.EUR.$c" > /dev/null
    for ext in bed bim fam log; do cp "$d/1000G.EUR.$c.$ext" "$DATA/.1000G.EUR.$c.$ext.tmp"; done
    for ext in bim fam log bed; do mv "$DATA/.1000G.EUR.$c.$ext.tmp" "$DATA/1000G.EUR.$c.$ext"; done
    echo "done chr$c: $(wc -l < "$DATA/1000G.EUR.$c.bim") SNPs, $(wc -l < "$DATA/1000G.EUR.$c.fam") people"
    rm -rf "$d"
}
export -f convert
# largest chromosomes first so the long ones start early
echo $CHROMS | tr ' ' '\n' | sort -n | xargs -P 8 -I{} bash -c 'set -euo pipefail; convert {}'
