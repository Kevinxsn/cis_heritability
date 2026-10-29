#!/usr/bin/env bash
# Build the `cish2` conda env from environment.yml on a compute node.
#   sbatch scripts/expanse/build_env.sbatch
#SBATCH --job-name=cish2_env
#SBATCH --account=sdp202
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/home/sux002/work/cis_heritability/logs/%x_%j.out
set -euo pipefail
source /home/sux002/work/cis_heritability/cish2/scripts/expanse/config.sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

MAMBA=/home/sux002/miniforge3/bin/mamba
if [ -d "$CISH2_ENV" ]; then
    echo "env exists at $CISH2_ENV -- updating"
    $MAMBA env update -y -n cish2 -f "$PROJ/cish2/environment.yml"
else
    $MAMBA env create -y -f "$PROJ/cish2/environment.yml"
fi

export PATH="$CISH2_ENV/bin:$PATH"
echo "== tool check"
which python Rscript gcta64 plink
python -c "import numpy, scipy, pandas, statsmodels, pandas_plink, yaml, pptx; print('python ok', numpy.__version__, pandas.__version__)"
Rscript -e 'suppressMessages({library(susieR); library(data.table); library(Rfast)}); cat("R ok susieR", as.character(packageVersion("susieR")), "\n")'
gcta64 2>&1 | head -5 || true
plink --version
