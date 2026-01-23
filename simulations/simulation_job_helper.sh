PLINK_PATH=$1
GENO_PREFIX=$2
TMPDIR=$3
CURRENT_DATE=$4


python simulations.py --plink_path $PLINK_PATH \
  --geno_prefix $GENO_PREFIX \
  --tmpdir $TMPDIR \
  --current_date $CURRENT_DATE \
