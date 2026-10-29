GENO_PREFIX=$1
TMPDIR=$2
CURRENT_DATE=$3


python simulation_low_ld_test.py \
  --geno_prefix $GENO_PREFIX \
  --tmpdir $TMPDIR \
  --current_date $CURRENT_DATE \
