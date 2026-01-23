PLINK_PATH=/projects/ps-renlab2/sux002/plink/plink
GENO_PREFIX=/projects/ps-renlab2/sux002/data/1000G.EUR.
TMPDIR=/projects/ps-renlab2/sux002/output/tmpdir_1000G_eur_regular
CURRENT_DATE=2026_01_22

echo $GENO_PREFIX
echo $TMPDIR
echo $CURRENT_DATE

bash simulation_job_helper.sh $PLINK_PATH $GENO_PREFIX $TMPDIR $CURRENT_DATE

