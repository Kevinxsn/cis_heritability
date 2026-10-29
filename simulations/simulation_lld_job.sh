PLINK_PATH=/projects/ps-renlab2/sux002/plink/plink
GENO_PREFIX=/projects/ps-renlab2/sux002/data/1000G.EUR.
TMPDIR=/projects/ps-renlab2/sux002/output/tmdir_1000G_eur_regular_low_ld
CURRENT_DATE=2026_3_05_low_ld

echo $GENO_PREFIX
echo $TMPDIR
echo $CURRENT_DATE

bash simulation_lld_job_helper.sh $GENO_PREFIX $TMPDIR $CURRENT_DATE

cd /projects/ps-renlab2/sux002/output/tmdir_1000G_eur_regular_low_ld
rm -rf -- *