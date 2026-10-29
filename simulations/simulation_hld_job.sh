PLINK_PATH=/projects/ps-renlab2/sux002/plink/plink
GENO_PREFIX=/projects/ps-renlab2/sux002/data/1000G.EUR.
TMPDIR=/projects/ps-renlab2/sux002/output/tmpdir_1000G_eur_regular_high_ld
CURRENT_DATE=2026_02_23_1

echo $GENO_PREFIX
echo $TMPDIR
echo $CURRENT_DATE

bash simulation_hld_job_helper.sh $GENO_PREFIX $TMPDIR $CURRENT_DATE

cd /projects/ps-renlab2/sux002/output/tmpdir_1000G_eur_regular_high_ld/
rm -rf -- *