PLINK_PATH=/projects/ps-renlab2/sux002/plink/plink
GENO_PREFIX=/projects/ps-renlab2/sux002/data/1000G.EUR.
TMPDIR=/projects/ps-renlab2/sux002/output/tmdir_1000G_eur_regular_high_al_freq
CURRENT_DATE=2026_03_05_high_al_freq

echo $GENO_PREFIX
echo $TMPDIR
echo $CURRENT_DATE

bash simulation_high_al_freq_job_helper.sh $GENO_PREFIX $TMPDIR $CURRENT_DATE

cd /projects/ps-renlab2/sux002/output/tmdir_1000G_eur_regular_high_al_freq
rm -rf -- *
