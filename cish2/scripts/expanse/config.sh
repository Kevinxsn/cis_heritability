#!/usr/bin/env bash
# Expanse locations for cish2.  Sourced by the sbatch scripts in this folder.
export PROJ=/home/sux002/work/cis_heritability
export DATA=/expanse/lustre/projects/ddp412/sux002/cis_heritability_data   # 1000G.EUR.<chr> + gtex/ (symlinks)
export VCF_DIR=/expanse/lustre/scratch/sux002/temp_project/1000g_phase3_vcf             # raw downloads, transient (scratch is purged)
export GTEX_SHARED=/expanse/lustre/projects/ddp412/tamariutabartell/gtex   # READ-ONLY, shared by the lab
export CISH2_ENV=/home/sux002/miniforge3/envs/cish2
export LOGS=$PROJ/logs
