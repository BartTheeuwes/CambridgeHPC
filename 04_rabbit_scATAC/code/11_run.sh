#!/bin/bash
#SBATCH -p skylake-himem
#SBATCH -A gottgens-sl2-cpu
#SBATCH -N 1
#SBATCH -n 30
#SBATCH --time 11:00:00
#SBATCH --job-name snakemake
#SBATCH --output 11_TF-%J.txt
Rscript "/rds/project/rds-SDzz0CATGms/users/bt392/04_Rabbit_ATAC_final/rerun/code/11_TF_targets.R"