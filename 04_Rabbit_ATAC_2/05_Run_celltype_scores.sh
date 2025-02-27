#!/bin/bash
#SBATCH -p skylake-himem
#SBATCH -A gottgens-sl2-cpu
#SBATCH -N 1
#SBATCH -n 12
#SBATCH --time 11:00:00
#SBATCH --job-name 05
#SBATCH --output 05_multiome-log-%J.txt
#SBATCH -e 05-err-%j-%N.err

Rscript /rds/project/rds-SDzz0CATGms/users/bt392/04_Rabbit_ATAC_2/05_celltype_scores.r