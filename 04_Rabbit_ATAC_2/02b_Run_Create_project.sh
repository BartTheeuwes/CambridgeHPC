#!/bin/bash
#SBATCH -p skylake-himem
#SBATCH -A gottgens-sl2-cpu
#SBATCH -N 1
#SBATCH -n 12
#SBATCH --time 11:00:00
#SBATCH --job-name 02b
#SBATCH --output 02b_Create_project-log-%J.txt
#SBATCH -e 02b-err-%j-%N.err

Rscript /rds/project/rds-SDzz0CATGms/users/bt392/04_Rabbit_ATAC_2/02b_Create_project.r