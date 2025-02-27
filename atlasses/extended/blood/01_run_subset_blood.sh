#!/bin/bash
#SBATCH -p skylake-himem # skylake
#SBATCH -A gottgens-sl3-cpu
#SBATCH -N 1
#SBATCH -n 17
#SBATCH --time 12:00:00
#SBATCH --job-name atlas
#SBATCH --output 01-%J.txt
#SBATCH --error 01-%J.err

Rscript /rds/project/rds-SDzz0CATGms/users/bt392/atlasses/extended/blood/01_subset_blood.R