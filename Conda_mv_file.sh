#!/bin/bash
#SBATCH -p skylake-himem
#SBATCH -A gottgens-sl3-cpu
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time 12:00:00
#SBATCH --job-name mv_miniconda
#SBATCH --output conda_mv/conda_mv-%J.txt

# move file
cd /home/bt392/miniconda3/$1/
find . -maxdepth 1 -type f -exec mv {} /rds/project/rds-SDzz0CATGms/users/bt392/miniconda3/$1/ \;
