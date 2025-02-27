#!/bin/bash
#SBATCH -p skylake-himem
#SBATCH -A gottgens-sl2-cpu
#SBATCH -N 1
#SBATCH -n 12
#SBATCH --time 11:00:00
#SBATCH --job-name milo
#SBATCH --output milo-log-%J.txt
Rscript tdTom_MILO.R --stages E7.5 E8.5 E9.5