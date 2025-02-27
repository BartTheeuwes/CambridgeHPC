#!/bin/bash
#SBATCH -p icelake-himem #-himem #skylake-himem #icelake #skylake #-himem #cclake
#SBATCH -A gottgens-sl2-cpu
#SBATCH -N 1
#SBATCH -n 76
#SBATCH --time 10:00:00
#SBATCH --job-name ATAC_trajectory
#SBATCH --output logs/ATAC_trajectory.txt

Rscript ATAC_trajectory.R