#!/bin/bash
#SBATCH -p cclake-himem #-himem #skylake-himem #icelake #skylake #-himem #cclake
#SBATCH -A gottgens-sl2-cpu
#SBATCH -N 1
#SBATCH -n 24
#SBATCH --time 10:00:00
#SBATCH --job-name RNA_trajectory
#SBATCH --output logs/RNA_trajectory.txt

Rscript RNA_trajectory.R