#!/bin/bash
#SBATCH -p icelake #-himem #skylake-himem #icelake #skylake #-himem #cclake
#SBATCH -A gottgens-sl2-cpu
#SBATCH -N 1
#SBATCH -n 64
#SBATCH --time 10:00:00
#SBATCH --job-name Trajectory
#SBATCH --output logs/Traj1_vs_traj2_ATAC.txt

Rscript traj1_vs_traj2_ATAC.R