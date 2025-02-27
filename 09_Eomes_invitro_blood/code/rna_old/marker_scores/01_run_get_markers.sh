#!/bin/bash
#SBATCH -p cclake-himem #-himem #skylake-himem #icelake #skylake #-himem #cclake
#SBATCH -A gottgens-sl3-cpu
#SBATCH -N 1
#SBATCH -n 24
#SBATCH --time 11:00:00
#SBATCH --job-name 01_get_markers

Rscript 01_get_markers.R