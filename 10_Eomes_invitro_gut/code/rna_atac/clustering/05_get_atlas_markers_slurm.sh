#!/bin/bash
#SBATCH -p sapphire #-himem #skylake-himem #icelake #skylake #-himem #cclake
#SBATCH -A gottgens-sl2-cpu
#SBATCH -N 1
#SBATCH -n 46
#SBATCH --time 8:00:00
#SBATCH --job-name markers
#SBATCH --output 05_get_atlas_markers.txt

Rscript 05_get_atlas_markers.R