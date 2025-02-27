#!/bin/bash
#SBATCH -p icelake-himem #-himem #skylake-himem #icelake #skylake #-himem #cclake
#SBATCH -A gottgens-sl2-cpu
#SBATCH -N 1
#SBATCH -n 72
#SBATCH --time 24:00:00
#SBATCH --job-name markers
#SBATCH --output 05_get_atlas_markers_ATAC.txt

Rscript 05_get_atlas_markers_ATAC.R