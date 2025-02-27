#!/bin/bash
#SBATCH -p skylake-himem #-himem #skylake-himem #icelake #skylake #-himem #cclake
#SBATCH -A gottgens-sl2-cpu
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --time 10:00:00
#SBATCH --job-name seurat
#SBATCH --output seurat_test.txt

Rscript seurat_test.r