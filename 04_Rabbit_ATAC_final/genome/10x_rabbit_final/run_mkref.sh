#!/bin/bash
#SBATCH -p skylake-himem
#SBATCH -A gottgens-sl2-cpu
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time 11:00:00
#SBATCH --job-name snakemake
#SBATCH --output snakemake-log-%J.txt
cellranger-atac mkref --config=config