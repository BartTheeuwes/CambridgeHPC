#!/bin/bash
#SBATCH -p sapphire #-himem #skylake-himem #icelake #skylake #-himem #cclake
#SBATCH -A gottgens-sl2-cpu
#SBATCH -N 1
#SBATCH -n 76
#SBATCH --time 07:59:00
#SBATCH --job-name 01_tradeSeq
#SBATCH --output 01_tradeSeq_output.txt

Rscript 01_tradeSeq_ATAC.R