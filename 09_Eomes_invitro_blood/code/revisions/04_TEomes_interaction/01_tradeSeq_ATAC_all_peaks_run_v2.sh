#!/bin/bash
#SBATCH -p icelake-himem #-himem #skylake-himem #icelake #skylake #-himem #cclake
#SBATCH -A gottgens-sl2-cpu
#SBATCH -N 1
#SBATCH -n 42
#SBATCH --time 07:59:00
#SBATCH --job-name 01_tradeSeq_all_peaks
#SBATCH --output logs/01_tradeSeq_all_peaks_output_v2.txt

Rscript 01_tradeSeq_ATAC_all_peaks_v2.R