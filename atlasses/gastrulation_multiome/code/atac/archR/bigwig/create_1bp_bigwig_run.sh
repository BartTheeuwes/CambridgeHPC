#!/bin/bash
#SBATCH -p icelake-himem  #-himem  #icelake #skylake #-himem #cclake
#SBATCH -A gottgens-sl2-cpu
#SBATCH -N 1
#SBATCH -n 72
#SBATCH --time 24:00:00
#SBATCH --job-name create_bigwig
#SBATCH --output log_bigwig.txt

Rscript create_1bp_bigwig.R