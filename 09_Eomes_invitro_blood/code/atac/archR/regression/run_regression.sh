#!/bin/bash
#SBATCH -p cclake-himem #-himem #skylake-himem #icelake #skylake #-himem #cclake
#SBATCH -A gottgens-sl3-cpu
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --time 11:59:00
#SBATCH --job-name regression

Rscript regression.R