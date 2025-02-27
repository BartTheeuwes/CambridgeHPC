#!/bin/bash
#SBATCH -p cclake #-himem #skylake-himem #icelake #skylake #-himem #cclake
#SBATCH -A gottgens-sl2-cpu
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --time 23:00:00
#SBATCH --job-name scVI
#SBATCH --output logs/scVI-log-%J.txt

python 