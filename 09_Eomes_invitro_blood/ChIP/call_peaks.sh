#!/bin/bash
#SBATCH -p cclake-himem #-himem #skylake-himem #icelake #skylake #-himem #cclake
#SBATCH -A gottgens-sl2-cpu
#SBATCH -N 1
#SBATCH -n 24
#SBATCH --time 10:00:00
#SBATCH --job-name peak_calling

macs2 callpeak -t $1 -c $2 -g mm -q $3 -n $4 --outdir $5