#!/bin/bash
#SBATCH -p skylake-himem #skylake-himem #icelake #skylake #-himem #cclake
#SBATCH -A gottgens-sl2-cpu
### Modify this according to your Ray workload.
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --mem 80G
#SBATCH --time 00:10:00
### outputs
#SBATCH --output logs/matrices-log-%J.txt
#SBATCH --error logs/matrices-err-%J.txt

Rscript prepare_matrices.r