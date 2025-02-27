#!/bin/bash
#SBATCH -p skylake-himem #skylake-himem #icelake #skylake #-himem #cclake
#SBATCH -A gottgens-sl2-cpu
### Modify this according to your Ray workload.
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --tasks-per-node=1
#SBATCH --time 22:00:00
### Modify this according to your Ray workload.
#SBATCH --cpus-per-task=30
#SBATCH --mem-per-cpu=5GB
#SBATCH --gpus-per-task=0
### outputs
#SBATCH --output logs/01-log-%J.txt
#SBATCH --error logs/01-err-%J.txt

python 01_regions_to_genes.py