#!/bin/bash
#SBATCH -p icelake-himem #skylake-himem #icelake #skylake #-himem #cclake
#SBATCH -A gottgens-sl2-cpu
### Modify this according to your Ray workload.
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --tasks-per-node=1
#SBATCH --time 22:00:00
### Modify this according to your Ray workload.
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=5GB
#SBATCH --gpus-per-task=0
### outputs
#SBATCH --output logs/scenic-%J-log.txt
#SBATCH --error logs/scenic-%J-err.txt

python 02a_scenic.py