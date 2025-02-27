#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=16
#SBATCH -p skylake-himem
#SBATCH -A gottgens-sl2-cpu
#SBATCH --time=1-11:59:59
#SBATCH -e cellranger-err-%j.%N.err
#SBATCH --output cellranger-log-%J.txt