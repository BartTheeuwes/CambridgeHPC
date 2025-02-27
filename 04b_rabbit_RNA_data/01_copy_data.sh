#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -p skylake
#SBATCH -A gottgens-sl2-cpu
#SBATCH --time=0-20:59:59

cp $1 $2

