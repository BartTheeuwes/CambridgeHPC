#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH -p skylake
#SBATCH -A gottgens-sl2-cpu
#SBATCH --time=0-20:59:59

tar -xvf $1 -C $2