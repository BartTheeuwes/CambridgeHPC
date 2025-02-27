#!/bin/bash
#SBATCH -p cclake-himem #-himem #skylake-himem #icelake #skylake #-himem #cclake
#SBATCH -A gottgens-sl2-cpu
#SBATCH -N 1
#SBATCH -n 24
#SBATCH --time 10:00:00
#SBATCH --job-name homer

/rds/project/rds-SDzz0CATGms/users/bt392/software/homer/bin/findMotifsGenome.pl $1 mm10 $2 -size $3 -S 25 -len 6,8,10 -p 24