#!/bin/bash
#SBATCH -p skylake
#SBATCH -A gottgens-sl2-cpu
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time 12:00:00
#SBATCH --job-name Rscript
#SBATCH --output 09-log-%J.txt
#SBATCH -e 09-err-%j-%N.err

#singularity exec /home/bt392/containers/jupyter.metztli.25h.sif Rscript /rds/project/rds-SDzz0CATGms/users/bt392/01_Eomes_RNA/09_determine_DEGs.R
Rscript /rds/project/rds-SDzz0CATGms/users/bt392/01_Eomes_RNA/09_determine_DEGs.R