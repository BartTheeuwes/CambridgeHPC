#!/bin/bash
#SBATCH -p skylake
#SBATCH -A gottgens-sl2-cpu
#SBATCH -N 1
#SBATCH -n 24
#SBATCH --time 12:00:00
#SBATCH --job-name cell_calling
#SBATCH --output 02-log-%J.txt


# Leave this file here UN-MODDIFIED, if you want to make 
# use of it, create a copy for yourself in your own area
# so that you can edit it as required.


# ## start an ipcluster instance and launch jupyter server from my container

singularity exec /home/bt392/containers/jupyter.metztli.25h.sif Rscript /rds/project/rds-SDzz0CATGms/users/bt392/03_Stat3_RNA/02_cell_calling.R

