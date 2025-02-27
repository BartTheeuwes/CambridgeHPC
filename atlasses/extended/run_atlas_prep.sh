#!/bin/bash
#SBATCH -p skylake-himem # skylake
#SBATCH -A gottgens-sl2-cpu
#SBATCH -N 1
#SBATCH -n 14
#SBATCH --time 12:00:00
#SBATCH --job-name atlas
#SBATCH --output Rscript-log-%J.txt


# Leave this file here UN-MODDIFIED, if you want to make 
# use of it, create a copy for yourself in your own area
# so that you can edit it as required.


# ## start an ipcluster instance and launch jupyter server from my container

singularity exec /home/bt392/containers/jupyter.metztli.25h.sif Rscript /rds/project/rds-SDzz0CATGms/users/bt392/atlasses/extended/extended_atlas_prep.R