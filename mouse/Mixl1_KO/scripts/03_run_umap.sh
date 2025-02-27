#!/bin/bash
#SBATCH -p clincloud-himem # skylake # clincloud-express
#SBATCH -A gottgens-ccld-sl2-cpu  # GOTTGENS-CCLD-SL2-CPU
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --time 24:00:00
#SBATCH --job-name tsne
#SBATCH --output tsne.txt


# Leave this file here UN-MODDIFIED, if you want to make
# use of it, create a copy for yourself in your own area
# so that you can edit it as required.

Rscript run_umap.R