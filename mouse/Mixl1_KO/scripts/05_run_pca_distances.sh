#!/bin/bash
#SBATCH -p clincloud # skylake # clincloud-express
#SBATCH -A gottgens-ccld-sl2-cpu  # GOTTGENS-CCLD-SL2-CPU
#SBATCH -N 1
#SBATCH -n 12
#SBATCH --time 9:00:00
#SBATCH --job-name distances
#SBATCH --output dynamics-log-%J.txt

python3.8 pca_distances2.py