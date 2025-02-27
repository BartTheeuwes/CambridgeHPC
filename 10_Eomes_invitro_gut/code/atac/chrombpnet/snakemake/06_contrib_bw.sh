#!/bin/bash
#SBATCH -p ampere  
#SBATCH -A GOTTGENS-SL2-GPU
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --gres=gpu:1
#SBATCH --time 3:00:00
#SBATCH --job-name gpu
#SBATCH --output logs/ChromBPNet_contribs_log-%j.txt

# Run before installing everything and everytime when using the environment!
. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load rhel8/default-amp              # REQUIRED - loads the basic environment
module load cuda/11.2
module load cudnn/8.1_cuda-11.2

chrombpnet contribs_bw  -m ../../../results/atac/chrombpnet/$1/fold_$2/models/chrombpnet_nobias.h5 \
    -r ../../../results/atac/chrombpnet/all_peaks.bed \
    -g /rds/project/rds-SDzz0CATGms/references/10x/refdata-cellranger-arc-mm10-2020-A-2.0.0/fasta/genome.fa \
    -c ../../../results/atac/chrombpnet/mm10.chrom.sizes \
    -op ../../../results/atac/chrombpnet/$1/fold_$2/predictions/ \ 
    -pc profile