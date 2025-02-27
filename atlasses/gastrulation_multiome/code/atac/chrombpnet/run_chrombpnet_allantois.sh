#!/bin/bash
#SBATCH -p ampere  
#SBATCH -A GOTTGENS-SL2-GPU
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 24
#SBATCH --gres=gpu:1
#SBATCH --time 15:00:00
#SBATCH --job-name gpu
#SBATCH --output ChromBPNet_allantois_log.txt

# Run before installing everything and everytime when using the environment!
. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load rhel8/default-amp              # REQUIRED - loads the basic environment
module load cuda/11.2
module load cudnn/8.1_cuda-11.2

sleep 1800 # sleep so that the background peaks finish before running this

chrombpnet pipeline \
        -ifrag ../../../results/atac/chrombpnet/allantois/allantois_fragments.tsv.gz \
        -d "ATAC" \
        -g /rds/project/rds-SDzz0CATGms/references/10x/refdata-cellranger-arc-mm10-2020-A-2.0.0/fasta/genome.fa \
        -c ../../../results/atac/chrombpnet/allantois/mm10.chrom.sizes \
        -p ../../../results/atac/chrombpnet/allantois/peaks.bed \
        -n ../../../results/atac/chrombpnet/allantois/output_negatives.bed \
        -fl ../../../results/atac/chrombpnet/fold_0.json \
        -b /rds/project/rds-SDzz0CATGms/users/bt392/software/chrombpnet_tutorial/bias_model/ENCSR868FGK_bias_fold_0.h5 \
        -o ../../../results/atac/chrombpnet/allantois/chrombpnet_model/