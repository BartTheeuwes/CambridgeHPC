#!/bin/bash
#SBATCH -p ampere  
#SBATCH -A GOTTGENS-SL2-GPU
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --gres=gpu:1
#SBATCH --time 10:00:00
#SBATCH --job-name gpu
#SBATCH --output bias_model_log-%j.txt

# Run before installing everything and everytime when using the environment!
. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load rhel8/default-amp              # REQUIRED - loads the basic environment
module load cuda/11.2
module load cudnn/8.1_cuda-11.2

# Activate ChromBPNet conda environment
eval "$(conda shell.bash hook)"
conda activate chrombpnet3

# Modify to not need to recreate bigwigs!
chrombpnet bias pipeline \
    -ifrag ../../../results/atac/chrombpnet/fragments/$1_fragments.tsv.gz \
    -d "ATAC" \
    -g /rds/project/rds-SDzz0CATGms/references/10x/refdata-cellranger-arc-mm10-2020-A-2.0.0/fasta/genome.fa \
    -c ../../../results/atac/chrombpnet/mm10.chrom.sizes \
    -p ../../../results/atac/chrombpnet/all_peaks.bed \
    -n ../../../results/atac/chrombpnet/all_negatives.bed \
    -fl ../../../results/atac/chrombpnet/fold_$2.json \
    -b 0.7 \
    -o ../../../results/atac/chrombpnet/$3/ \
    --tmpdir ../../../results/atac/chrombpnet/tmp     
