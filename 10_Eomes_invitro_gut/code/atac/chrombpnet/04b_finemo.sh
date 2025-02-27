#!/bin/bash
#SBATCH -p ampere  
#SBATCH -A GOTTGENS-SL2-GPU
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --gres=gpu:1
#SBATCH --time 10:00:00
#SBATCH --job-name gpu
#SBATCH --output ChromBPNet_log-%j.txt

# Run before installing everything and everytime when using the environment!
. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load rhel8/default-amp              # REQUIRED - loads the basic environment
module load cuda/11.2
module load cudnn/8.1_cuda-11.2

# Activate ChromBPNet conda environment
eval "$(conda shell.bash hook)"
conda activate finemo

finemo extract-regions-chrombpnet-h5 \
    -c ../../../results/atac/chrombpnet/PS_EOd/auxiliary/interpret_subsample/chrombpnet_nobias.profile_scores.h5 \
    -o ../../../results/atac/chrombpnet/PS_EOd/finemo/shap.npz
    
finemo call-hits \
    -r ../../../results/atac/chrombpnet/PS_EOd/finemo/shap.npz \
    -m ../../../results/atac/chrombpnet/PS_EOd/auxiliary/interpret_subsample/modisco_results_profile_scores.h5 \
    -p ../../../results/atac/chrombpnet/PS_EOd/auxiliary/interpret_subsample/chrombpnet_nobias.interpreted_regions.bed \
    -b 128 \
    -a 0.7 \
    -o ../../../results/atac/chrombpnet/PS_EOd/finemo/
   
finemo report \
    -r ../../../results/atac/chrombpnet/PS_EOd/finemo/shap.npz \
    -m ../../../results/atac/chrombpnet/PS_EOd/auxiliary/interpret_subsample/modisco_results_profile_scores.h5 \
    -p ../../../results/atac/chrombpnet/PS_EOd/auxiliary/interpret_subsample/chrombpnet_nobias.interpreted_regions.bed \
    -W 500 \
    -H ../../../results/atac/chrombpnet/PS_EOd/finemo/hits.tsv \
    -n \
    -s \
    -o ../../../results/atac/chrombpnet/PS_EOd/finemo/
   