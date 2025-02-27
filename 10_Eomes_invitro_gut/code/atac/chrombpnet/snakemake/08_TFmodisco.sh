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

modisco motifs -i H5PY \
    -n 1000000 \
    -op ../../../results/atac/predictions/$1/predictions/tfmodisco/ \
    -w 500 \
    [-v]

# OR: (more likely)
modisco motifs \
    -s ../../../results/atac/chrombpnet/$1/predictions/mean_shap_scores.onehot.h5 \
    -a ../../../results/atac/chrombpnet/$1/predictions/mean_shap_scores.h5 \
    -n 1000000 \
    -w 500 \
    -o ../../../results/atac/chrombpnet/$1/predictions/tfmodisco/ \
    -v
    
# https://github.com/kundajelab/neuro-variants/blob/main/src/run_modisco.jobscript.sh
# https://github.com/kundajelab/neuro-variants/blob/main/preprocess/5_run_modisco/domcke_2020/0_run_modisco.w500.sh
#