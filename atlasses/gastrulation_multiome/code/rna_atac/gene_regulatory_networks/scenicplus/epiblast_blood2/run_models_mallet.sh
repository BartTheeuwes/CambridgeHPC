#!/bin/bash
#SBATCH -p skylake-himem #skylake-himem #icelake #skylake #-himem #cclake
#SBATCH -A gottgens-sl2-cpu
### Modify this according to your Ray workload.
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --tasks-per-node=1
#SBATCH --time 12:00:00
### Modify this according to your Ray workload.
#SBATCH --cpus-per-task=30
#SBATCH --mem-per-cpu=5GB
#SBATCH --gpus-per-task=0
### outputs
#SBATCH --output logs/mallet-log-%J.txt
#SBATCH --error logs/mallet-err-%J.txt

python /rds/project/rds-SDzz0CATGms/users/bt392/software/github/pycisTopic/model_scripts/runModels_lda_mallet.py \
        -i /rds/project/rds-SDzz0CATGms/users/bt392/atlasses/gastrulation_multiome/results/rna_atac/gene_regulatory_networks/scenicplus/blood/pycistopic/cisTopicObject.pkl \
        -o /rds/project/rds-SDzz0CATGms/users/bt392/atlasses/gastrulation_multiome/results/rna_atac/gene_regulatory_networks/scenicplus/blood/pycistopic/models_500_mallet.pkl \
        -nt 2,4,10,15,25,35 \
        -c 30 \
        -it 500 \
        -a 50 \
        -abt True \
        -e 0.1 \
        -ebt False \
        -sp /home/bt392/ry/c/intermediate_models/ \
        -s 555 \
        -td /home/bt392/ry/c/