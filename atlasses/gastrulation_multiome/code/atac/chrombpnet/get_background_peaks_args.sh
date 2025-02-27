#!/bin/bash

# Activate ChromBPNet conda environment
#eval "$(conda shell.bash hook)"
#conda activate chrombpnet2

chrombpnet prep nonpeaks -g /rds/project/rds-SDzz0CATGms/references/10x/refdata-cellranger-arc-mm10-2020-A-2.0.0/fasta/genome.fa -p ../../../results/atac/chrombpnet/all_peaks.bed -c  ../../../results/atac/chrombpnet/mm10.chrom.sizes -fl ../../../results/atac/chrombpnet/fold_0.json -br ../../../results/atac/chrombpnet/blacklist.bed.gz -o ../../../results/atac/chrombpnet/all