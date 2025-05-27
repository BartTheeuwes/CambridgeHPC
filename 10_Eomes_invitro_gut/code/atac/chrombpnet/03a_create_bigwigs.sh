#!/bin/bash
#SBATCH -p icelake-himem  
#SBATCH -A GOTTGENS-SL2-CPU
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 64
#SBATCH --time 10:00:00
#SBATCH --job-name create_bigwigs
#SBATCH --output create_bigwigs_log-%j.txt

chrombpnet createBW \
    -g /rds/project/rds-SDzz0CATGms/references/10x/refdata-cellranger-arc-mm10-2020-A-2.0.0/fasta/genome.fa \
    -ifrag ../../../results/atac/chrombpnet/fragments/fragments_small.tsv.gz \
    -c ../../../results/atac/chrombpnet/mm10.chrom.sizes \
    -d "ATAC" \
    -o ../../../results/atac/chrombpnet/reads_to_bigwig_test2 \
    --tmpdir ../../../results/atac/chrombpnet/tmp
