#!/bin/bash
#SBATCH -p skylake
#SBATCH -A gottgens-sl2-cpu
#SBATCH -N 1
#SBATCH -n 12
#SBATCH --time 36:00:00
#SBATCH --job-name souporcell
#SBATCH --output souporcell-log-%J.txt
#SBATCH -e souporcell-err-%j-%N.err

sample_folder="/rds/project/rds-SDzz0CATGms/users/bt392/03_Stat3_RNA/SLX-21143_fastq/SLX-21143_SITTF3_HTJH3DSX2/outs/"
bam="possorted_genome_bam.bam"
barcodes="filtered_feature_bc_matrix/barcodes.tsv.gz"

output="/rds/project/rds-SDzz0CATGms/users/bt392/03_Stat3_RNA/souporcell/mm10_t2"

referencefasta="/rds/project/rds-SDzz0CATGms/users/bt392/03_Stat3_RNA/souporcell/genome.fa"

gunzip -c $sample_folder$barcodes.gz > $sample_folder$barcodes

singularity exec /rds/project/rds-SDzz0CATGms/users/bt392/software/souporcell_latest.sif souporcell_pipeline.py -i $sample_folder$bam  -b $sample_folder$barcodes  -f $referencefasta  -t 8  -k 2 -o $output