#!/bin/bash
#SBATCH -p skylake
#SBATCH -A gottgens-sl2-cpu
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time 24:00:00
#SBATCH --job-name souporcell
#SBATCH --output souporcell-log-%J.txt
#SBATCH -e souporcell.%N.%j.err

sample_folder=/rds/project/rds-SDzz0CATGms/users/bt392/01_Eomes_RNA/SLX-20795_fastq/SLX-20795_SITTA11_HKTG2DRXY/outs/
bam=possorted_genome_bam.bam
barcodes=filtered_feature_bc_matrix/barcodes.tsv.gz

output=/rds/project/rds-SDzz0CATGms/users/bt392/software/vcf/test/

referencefasta=/rds/project/rds-SDzz0CATGms/users/bt392/mouse/mm10_tomato/fasta/genome.fa

vcf=/rds/project/rds-SDzz0CATGms/users/bt392/software/vcf/C57BL6NJ_129S1.vcf


singularity exec /rds/project/rds-SDzz0CATGms/users/bt392/software/souporcell_latest.sif souporcell_pipeline.py -i ${sample_folder}${bam}  -b ${sample_folder}${barcodes}  -f ${referencefasta}  -t 8  -k 2 -o ${output}  --known_genotypes ${vcf}   --known_genotypes_sample_names 129S1_SvImJ C57BL_6NJ