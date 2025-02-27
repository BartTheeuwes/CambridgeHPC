#!/bin/bash
#SBATCH -p skylake # skylake-himem #cclake
#SBATCH -A gottgens-sl2-cpu
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time 24:00:00
#SBATCH --job-name Cuttag
#SBATCH --output Cuttag-log-%J.txt

TF=IgG
File=IgG_EKDL220006148-1a-N711-AK17192_HKHCTDSX3_L2

# Load modules and pointers
module load fastqc 
module load bowtie2-2.3.4.1-gcc-5.4.0-td2qanw
trimmomatic=/rds/project/rds-SDzz0CATGms/users/bt392/software/trimmomatic/
module load samtools

# fastqc on raw data
#fastqc -t 2 ${File}_1.fq.gz ${File}_2.fq.gz

# Trim adapter sequences
#mkdir trimmed_fq
#java -jar $trimmomatic/trimmomatic-0.38.jar PE -phred33 ${File}_1.fq.gz ${File}_2.fq.gz trimmed_fq/${TF}_1.fq.gz trimmed_fq/${TF}_unpaired_1.fq.gz trimmed_fq/${TF}_2.fq.gz trimmed_fq/${TF}_unpaired_2.fq.gz ILLUMINACLIP:$trimmomatic/adapters/NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# fastq trimmed files
#fastqc -t 2 trimmed_fq/${TF}_1.fq.gz trimmed_fq/${TF}_2.fq.gz

# Fastq mapping
#bowtie2 --local --very-sensitive-local --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700 -x /rds/project/rds-SDzz0CATGms/references/bioinformatics_resources/bowtie2_indexes_ucsc/mm10 -p 8 -1 trimmed_fq/${TF}_1.fq.gz -2 trimmed_fq/${TF}_2.fq.gz -S ${TF}.sam

# Sam to Bam
samtools sort -@ 8 ${TF}.sam > ${TF}_sorted.bam 

# Remove duplicates
java -jar /rds/project/rds-SDzz0CATGms/users/bt392/software/picard/picard.jar MarkDuplicates I=${TF}_sorted.bam O=${TF}_no_duplicates.bam REMOVE_DUPLICATES=true M=marked_dup_metrics.txt

# Sort bam for index
samtools index ${TF}_no_duplicates.bam

# Bam to bigwig
#https://deeptools.readthedocs.io/en/develop/content/feature/effectiveGenomeSize.html
bamCoverage -p 8 -b ${TF}_no_duplicates.bam -of bedgraph --binSize 10 --normalizeUsing RPGC --ignoreForNormalization chrX --effectiveGenomeSize 2494787188 --extendReads -o ${TF}.bedgraph
bamCoverage -p 8 -b ${TF}_no_duplicates.bam --binSize 10 --normalizeUsing RPGC --ignoreForNormalization chrX --effectiveGenomeSize 2494787188 --extendReads -o ${TF}.bw
