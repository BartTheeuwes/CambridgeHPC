#!/bin/bash
#SBATCH -p skylake # skylake-himem #cclake
#SBATCH -A gottgens-sl2-cpu
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time 24:00:00
#SBATCH --job-name bowtie2
#SBATCH --output bowtie2-log-%J.txt

TF=Pu1
File=Pu1_EKDL220006148-1a-AK11822-AK19404_HKHCTDSX3_L2

# Load modules and pointers
module load fastqc 
module load bowtie2-2.3.4.1-gcc-5.4.0-td2qanw
trimmomatic=/rds/project/rds-SDzz0CATGms/users/bt392/software/trimmomatic/
module load samtools
PATH=$PATH:/rds/project/rds-SDzz0CATGms/users/bt392/software/homer/bin
module load perl
module load bedtools

eval "$(conda shell.bash hook)"
conda activate deeptools

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
#samtools sort -@ 8 ${TF}.sam > ${TF}_sorted.bam 

# Remove duplicates
#java -jar /rds/project/rds-SDzz0CATGms/users/bt392/software/picard/picard.jar MarkDuplicates I=${TF}_sorted.bam O=${TF}_no_duplicates.bam REMOVE_DUPLICATES=true M=marked_dup_metrics.txt

# Sort bam for index
#samtools index ${TF}_no_duplicates.bam

# Bam to bigwig
#https://deeptools.readthedocs.io/en/develop/content/feature/effectiveGenomeSize.html
bamCoverage -p 8 -b ${TF}_no_duplicates.bam -of bedgraph --binSize 10 --normalizeUsing RPGC --ignoreForNormalization chrX --effectiveGenomeSize 2494787188 --extendReads -o ${TF}.bedgraph
bamCoverage -p 8 -b ${TF}_no_duplicates.bam --binSize 10 --normalizeUsing RPGC --ignoreForNormalization chrX --effectiveGenomeSize 2494787188 --extendReads -o ${TF}.bw

# peak calling
# Macs2 
base=/rds/project/rds-SDzz0CATGms/users/bt392/07_scRNA_CutTag/01_bulk/data/X204SC22051679-Z01-F001/raw_data/
macs2 callpeak -t ${TF}_no_duplicates.bam -c $base/IgG/IgG_no_duplicates.bam  -n ${TF} -p 1e-5 --outdir macs2

# Seacr
mkdir seacr
cd seacr
bash /rds/project/rds-SDzz0CATGms/users/bt392/software/seacr/SEACR_1.3.sh ../$TF.bedgraph $base/IgG/IgG.bedgraph norm relaxed $TF
cd ..

# Motif search using Homer
/rds/project/rds-SDzz0CATGms/users/bt392/software/homer/bin/findMotifsGenome.pl macs2/${TF}_summits.bed mm10 MotifOutput_macs2/ -size 150 -S 10
/rds/project/rds-SDzz0CATGms/users/bt392/software/homer/bin/findMotifsGenome.pl seacr/$TF.relaxed.bed mm10 MotifOutput_seacr/ -size 150 -S 10

# Clean up space
#rm ${TF}.sam
#rm ${TF}_sorted.bam 