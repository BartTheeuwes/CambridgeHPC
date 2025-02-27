 #!/bin/bash
#SBATCH -p skylake
#SBATCH -A gottgens-sl2-cpu
#SBATCH -N 1
#SBATCH -n 22
#SBATCH --time 24:00:00
#SBATCH --job-name souporcell
#SBATCH --output souporcell-log-%J.txt
#SBATCH -e souporcell.%N.%j.err


sample_folder="/rds/project/rds-SDzz0CATGms/users/bt392/06_Runx1_RNA/data/SLX21184/genotyping/SLX-21184_SITTE12_H5NKNDMXY/outs/"
bam="possorted_genome_bam.bam"
barcodes="barcodes.tsv"
genotypes="g129S1_SvImJ C57BL_6NJ"

output="/rds/project/rds-SDzz0CATGms/users/bt392/06_Runx1_RNA/data/genotyping/SITTE12_known_genotype_inferred_from_unkown_genotype_skip_remap/"

referencefasta="/rds/project/rds-SDzz0CATGms/references/10x/refdata-gex-mm10-2020-A/fasta/genome.fa"

vcf="/rds/project/rds-SDzz0CATGms/users/bt392/06_Runx1_RNA/data/genotyping/SITTD7/cluster_genotypes.vcf"

gunzip -c $sample_folder/filtered_feature_bc_matrix/$barcodes.gz > $sample_folder$barcodes

# with supplied genotype
#singularity exec /rds/project/rds-SDzz0CATGms/users/bt392/software/souporcell_latest.sif souporcell_pipeline.py -i $sample_folder$bam  -b $sample_folder$barcodes  -f $referencefasta  -t 22  -k 2 -o $output  --known_genotypes $vcf   --known_genotypes_sample_names $genotypes

# skip remap
singularity exec /rds/project/rds-SDzz0CATGms/users/bt392/software/souporcell_latest.sif souporcell_pipeline.py -i $sample_folder$bam  -b $sample_folder$barcodes  -f $referencefasta  -t 1   -k 2 -o $output  --known_genotypes $vcf   --known_genotypes_sample_names $genotypes --skip_remap SKIP_REMAP
