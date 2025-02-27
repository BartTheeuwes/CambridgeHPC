#!/bin/bash

# script to process bam file to fragments file
# input is the bam file provided by Blanca from her 2020 manuscript

# module load samtools
# module load htslib
# module load python
# import sinto

cd /bi/scratch/Stephen_Clark/gastrulation_multiome_10x/public_datasets/Pijuan-Sala_2020/bam

# sort bam file by position 
samtools sort -o embryo_revision1_possorted.bam -O bam -@ 8 embryo_revision1_sorted.bam  


sinto fragments -b embryo_revision1_possorted.bam -f embryo_revision1_fragments.tsv --barcode_regex "[^:]*" -p 8 --use_chrom "[1-9x-yX-Y]"

sort -k 1,1 -k2,2n embryo_revision1_fragments.tsv > embryo_revision1_fragments_sorted.tsv

bgzip embryo_revision1_fragments_sorted.tsv

tabix -p bed embryo_revision1_fragments_sorted.tsv.gz