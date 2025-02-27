input_bam="/Users/ricard/data/gastrulation_multiome_10x/multiome1/original/bam/atac_bam_chr19.bam"
output_bw="/Users/ricard/data/gastrulation_multiome_10x/multiome1/original/bam/atac_bam_chr19.bw"

samtools index $input_bam
bamCoverage --bam $input_bam -o $output_bw --binSize 100 --numberOfProcessors 2 --verbose 

# to-do: 
# (1) fragments are known to map contiguously, should be processed with read extension (--extendReads [INTEGER]).
# (2) use --minMappingQuality

