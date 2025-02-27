
samples=( "E7.75_rep1" "E8.5_CRISPR_T_KO" "E8.5_CRISPR_T_WT" )
samples=( "E7.75_rep1" )
samples=( "E8.5_CRISPR_T_KO" "E8.5_CRISPR_T_WT" )

# Extract a single chromosome for testing
# samtools view -b possorted_genome_bam.bam chr1 > possorted_genome_chr1.bam

for i in "${samples[@]}"; do
	echo $i
	cd $i
	sbatch -n 4 --mem 45G --wrap "samtools index -@ 4 possorted_genome_bam.bam"	
	cd ..
done

for i in "${samples[@]}"; do
	echo $i
	cd $i/outs
	sbatch -n 4 --mem 70G --wrap "samtools sort -@ 4 -m 10G -t CB -O BAM -o cellsorted_possorted_genome_bam.bam possorted_genome_bam.bam"
	cd ../..
done
