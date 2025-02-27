# decompress files
# gzip -d E7.5_rep1_atac_fragments.tsv.gz E7.5_rep2_atac_fragments.tsv.gz E8.0_rep1_atac_fragments.tsv.gz E8.0_rep2_atac_fragments.tsv.gz E8.5_rep1_atac_fragments.tsv.gz E8.5_rep2_atac_fragments.tsv.gz
gunzip -f E7.5_rep2_atac_fragments.tsv.gz E8.0_rep1_atac_fragments.tsv.gz E8.0_rep2_atac_fragments.tsv.gz E8.5_rep1_atac_fragments.tsv.gz E8.5_rep2_atac_fragments.tsv.gz

# merge files (avoids having to re-sort)
sort -m -k 1,1V -k2,2n E7.5_rep1_atac_fragments.tsv E7.5_rep2_atac_fragments.tsv E8.0_rep1_atac_fragments.tsv E8.0_rep2_atac_fragments.tsv E8.5_rep1_atac_fragments.tsv E8.5_rep2_atac_fragments.tsv > fragments.tsv

# doesnt work???
# sort-bed --check-sort E7.5_rep1_atac_fragments.tsv E7.5_rep2_atac_fragments.tsv E8.0_rep1_atac_fragments.tsv E8.0_rep2_atac_fragments.tsv E8.5_rep1_atac_fragments.tsv E8.5_rep2_atac_fragments.tsv > fragments.tsv

# Remove trash chromosomes
# grep chr fragments.tsv > fragments2.tsv

# block gzip compress the merged file
bgzip -@ 4 fragments.tsv

# index the bgzipped file
tabix --preset=bed fragments.tsv.gz





# R script to rename barcodes in the fragments file
library(data.table)
library(purrr)
basedir <- "/hps/nobackup2/research/stegle/users/ricard/gastrulation_multiome_10x/original/fragments"
outdir <- "/hps/nobackup2/research/stegle/users/ricard/gastrulation_multiome_10x/original/fragments/test"
samples <- c("E7.5_rep1", "E7.5_rep2", "E8.0_rep1", "E8.0_rep2", "E8.5_rep1", "E8.5_rep2")
for (i in samples) {
	dt <- fread(sprintf("%s/%s_atac_fragments.tsv",basedir,i), header=F) %>% .[,V4:=sprintf("%s#%s",i,V4)] %>% .[grep("chr",V1)]
	table(dt$V1)
	fwrite(dt, sprintf("%s/%s_atac_fragments.tsv",outdir,i), col.names=F, quote=F, sep="\t")
}

