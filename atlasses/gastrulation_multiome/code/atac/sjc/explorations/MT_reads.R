library(ggplot2)
library(cowplot)
library(data.table)
library(purrr)

source(here::here("settings.R"))



io$merged_metrics_file    <- file.path(io$basedir, "/processed/atac/signac/cell_metrics.tsv")
#io$outfile <- "/bi/scratch/Stephen_Clark/gastrulation_multiome_10x//processed/atac/signac/reads_per_chr.tsv.gz"
io$chrom_sizes <- "http://hgdownload.cse.ucsc.edu/goldenpath/mm10/bigZips/mm10.chrom.sizes"


chrom_sizes <- tempfile(fileext = ".txt")
download.file(io$chrom_sizes, chrom_sizes)

chrom_sizes_dt <- fread(chrom_sizes)

mt_size <- chrom_sizes_dt[, sum := sum(V2)][V1 == "chrM", V2/sum]

metrics <- fread(io$merged_metrics_file) %>% 
  .[, mt_fraction := atac_mitochondrial_reads / atac_raw_reads]

ggplot(metrics, aes(sample, mt_fraction, colour = sample, fill = sample)) +
  # geom_point(position = position_jitterdodge(jitter.width = 0.1),
  #            colour = "black", alpha = 0.1) +
  #geom_violin() +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  guides(fill = FALSE, colour = FALSE) +
  ylim(0, 0.1)


stop()
reads_per_chr <- fread(io$merged_fragment_file, header = FALSE) %>% 
  setnames(c("chr", "start", "end", "cell", "reads")) %>% 
  .[, .(sum = sum(reads), .N), .(chr, cell) ]
  
fwrite(reads_per_chr, io$outfile, sep = "\t", na = "NA")
#reads_per_chr <- fread(io$outfile)

reads_per_chr[, cell_total := sum(sum), cell]

mt_reads <- reads_per_chr[chr == "chrMT"]
