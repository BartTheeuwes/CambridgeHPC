metrics <- fread("/Users/ricard/data/gastrulation_multiome_10x/multiome1/original/per_barcode_metrics.csv")# %>%
  # .[is_cell==TRUE]

foo <- atac.metrics %>%
  .[,c("atac_mitochondrial_reads","atac_fragments","atac_raw_reads")] %>%
  .[,mt_ratio:=atac_mitochondrial_reads/atac_raw_reads] 

gghistogram(foo$mt_ratio, bins=60) +
  scale_x_continuous(limits=c(0,0.25))

metrics[grep("CGTGCTGCATATAACC",barcode)]

head(metrics$barcode)
