library(Seurat)
library(Signac)
library(purrr)
library(data.table)


library(ggplot2)
library(cowplot)

source(here::here("settings.R"))





io$anno      <- file.path(io$rawdata, "/processed/atac/signac/k3k27ac/h3k27ac_anno.rds")
io$signac    <- file.path(io$rawdata, "/processed/atac/signac/k3k27ac/signac_h3k27ac.rds")
io$plot_data <- file.path(io$rawdata, "/processed/atac/signac/k3k27ac/mean_h3k27ac.tsv.gz")

signac <- readRDS(io$signac)
signac

meta <- fread(io$metadata) %>% 
  .[, .(cell, celltype.denoised)]

anno <- readRDS(io$anno) %>% 
  as.data.table() %>% 
  setnames("seqnames", "chr") %>% 
  .[, locus := paste0(chr, "-", start, "-", end)] %>% 
  .[, .(locus, anno)]

atac <- signac@assays$h3k27ac@data %>% 
  as.data.table(keep.rownames = "locus") %>% 
  melt(id.vars = "locus", variable.name = "cell", value.name = "counts") %>% 
  merge(anno, by = "locus") %>% 
  .[, .(mean_counts = mean(counts)), .(cell, anno)] %>% 
  merge(meta, by = "cell")

fwrite(atac, io$plot_data, sep = "\t", na = "NA", quote = FALSE)
atac=fread(io$plot_data)

hline <- atac[, .(mean = mean(mean_counts)), anno]

ggplot(atac, aes(celltype.denoised, mean_counts, fill = celltype.denoised)) +
  #geom_violin(alpha = 0.75) +
  geom_boxplot(alpha = 0.75, outlier.shape = NA) +
  facet_wrap(~anno, ncol = 1) +
  geom_hline(data = hline, aes(yintercept = mean), linetype = "dashed") +
  theme_cowplot() +
  guides(fill = FALSE) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))





