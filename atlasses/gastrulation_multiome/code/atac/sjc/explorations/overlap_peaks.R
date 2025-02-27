library(Seurat)
library(Signac)
library(purrr)
library(data.table)

library(ggplot2)
library(cowplot)
library(ggrepel)



source(here::here("settings.R"))

io$peaks         <- file.path(io$rawdata, "/processed/atac/signac/signac_peaks_per_celltype_predicted.rds")
#io$peaks         <- "/bi/scratch/Stephen_Clark/multiome/raw/processed/atac/signac/old_jan21/signac_peaks_per_cell_type.rds"
io$signac_meta   <- file.path(io$rawdata, "processed/atac/signac/signac_metadata.tsv.gz")
io$anno_dir      <- "/bi/scratch/Stephen_Clark/annotations/gastrulation"
io$anno_files    <- file.path(io$anno_dir, c("H3K27ac_distal_E7.5_Ect_intersect12.bed",
                                             "H3K27ac_distal_E7.5_End_intersect12.bed",
                                             "H3K27ac_distal_E7.5_Mes_intersect12.bed"))

opts$extend_peaks <- 1e4 # add this many bp to either end of each locus 
#opts$peak_groups  <- as.character(0:30)

meta <- fread(io$signac_meta)
celltypes <- meta[, unique(celltype.predicted)]

h3k27ac <- map(io$anno_files, fread) %>% 
  rbindlist() %>% 
  setnames(c("chr", "start", "end","strand", "id", "anno")) %>% 
  .[, chr := paste0("chr", chr)] %>% 
  .[, c("start", "end") := .(start - opts$extend_peaks, end + opts$extend_peaks)] %>% 
  .[, width := end - start] %>% 
  setkey(chr, start ,end)

h3k27ac[, summary(width)]

peaks <- readRDS(io$peaks) %>% 
  as.data.table() %>% 
  setnames("seqnames", "chr") %>% 
  .[peak_called_in %in% celltypes] %>% 
  .[, peak := paste0(chr, "-", start, "-", end)] %>% 
  setkey(chr, start, end)  
  

overlap <- foverlaps(peaks, h3k27ac)



proportions <- overlap[, .N, .(peak_called_in, anno)]

ggplot(proportions, aes(peak_called_in, N, fill = anno)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

rev_overlap <- foverlaps(h3k27ac, peaks)
rev_props <- rev_overlap[, .N, .(peak_called_in, anno)]

ggplot(rev_props, aes(anno, N, fill = peak_called_in)) +
  geom_bar(stat = "identity", position = "fill")+
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# if the peaks are called on clusters we will need to asign labels to clusters

# asign_cluster <- function(celltypes, min_fraction = 0.5){
#   celltypes <- celltypes[!is.na(celltypes)]
#   n <- length(celltypes)
#   order <- table(celltypes)%>%
#     .[order(-rank(.))]
#   if (order[[1]] > n * min_fraction) {
#     return(names(order)[1])
#   } else {
#     return("no celltype assigned")
#   }
# }
# 
# celltypes <- meta[!is.na(celltype.denoised), .(seurat_clusters, celltype.denoised)] %>% 
#   .[, .(celltype = asign_cluster(celltype.denoised)), seurat_clusters]
# 
# umap <- meta[, .(cell, seurat_clusters, UMAP_1, UMAP_2)] %>% 
#   merge(celltypes, by = "seurat_clusters")
# 
# ggplot(umap, aes(UMAP_1, UMAP_2, colour = celltype)) +
#   geom_point(alpha = 0.5) +
#   theme_cowplot()







