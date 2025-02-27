library(Seurat)
library(Signac)
library(purrr)
library(data.table)

library(GenomeInfoDb)
library(GenomicRanges)



source(here::here("settings.R"))



io$signac        <- file.path(io$rawdata, "/processed/atac/signac/signac_over_nmp_anno.rds")
io$outfile        <- file.path(io$rawdata, "/processed/atac/signac/pseudobulk/NMPs.tsv.gz")



opts$group_by    <- "celltype.predicted"#"seurat_clusters"#"celltype.predicted"#"seurat_clusters"#"celltype.mapped" # split data by this before peak calling (or NULL)

dir.create(dirname(io$outfile), recursive = TRUE)





signac <- readRDS(io$signac)
signac
signac@meta.data$sample %>% unique()

Idents(signac) <- opts$group_by
signac

celltypes <- unique(Idents(signac)) %>% 
  .[!is.na(.)]


anno <- as.data.table(signac@assays$anno@ranges)

anno_sets <- anno[, unique(anno)]

c <- celltypes[1]
a <- anno_sets[1]


pseudobulk <- map(celltypes, function(c){
  map(anno_sets, function(a){
    
    feats <- anno[anno == a, paste0(seqnames, "-", start, "-", end)]
    
    sub <- subset(signac, idents = c, features = feats)
    mat <- sub@assays$anno@counts
    dt <- data.table(id = rownames(mat), counts = rowSums(mat)) %>% 
      .[, c("anno", "celltype") := .(a, c)]
  }) %>% 
    rbindlist()
}) %>% 
  rbindlist()

# normalise
pseudobulk[, celltype_sum := sum(counts), celltype]
pseudobulk[, log2_cpm := log2(counts/celltype_sum)]

fwrite(pseudobulk, io$outfile, sep = "\t", na = "NA", quote = FALSE)


