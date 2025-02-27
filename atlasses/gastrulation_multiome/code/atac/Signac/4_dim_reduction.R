library(Seurat)
library(Signac)
library(purrr)
library(data.table)
library(future)
library(cowplot)

set.seed(1234)



source(here::here("settings.R"))



io$signac        <- file.path(io$rawdata, "/processed/atac/signac/signac.rds")
io$metadata_out  <- file.path(io$rawdata, "/processed/atac/signac/signac_metadata.tsv.gz")
io$plots_out     <- "/bi/home/clarks/plots/10X_multiome/signac"

opts$cores       <- 8
opts$mem         <- 5 # GB
opts$assay       <- "peaks"#"bins"#"peaks"


io$plots_out <- file.path(io$plots_out, opts$assay)
dir.create(io$plots_out, recursive = TRUE)


plan("multiprocess", workers = opts$cores)
options(future.globals.maxSize = opts$mem * 1024 ^ 3)
plan()

signac <- readRDS(io$signac)
signac

signac@meta.data$sample %>% unique()

# switch assay
DefaultAssay(signac) <- opts$assay
signac

# run dimensionality reduction
signac <- FindTopFeatures(signac, min.cutoff = 5)
signac <- RunTFIDF(signac)
signac <- RunSVD(signac)



# correlation between sequencing depth and components
depthcor <- DepthCor(signac)
save_plot(paste0(io$plots_out, "/depth_cor.pdf"), depthcor)

# run umap
signac <- RunUMAP(object = signac, reduction = 'lsi', dims = 2:30)
signac <- FindNeighbors(object = signac, reduction = 'lsi', dims = 2:30)
signac <- FindClusters(object = signac, verbose = FALSE, algorithm = 3)



# umap plots
umap <- DimPlot(object = signac, label = TRUE) + NoLegend()
save_plot(paste0(io$plots_out, "/umap.pdf"), umap)

# save object
saveRDS(signac, io$signac)

# save metadata
sig_meta <- signac@meta.data %>% 
  cbind(signac[["umap"]]@cell.embeddings) %>% 
  setDT(keep.rownames = "cell")

fwrite(sig_meta, io$metadata_out, sep = "\t", na="NA")

# add cell types
Idents(signac)
Idents(signac) <- signac$celltype.mapped

umap_cells <- DimPlot(object = signac, label = TRUE, label.size = 2, group.by = "celltype.mapped", repel = TRUE) + NoLegend()

save_plot(paste0(io$plots_out, "/umap_celltypes.pdf"), umap_cells)


# save object
signac@meta.data$sample %>% unique()
saveRDS(signac, io$signac)

