library(data.table)
library(purrr)
library(Signac)
library(Seurat)


source(here::here("settings.R"))



io$signac           <- file.path(io$rawdata, "/processed/atac/signac/archr_signac.rds")

io$plots_out     <- "/bi/home/clarks/plots/10X_multiome/signac/from_archr"
dir.create(io$plots_out)

signac <- readRDS(io$signac)

# add cell types
Idents(signac)
Idents(signac) <- signac$celltype.predicted

umap_cells <- DimPlot(
  object = signac, 
  label = TRUE, 
  label.size = 2, 
  group.by = "celltype.mapped", 
  repel = TRUE
  ) + 
  NoLegend()
umap_cells

cowplot::save_plot(paste0(io$plots_out, "/umap_celltypes.pdf"), umap_cells)
