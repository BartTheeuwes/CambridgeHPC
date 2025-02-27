library(SnapATAC)
library(purrr)
library(data.table)
library(Seurat)


source(here::here("settings.R"))

io$snap_file   <- file.path(io$basedir, "/processed/atac/snapatac/snapatac.rds")
io$outfile     <- file.path(io$basedir, "/processed/atac/snapatac/seurat_from_snap.rds")

snap <- readRDS(io$snap_file)
createGmatFromMat(obj, input.mat, genes, do.par, num.cores)

seurat <- snapToSeurat(snap)

saveRDS(seurat, io$outfile)
