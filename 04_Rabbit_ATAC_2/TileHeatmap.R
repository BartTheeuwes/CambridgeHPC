suppressPackageStartupMessages({
    library(ArchR) 
    library(data.table)
    library(purrr)
    library(parallel)
    library(dplyr)
    library(Matrix)
    library(SingleCellExperiment)
    library(Seurat)
    library(reshape2)
    library(scater)
    library(viridis)
    library(gridExtra)
    library(ggpubr)
})

options(repr.plot.width=15, repr.plot.height=8)


# I/O
io = list()
io$basedir='/rds/project/rds-SDzz0CATGms/users/bt392/04_Rabbit_ATAC_2'
io$output.directory <- file.path(io$basedir,"ArchR")
io$plot.dir = file.path(io$output.directory,'Plots')
io$plot.genes = file.path(io$plot.dir,'Genes')
io$markers = file.path(io$output.directory,'Markers')

dir.create(file.path(io$plot.genes), showWarnings = FALSE)
dir.create(file.path(io$markers), showWarnings = FALSE)


setwd(io$output.directory)

io$archR.directory = file.path(io$output.directory, 'Project/')
ArchRProject.filt = loadArchRProject(io$archR.directory)

# Tiles
markersTile <- getMarkerFeatures(
    ArchRProj = ArchRProject.filt, 
    useMatrix = "TileMatrix", 
    groupBy = "Clusters",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

saveRDS(markersTile, file.path(io$markers, 'markersTile.rds'))

heatmapGS <- plotMarkerHeatmap(
  seMarker = markersTile, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25", 
  transpose = TRUE
)

plotPDF(heatmapGS, name = "Tile-Marker-Heatmap", width = 12, height = 9, ArchRProj = ArchRProject.filt, addDOC = FALSE)