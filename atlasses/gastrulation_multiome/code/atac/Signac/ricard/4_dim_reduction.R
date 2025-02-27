#################
## Description ##
#################

# loads a metadata file and incorporates it into the Seurat object

####################
## Load libraries ##
####################

library(Seurat)
library(Signac)

#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
} else {
  stop("Computer not recognised")
}

# I/O
io$metadata <- paste0(io$basedir, "/sample_metadata.txt.gz")
io$outdir <- paste0(io$basedir, "/results/atac/signac/dimensionality_reduction")

# Options
opts$cores       <- 8
opts$mem         <- 5 # GB
opts$assay       <- "peaks"

# Multiprocessing
plan("multiprocess", workers = opts$cores)
options(future.globals.maxSize = opts$mem * 1024 ^ 3)
plan()

#################
## Load Signac ##
#################

signac <- readRDS(io$signac)
signac

# switch assay
DefaultAssay(signac) <- opts$assay
signac

################
## Parse data ##
################

signac <- FindTopFeatures(signac, min.cutoff = "q10")
signac <- RunTFIDF(signac)

###############
## LSI (SVD) ##
###############

signac <- RunSVD(signac)

# correlation between sequencing depth and components
head(cor(signac@reductions$lsi@cell.embeddings, signac@meta.data$nFeature_bins))
# depthcor <- DepthCor(signac)
# save_plot(paste0(io$plots_out, "/depth_cor.pdf"), depthcor)

##########
## UMAP ##
##########

signac <- RunUMAP(object = signac, reduction = 'lsi', dims = 2:30, n.neighbors = 30L, min.dist = 0.3)

##########
## Plot ##
##########

Idents(signac) <- signac$celltype.predicted
p <- DimPlot(signac, group.by = 'celltype.predicted', label = TRUE) + NoLegend() + NoAxes()

pdf(sprintf("%s/signac_umap_celltype.pdf",io$outdir))
print(p)
dev.off()

##########
## Save ##
##########

# Save Signac object
# saveRDS(signac, io$signac)

# Save coordinates
umap.dt <- signac@reductions$umap@cell.embeddings %>% as.data.table %>% 
  .[,cell:=colnames(signac)] %>%
  setnames(c("x","y","cell"))
fwrite(umap.dt, paste0(io$outdir,"/umap.txt.gz"))
