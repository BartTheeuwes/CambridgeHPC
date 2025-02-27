suppressPackageStartupMessages({
  library("reticulate")
  library("SingleCellExperiment")
})

#####################
## Define settings ##
#####################

source("/homes/ricard/gastrulation_multiome_10x/settings.R")
io$outfile <- paste0(io$basedir,"/processed/rna/anndata.h5ad")

#####################################
## Reticulate connection to scanpy ##
#####################################

sc <- import("scanpy")

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(io$metadata) %>% 
  .[pass_rnaQC==TRUE & doublet_call==FALSE & !is.na(celltype.mapped)] %>%
  .[,c("cell", "sample", "stage", "nFeature_RNA", "mitochondrial_percent_RNA", "celltype.mapped")]

###############
## Load data ##
###############

# Load RNA expression data as SingleCellExperiment object
sce <- load_SingleCellExperiment(io$rna.sce, cells=sample_metadata$cell, normalise = FALSE)

# Add sample metadata as colData
colData(sce) <- sample_metadata %>% tibble::column_to_rownames("cell") %>% DataFrame


#############################################
## Convert SingleCellExperiment to AnnData ##
#############################################

adata_sce <- sc$AnnData(
    X   = t(counts(sce)),
    obs = as.data.frame(colData(sce)),
    var = data.frame(gene=rownames(sce), row.names=rownames(sce))
)
# adata_sce$obsm$update(umap = reducedDim(sce, "umap"))

adata_sce


# Add cell type colors
adata_sce$uns$update(celltype.mapped_colors = opts$celltype.colors[sort(unique(as.character(adata_sce$obs$celltype.mapped)))])
adata_sce$uns["celltype.mapped_colors"]

##########
## Save ##
##########

adata_sce$write_h5ad(io$outfile)
