# ## Mapping to extended atlas

# Load packages
suppressPackageStartupMessages({
    library(Seurat)
    library(scran)
    library(scater)
    library(batchelor)
    library(SingleCellExperiment)
    library(scater)
    library(Matrix)
    library(reshape2)
    library(data.table)
    library(BiocParallel)
    library(dplyr)
    library(gridExtra)
})

    ncores = 4
    mcparam = MulticoreParam(workers = ncores)
    register(mcparam)
    BPPARAM = SerialParam()

options(repr.plot.width=15, repr.plot.height=8)

main = "/rds/project/rds-SDzz0CATGms/users/bt392/01_Eomes_RNA/"
out_dir = paste0(main, '07_atlas_mapping/')
plot_dir = paste0(main, '07_atlas_mapping/plots/')
dir.create(out_dir, showWarnings = FALSE)
dir.create(plot_dir, showWarnings = FALSE)

source(paste0(main, "mapping_functions.R"))

load_data(normalise=TRUE)
sce_query = sce[Matrix::rowSums(counts(sce)) > 0,]
meta_query = meta

load_atlas_data(normalise=TRUE)
sce_atlas = sce[Matrix::rowSums(counts(sce)) > 0,]
meta_atlas = meta

# genes shared across datasets
shared_genes = merge(data.frame('V1'=names(sce_atlas)), data.frame('V1'=names(sce_query)), by='V1')

# filter cells & genes chimera dataset
sce_query = sce_query[shared_genes$V1,!c(meta_query$doublet | meta_query$stripped)]
meta_query = meta_query[!c(meta_query$doublet | meta_query$stripped),]

# filter cells & genes atlas dataset
sce_atlas = sce_atlas[shared_genes$V1,!c(meta_atlas$doublet | meta_atlas$stripped)]
meta_atlas = meta_atlas[!c(meta_atlas$doublet | meta_atlas$stripped),]

# +
#Downsample to 10k cells per sample (or maximum number of cells in sample)
set.seed(42)
keep = lapply(unique(meta_atlas$stage), function(x){
  if(x == "mixed_gastrulation"){
    return(c())
  } else if(sum(meta_atlas$stage == x) < 10000) {
    return(which(meta_atlas$stage == x))
  } else {
    hits = which(meta_atlas$stage == x)
    return(sample(hits, 10000))
  }
})
keep = do.call(c, keep)

sce_atlas = logNormCounts(sce_atlas[,keep])
meta_atlas = meta_atlas[keep,]
# -

mappings = lapply(unique(meta_query$sample), function(x){
    mapWrap(
        atlas_sce = sce_atlas, 
        atlas_meta = meta_atlas,
        map_sce = sce_query[, meta_query$sample == x], 
        map_meta = meta_query[meta_query$sample == x,], 
        # genes = marker_genes, 
        npcs = 50, 
        k = 30)$mapping
})

mappings = do.call(rbind, mappings)

write.csv(mappings, file =  paste0(out_dir,"mapped.csv"),row.names=FALSE)
