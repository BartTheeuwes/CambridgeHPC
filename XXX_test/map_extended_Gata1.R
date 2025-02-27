## ----load, message = FALSE-------------------------------------------------------------------------------------
library(Matrix)
library(scran)
library(scater)
library(reshape2)
library(data.table)


# use https://github.com/MarioniLab/EmbryoTimecourse2018/blob/master/analysis_scripts/
# atlas/core_functions.R

setwd("/nfs/research/marioni/magda/chimera")
source("chimera_core_functions_extract.R")

load_atlas_data()

atlas_sce = sce
meta <- read.table("/nfs/research/marioni/magda/chimera/fromJonny/atlas_data/meta_submission.tab",
                   fill=TRUE,header=TRUE)
atlas_meta = meta

atlas_meta$celltype[atlas_meta$doublet] = "Doublet"
atlas_meta$celltype[atlas_meta$stripped] = "Stripped"
complete_atlas = atlas_sce
complete_meta = atlas_meta


#Downsample to 10k cells per sample (or maximum number of cells in sample)
set.seed(42)
keep = lapply(unique(atlas_meta$stage), function(x){
  if(x == "mixed_gastrulation"){
    return(c())
  } else if(sum(atlas_meta$stage == x) < 10000) {
    return(which(atlas_meta$stage == x))
  } else {
    hits = which(atlas_meta$stage == x)
    return(sample(hits, 10000))
  }
})
keep = do.call(c, keep)

atlas_sce = logNormCounts(atlas_sce[,keep])
atlas_meta = atlas_meta[keep,]


## ---- warning = FALSE, message = FALSE-------------------------------------------------------------------------
chimeraGata1 <- readRDS("data/chimeraGata1.rds")
names(colData(chimeraGata1))[names(colData(chimeraGata1)) == "Sample"] <- "sample"
chimeraGata1$cell <- colnames(counts(chimeraGata1))
sce_chimera <- SingleCellExperiment(assays = list(counts=counts(chimeraGata1)))
sce_chimera <- computeSumFactors(sce_chimera)
sce_chimera <- logNormCounts(sce_chimera,)
meta_chimera <- colData(chimeraGata1)
#xx <- sample(1:ncol(sce_chimera),100)
#sce_chimera <- sce_chimera[,xx]
#meta_chimera <- meta_chimera[xx,]


## ---- warning = FALSE, message = FALSE-------------------------------------------------------------------------
mappings = lapply(unique(meta_chimera$sample), function(x){
  mapWrap2(atlas_sce, atlas_meta, sce_chimera[-nrow(sce_chimera), meta_chimera$sample == x], meta_chimera[meta_chimera$sample == x,],
           target_name=paste0("results/Gata1/Gata1_sample",x))
})
mappings = do.call(rbind, mappings)

saveRDS(mappings, file = "/nfs/research/marioni/magda/chimera/data/Gata1_mapping/mapping_Gata1.rds")


mappings = readRDS("/nfs/research/marioni/magda/chimera/data/Gata1_mapping/mapping_Gata1.rds")

chimeraGata1$stage.mapped = as.character(mappings$stage.mapped[match(chimeraGata1$cell, mappings$cell)])
chimeraGata1$celltype.mapped = as.character(mappings$celltype.mapped[match(chimeraGata1$cell, mappings$cell)])
chimeraGata1$closest.cell = mappings$closest.cell[match(chimeraGata1$cell, mappings$cell)]
chimeraGata1$distance.to.closest.cell = mappings$distance.to.closest.cel[match(chimeraGata1$cell, mappings$cell)]
chimeraGata1$corr.to.closest.cell = mappings$corr.closest.cell[match(chimeraGata1$cell, mappings$cell)]
saveRDS(chimeraGata1,file="data/chimeraGata1.rds")


## ----sessinf---------------------------------------------------------------------------------------------------
sessionInfo()

