source("/Users/ricard/gastrulation_multiome_10x/settings.R")
seurat1 <- readRDS("/Users/ricard/data/gastrulation_multiome_10x/multiome1/processed/seurat.rds")
seurat2 <- readRDS("/Users/ricard/data/gastrulation_multiome_10x/multiome2/processed/seurat.rds")

metadata1 <- fread("/Users/ricard/data/gastrulation_multiome_10x/multiome1/sample_metadata.csv") %>%
  .[,cell:=stringr::str_replace_all(cell,"-1","")] %>%
  .[,cell:=paste0("multiome1_",cell)] %>%
  .[,batch:="E8.5_rep1"]
metadata2 <- fread("/Users/ricard/data/gastrulation_multiome_10x/multiome2/sample_metadata.csv") %>%
  .[,stage:="E8.5"] %>%
  .[,batch:="E8.5_rep2"]

# cols <- c("cell", "batch", "nCount_RNA", "nFeature_RNA", "nCount_ATAC", "nFeature_ATAC", "percent.mt", "percent.ribo", "pass_rnaQC", "celltype.mapped", 
#           "celltype.score", "closest.cell", "stage", "S.Score", "G2M.Score", "Phase", "hybrid_score", "hybrid_call")

cols <- c("cell", "batch", "nCount_RNA", "nFeature_RNA", "nCount_ATAC", "nFeature_ATAC", "percent.mt", "percent.ribo", "pass_rnaQC","stage")

metadata1.new <- metadata1[,..cols]
sum(duplicated(metadata1.new$cell))

metadata2.new <- metadata2[,..cols]
sum(duplicated(metadata2.new$cell))

metadata.new <- rbind(metadata1.new,metadata2.new)
sum(duplicated(metadata.new$cell))

fwrite(metadata.new, io$metadata, sep="\t", quote=F)

