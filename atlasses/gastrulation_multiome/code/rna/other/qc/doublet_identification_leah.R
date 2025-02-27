library(scran)
library(scds)

#####################
## Define settings ##
#####################

source("/Users/ricard/gastrulation_multiome_10x/settings.R")
io$outdir <- paste0(io$basedir,"/results/rna/doublets")

# opts$k <- 20

opts$samples <- c(
    "E7.5_rep1",
    "E7.5_rep2",
    "E8.5_rep1",
    "E8.5_rep2"
)

############################
## Update sample metadata ##
############################

io$metadata <- paste0(io$basedir,"/results/rna/mapping/sample_metadata_after_mapping.txt.gz")
sample_metadata <- fread(io$metadata) %>% 
    .[pass_QC==TRUE & sample%in%opts$samples]

###############
## Load data ##
###############

# load SingleCellExperiment object
# io$sce <- "/Users/ricard/data/gastrulation_multiome_10x/multiome2/processed/SingleCellExperiment.rds"
sce <- load_SingleCellExperiment(io$sce, cells=sample_metadata$cell, normalise = TRUE)
dim(sce)

#############################
## Calculate doublet score ##
#############################

# Annotates doublets/multiplets using the hybrid approach in scds
sce <- cxds_bcds_hybrid(sce, estNdbl=TRUE)

dt <- colData(sce) %>%
    .[,c("cxds_score", "cxds_call", "bcds_score", "bcds_call", "hybrid_score", "hybrid_call")] %>%
    as.data.frame %>% tibble::rownames_to_column("cell") %>% as.data.table

dt[,hybrid_call:=hybrid_score>1]

table(dt$hybrid_call)
# Find neighbours using Seurat
# seurat <- NormalizeData(seurat)
# seurat <- FindVariableFeatures(seurat)
# seurat <- ScaleData(seurat)
# seurat <- RunPCA(seurat, verbose=FALSE)
# seurat <- FindNeighbors(seurat, k.param=opts$k, reduction = "pca")
# tmp <- data.table(
#     hybrid_map = rowSums(seurat@graphs$RNA_nn[,dt[hybrid_call==TRUE,cell]])/opts$k,
#     cell = rownames(seurat@graphs$RNA_nn)
# ) %>% merge(dt, by="cell", all.y=TRUE)


############################
## Update sample metadata ##
############################

sample_metadata <- fread(io$metadata) %>% 
    # .[,doublet_score:=NULL] %>%
    # .[,c("celltype.mapped.x","celltype.score.x","closest.cell.x"):=NULL] %>%
    .[,c("hybrid_score","hybrid_call"):=NULL] %>%
    merge(dt[,c("cell","hybrid_score","hybrid_call")], by="cell", all.x=TRUE)

head(sample_metadata)

# Save output
fwrite(sample_metadata, io$metadata, sep="\t", na="NA", quote=F)

#########################################
## Plot UMAP coloured by doublet score ##
#########################################

umap.dt <- fread(io$precomputed.umap) %>%
    merge(dt, by="cell")

p <- ggplot(umap.dt, aes(x=UMAP1, y=UMAP2, fill=hybrid_score)) +
    geom_point(size=1.5, shape=21, stroke=0.2) +
    scale_fill_gradientn(colours = terrain.colors(10)) +
    theme_classic() +
    theme(
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()
    )

# pdf(paste0(io$outdir,"/rna_umap.pdf"), width=5, height=3, useDingbats = F)
print(p)
# dev.off()

p <- ggplot(umap.dt, aes(x=UMAP1, y=UMAP2, fill=hybrid_call)) +
    geom_point(size=1.5, shape=21, stroke=0.2) +
    # scale_fill_gradientn(colours = terrain.colors(10)) +
    theme_classic() +
    theme(
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()
    )

# pdf(paste0(io$outdir,"/rna_umap.pdf"), width=5, height=3, useDingbats = F)
print(p)
# dev.off()

