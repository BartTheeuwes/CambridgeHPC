library(Seurat)

#####################
## Define settings ##
#####################

# define I/O
source("/Users/ricard/gastrulation_multiome_10x/settings.R")
io$outdir <- paste0(io$basedir,"/results/rna/cell_cycle")

# Define options

############################
## Update sample metadata ##
############################

sample_metadata <- sample_metadata %>%
  .[pass_rnaQC==TRUE]

###############
## Load data ##
###############

# Load RNA expression data as Seurat object
seurat <- readRDS(io$seurat)[,sample_metadata$cell]
dim(seurat)

# Add metadata to the Seurat object
seurat@meta.data <- sample_metadata %>% as.data.frame %>% tibble::column_to_rownames("cell") %>%
  .[colnames(seurat),]
head(seurat@meta.data)

###############
## Normalise ##
###############

# SCTransform
# seurat <- SCTransform(seurat)

# LogNormalize
seurat <- NormalizeData(seurat, normalization.method = "LogNormalize")

# Scale
seurat <- ScaleData(seurat, do.scale = F)

##########################
## Cell cycle inference ##
##########################


# Rename genes
row.names <- toupper(rownames(seurat))
rownames(seurat@assays$RNA@counts) <- row.names
rownames(seurat@assays$RNA@data) <- row.names
rownames(seurat@assays$RNA@scale.data) <- row.names
# seurat[["RNA"]]@meta.features <- data.frame(row.names=foo)

s.genes <- fread("/Users/ricard/data/mm10_regulation/cellcycle/s_genes.txt", header = FALSE)$V1
g2m.genes <- fread("/Users/ricard/data/mm10_regulation/cellcycle/g2m_genes.txt", header = FALSE)$V1

s.genes <- s.genes[s.genes%in%rownames(seurat)]
g2m.genes <- g2m.genes[g2m.genes%in%rownames(seurat)]

seurat <- CellCycleScoring(seurat, s.features = s.genes, g2m.features = g2m.genes)

##########
## Save ##
##########

dt <- seurat@meta.data[,c("S.Score","G2M.Score","Phase")] %>%
  as.data.table(keep.rownames = T) %>%
  setnames("rn","cell")

# fwrite(dt, paste0(io$outdir,"/cell_cycle_seurat.txt.gz"))  

################################
## Merge with sample metadata ##
################################

sample_metadata <- fread(io$metadata) %>% 
  merge(dt,by="cell",all.x=TRUE)
head(sample_metadata)

fwrite(sample_metadata, io$metadata, sep="\t", na="NA", quote=F)
