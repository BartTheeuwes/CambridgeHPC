# library(Seurat)
library(SingleCellExperiment)
library(scater)
library(scran)

# Load default settings
source("/Users/ricard/gastrulation_multiome_10x/settings.R")
# source("/homes/ricard/gastrulation_multiome_10x/settings.R")

# Define I/O
# io$outfile <- io$sce
io$outfile <- "/Users/ricard/data/gastrulation_multiome_10x/multiome2/processed/SingleCellExperiment.rds"

# Load sample metadata
sample_metadata <- fread(io$metadata)

# Load seurat 
seurat <- readRDS(io$seurat)

# Convert to SingleCellExperiment
sce <- Seurat::as.SingleCellExperiment(seurat)

# Add metadata
# stopifnot(sample_metadata$cell%in%colnames(sce))
# stopifnot(colnames(sce)%in%sample_metadata$cell)
# sample_metadata <- sample_metadata %>% setkey(cell) %>% .[colnames(sce)]
# stopifnot(sample_metadata$cell == colnames(sce))
# colData(sce) <- sample_metadata %>% as.data.frame %>% tibble::column_to_rownames("cell") %>%
#   .[colnames(sce),] %>% DataFrame()

# Normalise using scran
clusts = as.numeric(quickCluster(sce, method = "igraph", min.size = 100, BPPARAM = mcparam))
min.clust = min(table(clusts))/2
new_sizes = c(floor(min.clust/3), floor(min.clust/2), floor(min.clust))
sce = computeSumFactors(sce, clusters = clusts, sizes = new_sizes, max.cluster.size = 3000)

to.plot <- data.frame(X = Matrix::colSums(counts(sce)), Y = sizeFactors(sce))
ggplot(to.plot, mapping = aes(x = X, y = Y)) +
  geom_point() +
  labs(x = "Number of UMIs", y = "Size Factor") +
  theme_classic()

sce <- logNormCounts(sce)

saveRDS(sce, io$outfile)
