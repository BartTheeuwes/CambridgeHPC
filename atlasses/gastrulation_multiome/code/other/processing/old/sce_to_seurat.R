library(Seurat)
library(SingleCellExperiment)
library(scater)

#####################
## Define settings ##
#####################

# Load default settings
source("/Users/ricard/gastrulation_multiome_10x/settings.R")

# Define I/O
io$outfile <- io$seurat
# io$outfile <- "/Users/ricard/data/gastrulation_multiome_10x/multiome2/processed/SingleCellExperiment.rds"

###############
## Load data ##
###############

# Load sample metadata
sample_metadata <- fread(io$metadata)

# Load SingleCellExperiment 
sce <- readRDS(io$sce)

################
## Conversion ##
################

# Convert to Seurat

# Add metadata
# stopifnot(sample_metadata$cell%in%colnames(sce))
# stopifnot(colnames(sce)%in%sample_metadata$cell)
# sample_metadata <- sample_metadata %>% setkey(cell) %>% .[colnames(sce)]
# stopifnot(sample_metadata$cell == colnames(sce))
# colData(sce) <- sample_metadata %>% as.data.frame %>% tibble::column_to_rownames("cell") %>%
#   .[colnames(sce),] %>% DataFrame()

##########
## Save ##
##########

saveRDS(seurat, io$outfile)
