library(SingleCellExperiment)
library(scran)
library(zellkonverter)

#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
} else {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
}

setwd("/Users/ricard/gastrulation_multiome_10x/rna/scanpy")

# Define I/O
io$metadata <- paste0(io$basedir,"/results/rna/doublets/sample_metadata_after_doublets.txt.gz")
io$outfile <- io$anndata

# Define options
opts$samples <- c(
  "E7.5_rep1",
  "E7.5_rep2",
  "E8.5_rep1",
  "E8.5_rep2"
)

###############
## Load data ##
###############

# Load cell metadata
sample_metadata <- fread(io$metadata) %>%
  .[pass_QC==TRUE & sample%in%opts$samples]

# Load SingleCellExperiment
sce <- load_SingleCellExperiment(io$sce,  cells=sample_metadata$cell)

print(object.size(sce@assays@data$counts), units='auto')
print(object.size(sce@assays@data$logcounts), units='auto')
sce@assays@data$logcounts[1:5,1:5]
sce@assays@data$counts[1:5,1:5]

# Add sample metadata as colData
colData(sce) <- sample_metadata %>% tibble::column_to_rownames("cell") %>% DataFrame

################
## Conversion ##
################

anndata.object <- SCE2AnnData(sce, X_name = "counts", skip_assays = FALSE)

##########
## Save ##
##########

# writeH5AD(anndata.object, file = io$outfile)
anndata.object$write_h5ad(io$outfile)
