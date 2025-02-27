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
io$outfile <- io$signac

###################
## Load metadata ##
###################

sample_metadata <- fread(io$metadata) %>%
  .[pass_atacQC==TRUE]

#################
## Load Signac ##
#################

Signac <- readRDS(io$signac)

##########################
## Filter Signac object ##
##########################

# Select cells 
mean(colnames(Signac) %in% sample_metadata$cell)
cells <- intersect(sample_metadata$cell,colnames(Signac))

sample_metadata.filt <- sample_metadata %>% .[cell%in%cells] %>% setkey(cell) %>% .[cells]
Signac.filt <- Signac[,sample_metadata.filt$cell]

# Merge metadata
metadata.to.signac <- Signac.filt@meta.data %>% tibble::rownames_to_column("cell") %>% as.data.table %>%
  .[,c("cell","orig.ident", "nCount_bins", "nFeature_bins")] %>%
  merge(sample_metadata.filt, by=c("cell"))
rownames(metadata.to.signac) <- metadata.to.signac$cell

# Add metadata 
Signac.filt@meta.data <- metadata.to.signac

# Sanity check
stopifnot(colnames(Signac.filt@assays$bins@data) == rownames(Signac.filt@meta.data))

##########
## Save ##
##########

saveRDS(Signac.filt, io$outfile)