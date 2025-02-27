
####################
## Load libraries ##
####################

suppressPackageStartupMessages(library(Signac))

#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/public_datasets/Pijuan-Sala_2020/settings.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/public_datasets/Pijuan-Sala_2020/settings.R")
} else {
  stop("Computer not recognised")
}

# Define I/O
io$outdir <- paste0(io$basedir,"/results/results/differential")

# Define options
opts$test <- TRUE

#####################
## Update metadata ##
#####################

if (opts$test) sample_metadata <- head(sample_metadata,n=100)

###############
## Load data ##
###############

# Load Signac object
seurat <- readRDS(io$seurat)[,sample_metadata$cell]
dim(seurat)
