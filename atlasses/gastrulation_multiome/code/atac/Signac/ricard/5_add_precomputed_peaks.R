#################
## Description ##
#################

# Add pre-computed peaks into the Signac object

####################
## Load libraries ##
####################

library(Seurat)
library(Signac)
library(GenomicRanges)
library(future)

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
io$metadata      <- paste0(io$basedir, "/sample_metadata.txt.gz")
# io$signac        <- file.path(io$basedir, "/processed/atac/signac/signac.rds")
io$outdir        <- paste0(io$basedir, "/results/atac/signac/dimensionality_reduction")
io$signac_out    <- io$signac
io$peaks         <- paste0(io$archR.directory,"/signac_peaks_per_cluster.rds")

# Options
opts$cores       <- 3
opts$mem         <- 5 # GB
opts$assay       <- "peaks"
opts$group_by    <- "celltype.predicted"
opts$remove.bins.assay <- TRUE

# Multiprocessing
plan("multiprocess", workers = opts$cores)
options(future.globals.maxSize = opts$mem * 1024 ^ 3)
plan()

#################
## Load Signac ##
#################

signac <- readRDS(io$signac)

# subset current assay to avoid out-of-memory problems
signac <- signac[1:100,]

signac

###################
## Load peak set ##
###################

io$peaks <- io$archR.peakSet.granges
peaks.gr <- readRDS(io$peaks)

##################################
## quantify counts in each peak ##
##################################

# returns a sparse matrix
counts.mtx <- FeatureMatrix(
  fragments = Fragments(signac),
  features = peaks.gr,
  cells = colnames(signac)
)

# create a new assay using the peak set and add it to the Seurat object
signac[[opts$assay]] <- CreateChromatinAssay(
  counts = counts.mtx,
  fragments = Fragments(signac),
  annotation = Annotation(signac)
)

# change assay
DefaultAssay(signac) <- opts$assay

if (opts$remove.bins.assay & "bins"%in%Assays(signac)) {
  signac[['bins']] <- NULL
}

##########
## Save ##
##########

saveRDS(signac, io$signac_out)
