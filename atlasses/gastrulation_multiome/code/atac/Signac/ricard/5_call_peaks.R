#################
## Description ##
#################

# Peak calling using MACS2

####################
## Load libraries ##
####################

library(Seurat)
library(Signac)
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
io$peaks_out     <- paste0(io$signac.directory,"/signac_peaks_per_cluster.rds")

# Options
opts$cores       <- 1
opts$mem         <- 5 # GB
opts$assay       <- "peaks"
opts$group_by    <- "celltype.predicted"

# Multiprocessing
plan("multiprocess", workers = opts$cores)
options(future.globals.maxSize = opts$mem * 1024 ^ 3)
plan()

#################
## Load Signac ##
#################

signac <- readRDS(io$signac)
signac

# switch assay
# DefaultAssay(signac) <- opts$assay
# signac

###############
## Run MACS2 ##
###############

macs2 <- system("which macs2", intern = TRUE)
macs2

# Returns the '.narrowPeak' MACS output as a 'GRanges' object
peaks <- CallPeaks(
  object = signac,
  group.by = opts$group_by,
  macs2.path = macs2
)

##################
## Filter peaks ##
##################

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_mm10, invert = TRUE)

saveRDS(peaks, io$peaks_out)


##################################
## quantify counts in each peak ##
##################################

macs2_counts <- FeatureMatrix(
  fragments = Fragments(signac),
  features = peaks,
  cells = colnames(signac)
)

# create a new assay using the MACS2 peak set and add it to the Seurat object
signac[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = Fragments(signac),
  annotation = Annotation(signac)
)

# change assay
DefaultAssay(signac) <- "peaks"

# signac[["bins"]] <- NULL
signac@meta.data$sample %>% unique()

##########
## Save ##
##########

saveRDS(signac, io$signac_out)



