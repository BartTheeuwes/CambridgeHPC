suppressPackageStartupMessages(library(ArchR))

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
setwd(io$archR.directory)


################
## Define I/O ##
################

io$archR.directory <- paste0(io$basedir,"/processed/atac/archR_integrated")
io$metadata <- paste0(io$archR.directory,"/sample_metadata.txt.gz")
# io$archR.projectMetadata <- paste0(io$archR.directory,"/projectMetadata.rds")
# io$archR.peakSet.granges <- paste0(io$archR.directory,"/peakSet.rds")
# io$archR.bgdPeaks <- paste0(io$archR.directory,"/Background-Peaks.rds")
# io$archR.peakSet.bed <- paste0(io$archR.directory,"/PeakCalls/bed/peaks_archR_macs2.bed.gz")
# io$archR.pseudobulk.peaks <- paste0(io$archR.directory,"/pseudobulk/pseudobulk_PeakMatrix_summarized_experiment.rds")
# io$archR.peak.variability <- paste0(io$basedir,"/results/atac/archR/variability/peak_variability.txt.gz")
# io$archR.deviations.se <- paste0(io$basedir,"/results/atac/archR/chromvar/cisbp/deviations_summarized_experiment.rds")

####################
## Define options ##
####################


opts$stages <- c(
  "E7.5",
  "E8.25",
  "E8.5"
)

opts$samples <- c(
  "E7.5_rep1",
  "E7.5_rep2",
  "E8.25_PijuanSala",
  "E8.5_rep1",
  "E8.5_rep2"
)


# ArchR options
addArchRThreads(threads = 1) 
addArchRGenome("mm10")

########################
## Load ArchR project ##
########################

ArchRProject <- loadArchRProject(io$archR.directory)

# Load ArchR projectMetadata
# ArchRProject@projectMetadata <- readRDS(io$archR.projectMetadata)

# Load peaks
# ArchRProject <- addPeakSet(ArchRProject, peakSet = readRDS(io$archR.peakSet.granges), force = TRUE)

# Add background peaks
# if (!"bgdPeaks" %in% metadata(getPeakSet(ArchRProject))$bgdPeaks) {
#   metadata(ArchRProject@peakSet)$bgdPeaks <- io$archR.bgdPeaks
# }

# Load motif annotations over peaks
# ArchRProject@peakAnnotation <- readRDS(sprintf("%s/Annotations/peakAnnotation.rds",io$archR.directory))
