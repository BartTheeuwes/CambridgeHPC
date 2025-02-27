suppressPackageStartupMessages(library(ArchR))

################
## Define I/O ##
################

setwd(io$archR.directory)

####################
## Define options ##
####################

addArchRGenome("mm10")
addArchRThreads(threads = 1) 

########################
## Load ArchR project ##
########################

ArchRProject <- loadArchRProject(io$archR.directory)

# Load ArchR projectMetadata
# if (file.exists(io$archR.projectMetadata)) {
# 	ArchRProject@projectMetadata <- readRDS(io$archR.projectMetadata)
# }

# Load peaks
# if (file.exists(io$archR.peakSet.granges)) {
# 	ArchRProject <- addPeakSet(ArchRProject, peakSet = readRDS(io$archR.peakSet.granges), force = TRUE)
# }

# Load motif annotations over peaks
# io$archR.peakAnnotation <- sprintf("%s/Annotations/peakAnnotation.rds",io$archR.directory)
# if (file.exists(io$archR.peakSet.granges)) {
# 	ArchRProject@peakAnnotation <- readRDS(io$archR.peakAnnotation)
# }

# Add background peaks
# io$archR.bgdPeaks <- file.path(io$archR.directory, "Background-Peaks.rds")
# if (!"bgdPeaks" %in% metadata(getPeakSet(ArchRProject))$bgdPeaks) {
# 	if (file.exists(io$archR.bgdPeaks)) metadata(ArchRProject@peakSet)$bgdPeaks <- io$archR.bgdPeaks
# }
