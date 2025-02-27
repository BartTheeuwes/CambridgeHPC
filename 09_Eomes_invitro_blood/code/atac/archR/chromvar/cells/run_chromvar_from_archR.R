# https://www.archrproject.com/bookdown/chromvar-deviatons-enrichment-with-archr.html
# Note: this requires the previous execution of (...)/add_motif_annotation/archR_add_motif_annotation.R
here::i_am("atac/archR/chromvar/cells/run_chromvar_from_archR.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

suppressMessages(library(GenomicRanges))
suppressMessages(library(ArchR))

rhdf5::h5disableFileLocking()

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--archr_directory',    type="character",    help='ArchR directory')
p$add_argument('--motif_annotation',  type="character",                help='Motif annotation') 
p$add_argument('--metadata',  type="character",                help='Metadata') 
p$add_argument('--threads',            type="integer",    default=1,    help='Number of cores')
p$add_argument('--force', action="store_true", 				help='Force')
p$add_argument('--outdir',  type="character",                help='Output directory') 
args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# test=FALSE
# if(test==TRUE){
    # io$basedir <- file.path(io$basedir,"test")
    # args <- list()
    # args$archr_directory <- file.path(io$basedir,"processed/atac/archR")
    # args$motif_annotation <- "CISBP"
    # args$metadata <- file.path(io$basedir,"results/atac/archR/qc/sample_metadata_after_qc.txt.gz")
    # args$threads <- 1
    # args$force <- FALSE
    # args$outdir <- file.path(io$basedir,"results/atac/archR/chromvar")
# }
## END TEST ##

dir.create(args$outdir, recursive=TRUE, showWarnings = FALSE)

#####################
## Load metadata ##
#####################

sample_metadata <- fread(args$metadata) %>%
  .[pass_atacQC==TRUE & doublet_call==FALSE | pass_atacQC==TRUE & is.na(doublet_call)]

########################
## Load ArchR Project ##
########################

#print('test')
setwd(args$archr_directory)
ArchRProject <- loadArchRProject(args$archr_directory)[sample_metadata$cell]

#source(here::here("atac/archR/load_archR_project.R"))
#ArchRProject <- ArchRProject[sample_metadata$cell]

addArchRGenome("mm10")
addArchRThreads(threads = args$threads)

#######################################
## Sanity checks on the ArchR object ##
#######################################

# Load motif annotations over peaks
tmp <- file.path(args$archr_directory,"Annotations/peakAnnotation.rds")
if (file.exists(tmp)) {
	ArchRProject@peakAnnotation <- readRDS(tmp)
}
stopifnot(TRUE %in% grepl(paste(args$motif_annotation,collapse="|"), names(ArchRProject@peakAnnotation), ignore.case = TRUE))

##########################
## Add background peaks ##
##########################

metadata(ArchRProject@peakSet)$bgdPeaks <- file.path(args$archr_directory, "Background-Peaks.rds")

stopifnot("bgdPeaks"%in%names(metadata(ArchRProject@peakSet)))

###################################
## Compute deviations for motifs ##
###################################

# The function computeDeviations returns a SummarizedExperiment with two "assays":  
# - The first matrix (accessible via `deviations(dev)` or `assays(dev)$deviations)` will give the bias corrected deviation in accessibility for each set of peaks (rows) for each cell or sample (columns). This metric represent how accessible the set of peaks is relative to the expectation based on equal chromatin accessibility profiles across cells/samples, normalized by a set of background peak sets matched for GC and average accessibility. 
# - The second matrix `deviationScores(dev)` or `assays(deviations)$z` gives the deviation Z-score, which takes into account how likely such a score would occur if randomly sampling sets of beaks with similar GC content and average accessibility.

print("Running chromVAR implementation in ArchR...")

matrix_name <- paste0("DeviationMatrix_",args$motif_annotation)

# file to save to
outfile = file.path(args$outdir,sprintf("chromVAR_deviations_%s_archr.rds",args$motif_annotation))

# Rename motif annotation to correct name
args$motif_annotation = names(ArchRProject@peakAnnotation)[grepl(args$motif_annotation, names(ArchRProject@peakAnnotation), ignore.case = TRUE)]

if (matrix_name %in% getAvailableMatrices(ArchRProject)) {
	if (args$force) {
		ArchRProject <- addDeviationsMatrix(
		  ArchRProject, 
		  matrixName = matrix_name,
		  peakAnnotation = args$motif_annotation,
		  out = c("z", "deviations"),
		  binarize = FALSE,
		  force = TRUE
		)
	} else {
		stop(sprintf("%s already found in the ArchR object. Pass the --force argument to replace the existing data matrix...",matrix_name))
	}

} else {
	ArchRProject <- addDeviationsMatrix(
	  ArchRProject, 
	  matrixName = matrix_name,
	  peakAnnotation = args$motif_annotation,
	  out = c("z", "deviations"),
	  binarize = FALSE
	)
}

##########
## Save ##
##########

# Fetch deviations matrix as a SummarisedExperiment
addArchRThreads(threads = 1)
chromvar_deviations.se <- getMatrixFromProject(ArchRProject, matrix_name)

# Save SummarisedExperiment object
saveRDS(chromvar_deviations.se, outfile)