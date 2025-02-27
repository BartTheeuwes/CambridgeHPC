# https://www.archrproject.com/bookdown/chromvar-deviatons-enrichment-with-archr.html
# Note: this requires the previous execution of (...)/add_motif_annotation/archR_add_motif_annotation.R
suppressMessages(library(GenomicRanges))
suppressMessages(library(argparse))

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--motif_annotation',  type="character",                help='Motif annotation') 
p$add_argument('--metadata',  type="character",                help='Metadata') 
p$add_argument('--threads',            type="integer",    default=1,    help='Number of cores')
p$add_argument('--force', action="store_true", 				help='Force')
p$add_argument('--outdir',  type="character",                help='Output directory') 
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

source(here::here("settings.R"))
source(here::here("utils.R"))

## START TEST ##
# args <- list()
# args$motif_annotation <- "CISBP"
# args$metadata <- file.path(io$basedir,"results_new/atac/archR/qc/sample_metadata_after_qc.txt.gz")
# args$threads <- 2
# args$force <- FALSE
# args$outdir <- file.path(io$basedir,"results_new/atac/archR/chromvar")
## END TEST ##

#####################
## Load metadata ##
#####################

sample_metadata <- fread(args$metadata) %>%
  .[pass_atacQC==TRUE & doublet_call==FALSE]

########################
## Load ArchR Project ##
########################

source(here::here("atac/archR/load_archR_project.R"))

addArchRThreads(threads = args$threads) 
stopifnot(args$motif_annotation%in%names(ArchRProject@peakAnnotation))

# Subset ArchR
stopifnot(sample_metadata$cell%in%rownames(ArchRProject))
ArchRProject.filt <- ArchRProject[sample_metadata$cell,]

##################
## Add Peak Set ##
##################

# Load Peak Set
# peaks.dt <- fread(io$archR.peakSet.bed) %>%
#   setnames(c("chr","start","end")) %>%
#   .[,id:=sprintf("%s:%s-%s",chr,start,end)]
# peaks.granges <- makeGRangesFromDataFrame(peaks.dt, keep.extra.columns = T)
# ArchRProject.filt <- addPeakSet(ArchRProject.filt, peaks.granges)

##########################
## Add background peaks ##
##########################

# This function will compute background peaks controlling for total accessibility and GC-content
# changes in the ArchR project: (1) it creates Background-Peaks.rds and (2) adds "bgdPeaks" entry to "metadata(getPeakSet(ArchRProject))"

# Background peaks are chosen by sampling peaks based on similarity in GC content and # of fragments across samples using the Mahalanobis distance. 
# The "w" paramter controls how similar background peaks should be. The bs parameter controls the precision with which the similarity is computed; 
# increasing bs will make the function run slower.
# Returns a matrix with one row per peak and one column per iteration. values in a row represent indices of background peaks for the peak with that index

# if (!file.exists(metadata(ArchRProject.filt@peakSet)$bgdPeaks)) {
#   ArchRProject.filt <- addBgdPeaks(
#     ArchRProject.filt,
#     nIterations = 50,
#     w = 0.1,
#     binSize = 50,
#     method = "chromVAR",
#     seed = 42,
#     # outFile = file.path(getOutputDirectory(ArchRProj), "Background-Peaks.rds"),
#     force = TRUE
#   )
# }

###################################
## Compute deviations for motifs ##
###################################

# The function computeDeviations returns a SummarizedExperiment with two "assays":  
# - The first matrix (accessible via `deviations(dev)` or `assays(dev)$deviations)` will give the bias corrected deviation in accessibility for each set of peaks (rows) for each cell or sample (columns). This metric represent how accessible the set of peaks is relative to the expectation based on equal chromatin accessibility profiles across cells/samples, normalized by a set of background peak sets matched for GC and average accessibility. 
# - The second matrix `deviationScores(dev)` or `assays(deviations)$z` gives the deviation Z-score, which takes into account how likely such a score would occur if randomly sampling sets of beaks with similar GC content and average accessibility.

matrix_name <- paste0("DeviationMatrix_",args$motif_annotation)

if (matrix_name %in% getAvailableMatrices(ArchRProject)) {
	if (args$force) {
		ArchRProject.filt <- addDeviationsMatrix(
		  ArchRProject.filt, 
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
	ArchRProject.filt <- addDeviationsMatrix(
	  ArchRProject.filt, 
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
chromvar_deviations.se <- getMatrixFromProject(ArchRProject.filt, matrix_name)

# Save SummarisedExperiment object
saveRDS(chromvar_deviations.se, file.path(args$outdir,sprintf("chromVAR_deviations_%s_archr.rds",args$motif_annotation)))

##########
## TEST ##
##########

# if (grepl("ricard",Sys.info()['nodename'])) {
#   R.utils::sourceDirectory("/Users/ricard/git/ArchR/R/", verbose=T, modifiedOnly=FALSE)
# } else if (grepl("ebi",Sys.info()['nodename'])) {
#   R.utils::sourceDirectory("/homes/ricard/git/ArchR/R/", verbose=T, modifiedOnly=FALSE)
# } else {
#   stop("Computer not recognised")
# }

# ArchRProj = ArchRProject.filt
# peakAnnotation = "Motif_JASPAR2020_human"
# matches = NULL
# bgdPeaks = getBgdPeaks(ArchRProj, method = "chromVAR")
# matrixName = NULL
# out = c("z", "deviations")
# binarize = FALSE
# threads = 1
# verbose = TRUE
# parallelParam = NULL
# force = FALSE
# logFile = createLogFile("addDeviationsMatrix")
