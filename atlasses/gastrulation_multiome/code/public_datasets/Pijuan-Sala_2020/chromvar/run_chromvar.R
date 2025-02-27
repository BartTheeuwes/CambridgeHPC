
####################
## Load libraries ##
####################

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Signac))
suppressPackageStartupMessages(library(chromVAR))
suppressPackageStartupMessages(library(BSgenome.Mmusculus.UCSC.mm10))
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(chromVARmotifs))
suppressPackageStartupMessages(library(motifmatchr))

#####################
## Define settings ##
#####################


# Load default settings
if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/public_datasets/Pijuan-Sala_2020/settings.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/public_datasets/Pijuan-Sala_2020/settings.R")
} else {
  stop("Computer not recognised")
}

# Define I/O
io$outdir <- paste0(io$basedir,"/results/chromVAR")

# Define options
opts$test <- FALSE


#####################
## Update metadata ##
#####################

if (opts$test) sample_metadata <- head(sample_metadata,n=100)


###############
## Load data ##
###############

# Load Signac object
seurat <- readRDS(io$seurat)[,sample_metadata$cell]
if (opts$test) {
  seurat <- seurat[1:1000,] 
  seurat@assays$peaks@ranges <- seurat@assays$peaks@ranges[1:1000]
}


##############
## ChromVar ##
##############

## Create SummarizedExperiment object

peak.matrix <- GetAssayData(seurat, slot = "counts")
peak.matrix <- peak.matrix[rowSums(peak.matrix)>0, , drop = FALSE]
seurat <- seurat[rownames(peak.matrix),]
seurat@assays$peaks@ranges <- seurat@assays$peaks@ranges[which(rowSums(peak.matrix)>0)]
peak.ranges <- Signac::granges(seurat)
chromvar.obj <- SummarizedExperiment(
  assays = list(counts = peak.matrix),
  rowRanges = peak.ranges
)


## Compute GC content for peaks

chromvar.obj <- addGCBias(chromvar.obj, genome = BSgenome.Mmusculus.UCSC.mm10)

## Get Background peaks

# Background peaks are chosen by sampling peaks based on similarity in GC content and # of fragments across samples using the Mahalanobis distance. 
# The w paramter controls how similar background peaks should be. The bs parameter controls the precision with which the similarity is computed; 
# increasing bs will make the function run slower.
# Returns a matrix with one row per peak and one column per iteration. values in a row represent indices of background peaks for the peak with that index

bg <- getBackgroundPeaks(chromvar.obj, niterations = 50, w = 0.1, bs = 50)

## Collect motifs

# The function matchMotifs from the motifmatchr package finds which peaks contain which motifs. Returns:
# - out=="matches": returns a SummarizedExperiment with a sparse matrix with values set to TRUE for a match (default)
# - out=="scores": returns a SummarizedExperiment with a matches matrix as well as matrices with the maximum motif score and total motif counts
# - out=="positions": a GenomicRangesList with all the positions of matches

# Background nucleotide frequencies can be set to "subject" to use the subject sequences or ranges for computing the nucleotide frequencies, "genome" for using the genomice frequencies (in which case a genome must be specified), "even" for using 0.25 for each base, or a numeric vector with A, C, G, and T frequencies.
data("mouse_pwms_v2")
motifs = mouse_pwms_v2
motif_ix <- matchMotifs(motifs, seurat@assays$peaks@ranges, genome = BSgenome.Mmusculus.UCSC.mm10)


# If instead of using known motifs, you want to use all kmers of a certain length, the matchKmers function can be used. For more about using kmers as inputs, see the the Annotations section of the documentation website.
# kmer_ix <- matchKmers(
#   k = 6, 
#   subject = seurat@assays$peaks@ranges, 
#   genome = BSgenome.Mmusculus.UCSC.mm10,
#   out = "matches"
# )


## Compute deviations

# The function computeDeviations returns a SummarizedExperiment with two "assays":  
# - The first matrix (accessible via `deviations(dev)` or `assays(dev)$deviations)` will give the bias corrected deviation in accessibility for each set of peaks (rows) for each cell or sample (columns). This metric represent how accessible the set of peaks is relative to the expectation based on equal chromatin accessibility profiles across cells/samples, normalized by a set of background peak sets matched for GC and average accessibility. 
# - The second matrix `deviationScores(dev)` or `assays(deviations)$z` gives the deviation Z-score, which takes into account how likely such a score would occur if randomly sampling sets of beaks with similar GC content and average accessibility.

deviations <- computeDeviations(chromvar.obj,
  annotations = motif_ix,
  background_peaks = bg
)

##########
## Save ##
##########

saveRDS(chromvar.obj, paste0(io$outdir,"/chromvar_object.rds"))
saveRDS(deviations, paste0(io$outdir,"/chromvar_deviations.rds"))
