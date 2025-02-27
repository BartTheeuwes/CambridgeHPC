library(chromVARmotifs)
library(BSgenome.Mmusculus.UCSC.mm10)
library(motifmatchr)

#####################
## Define settings ##
#####################

# Load default settings
if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
  source("/Users/ricard/gastrulation_multiome_10x/atac/archR/add_motif_annotation/utils.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
  source("/homes/ricard/gastrulation_multiome_10x/atac/archR/add_motif_annotation/utils.R")
} else {
  stop("Computer not recognised")
}

opts$test <- FALSE

########################
## Load ArchR Project ##
########################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/atac/archR/load_archR_project.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/atac/archR/load_archR_project.R")
} else {
  stop("Computer not recognised")
}

######################
## Define variables ##
######################

name = "Motif_cisbp"
cutOff = 5e-5
width = 7

# name = "Motif_cisbp_lenient"
# cutOff = 1e-4
# width = 7

#####################
## Run motifmatchr ##
#####################

# Get PWM List adapted from chromVAR!
data("human_pwms_v2")
# motifs <- mouse_pwms_v2
motifs <- human_pwms_v2
obj <- .summarizeChromVARMotifs(motifs)
motifs <- obj$motifs
motifSummary <- obj$motifSummary
  
if (opts$test) {
  motifs <- motifs[1:10]
  motifSummary <- motifSummary[1:10,]
}

stopifnot(names(motifs)==rownames(motifSummary))


######################
## Get motif scores ##
######################

motifScores <- matchMotifs(
  pwms = motifs,
  subject = ArchRProject@peakSet,
  genome = BSgenome.Mmusculus.UCSC.mm10, 
  out = "scores", 
  p.cutoff = cutOff, 
  w = width
)

##########
## Save ##
##########

class(assay(motifScores,"motifScores"))
class(assay(motifScores,"motifMatches"))
class(assay(motifScores,"motifCounts"))
  
# Save output
# io$outfile <- sprintf("%s/Annotations/%s_cutOff%s_width%s-Scores.rds",io$archR.directory,name,cutOff,width)
io$outfile <- sprintf("%s/Annotations/%s-Scores.rds",io$archR.directory,name)
saveRDS(motifScores, io$outfile, compress = TRUE)
