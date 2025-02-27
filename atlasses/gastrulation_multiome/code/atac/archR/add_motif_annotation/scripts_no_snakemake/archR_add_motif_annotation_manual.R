library(chromVARmotifs)
library(BSgenome.Mmusculus.UCSC.mm10)
library(motifmatchr)

#####################
## Define settings ##
#####################

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))
source(here::here("atac/archR/add_motif_annotation/utils.R"))

opts$test <- TRUE

########################
## Load ArchR Project ##
########################

source(here::here("atac/archR/load_archR_project.R"))

######################
## Define variables ##
######################

name = "Motif_cisbp_lenient"
cutOff = 1e-4
width = 7

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

#########################
## Get motif positions ##
#########################

motifPositions <- matchMotifs(
  pwms = motifs,
  subject = ArchRProject@peakSet,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  out = "positions",
  p.cutoff = cutOff,
  w = width
)

# Save
# io$outfile <- sprintf("%s/Annotations/%s_cutOff%s_width%s-Scores.rds",io$archR.directory,name,cutOff,width)
# io$outfile <- sprintf("%s/Annotations/%s-Positions-In-Peaks.rds",io$archR.directory,name)
saveRDS(motifScores, file.path(getOutputDirectory(ArchRProject), paste0(name,"-Scores.rds")), compress = TRUE)

######################
## Get motif scores ##
######################

motifScores <- matchMotifs(
  pwms = motifs[1:10],
  subject = ArchRProject@peakSet,
  genome = BSgenome.Mmusculus.UCSC.mm10, 
  out = "scores", 
  p.cutoff = cutOff, 
  w = width
)

class(assay(motifScores,"motifScores"))
class(assay(motifScores,"motifMatches"))
class(assay(motifScores,"motifCounts"))

# Save
# io$outfile <- sprintf("%s/Annotations/%s_cutOff%s_width%s-Scores.rds",io$archR.directory,name,cutOff,width)
# io$outfile <- sprintf("%s/Annotations/%s-Scores.rds",io$archR.directory,name)
saveRDS(motifScores, file.path(getOutputDirectory(ArchRProject), paste0(name,"-Scores.rds")), compress = TRUE)
  
