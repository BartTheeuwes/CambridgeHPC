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

# Options
# opts$samples <- c(
#   "E7.5_rep1",
#   "E7.5_rep2",
#   "E8.5_rep1",
#   "E8.5_rep2",
#   "E8.25_PijuanSala"
# )
# opts$min.fragments <- 1000
# opts$filterTSS.score <- 4

# I/O
io$output.directory <- paste0(io$basedir,"/processed/atac/archR_integrated")

# io$arrowFiles <- c(
#   "E7.5_rep1" = paste0(io$basedir,"/processed/atac/archR/ArrowFiles/E7.5_rep1.arrow"),
#   "E7.5_rep2" = paste0(io$basedir,"/processed/atac/archR/ArrowFiles/E7.5_rep2.arrow"),
#   "E8.5_rep1" = paste0(io$basedir,"/processed/atac/archR/ArrowFiles/E8.5_rep1.arrow"),
#   "E8.5_rep2" = paste0(io$basedir,"/processed/atac/archR/ArrowFiles/E8.5_rep2.arrow"),
#   "E8.25_PijuanSala" = paste0(io$pijuansala.archR.directory,"/ArrowFiles/E8.25_PijuanSala.arrow")
# )
io$arrowFiles <- c(
  "E7.5_rep1" = paste0(io$basedir,"/processed/atac/archR_integrated/ArrowFiles/E7.5_rep1.arrow"),
  "E7.5_rep2" = paste0(io$basedir,"/processed/atac/archR_integrated/ArrowFiles/E7.5_rep2.arrow"),
  "E8.5_rep1" = paste0(io$basedir,"/processed/atac/archR_integrated/ArrowFiles/E8.5_rep1.arrow"),
  "E8.5_rep2" = paste0(io$basedir,"/processed/atac/archR_integrated/ArrowFiles/E8.5_rep2.arrow"),
  "E8.25_PijuanSala" = paste0(io$basedir,"/processed/atac/archR_integrated/ArrowFiles/E8.25_PijuanSala.arrow")
)

# ArchR options
addArchRThreads(threads = 1) 
addArchRGenome("mm10")

############################
## create an ArchRProject ##
############################

ArchRProject <- ArchRProject(
  ArrowFiles = io$arrowFiles, 
  outputDirectory = io$output.directory,
  copyArrows = FALSE
)

##########
## Save ##
##########

saveArchRProject(ArchRProject)
# saveArchRProject(ArchRProject, io$output.directory)