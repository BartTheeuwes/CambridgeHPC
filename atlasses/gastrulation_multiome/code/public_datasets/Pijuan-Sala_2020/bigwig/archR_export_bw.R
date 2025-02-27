library(ArchR)
library(ggpubr)

#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/public_datasets/Pijuan-Sala_2020/archR/load_archR_project.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/public_datasets/Pijuan-Sala_2020/archR/load_archR_project.R")
} else {
  stop("Computer not recognised")
}
 
###################
## Export bigwig ##
###################

# This function will group, summarize and export a bigwig for each group in an ArchRProject.
getGroupBW(
  ArchRProj = ArchRProject.filt,
  groupBy = "celltype",
  # groupBy = "Sample",
  normMethod = "ReadsInTSS",
  tileSize = 100,
  maxCells = 1000,
  ceiling = 4
)