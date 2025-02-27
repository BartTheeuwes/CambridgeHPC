########################
## Load ArchR project ##
########################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/atac/archR/load_archR_project.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/atac/archR/load_archR_project.R")
} else {
  stop("Computer not recognised")
}

#####################
## Define settings ##
#####################

# I/O
# io$metadata <- paste0(io$basedir,"/sample_metadata.txt.gz")
# io$outdir <- paste0(io$basedir,"/results/atac/archR/gene_scores")

#####################
## Add tile matrix ##
#####################

addTileMatrix(
  input = ArchRProject,
  tileSize = 500,
  binarize = TRUE,
  excludeChr = c("chrM", "chrY"),
  force = TRUE
)
