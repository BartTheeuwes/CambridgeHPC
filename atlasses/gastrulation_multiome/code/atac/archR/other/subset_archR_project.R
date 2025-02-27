
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

io$outdir <- paste0(io$basedir,"/processed/atac/archR_subset")

# opts$samples <- c(
#   "E7.5_rep1"
#   "E7.5_rep2",
#   "E8.0_rep1",
#   "E8.0_rep2",
#   "E8.5_rep1",
#   "E8.5_rep2"
# )


########################
## Load cell metadata ##
########################

sample_metadata <- fread(io$metadata) %>%
  .[pass_atacQC==TRUE]

##################
## Subset ArchR ##
##################

ArchRProject.subset <- subsetArchRProject(
  ArchRProj = ArchRProject,
  cells = sample_metadata$cell,
  outputDirectory = io$outdir,
  dropCells = TRUE,
  force = TRUE
)
