
suppressPackageStartupMessages(library(ArchR))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(GenomicRanges))

here::i_am("atac/archR/peak_calling/add_custom_peakSet_archR.R")

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--peak_set',    type="character",    help='peak set file')
p$add_argument('--peakset_name',     type="character",  help='Name of peakset matrix')
p$add_argument('--outdir',    type="character",    help='Output directory')
p$add_argument('--threads',     type="integer",    default=1,    help='Number of threads')

args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

source(here::here("settings.R"))
source(here::here("utils.R"))

## START TEST ##
# args$peak_set <- "/rds/project/rds-SDzz0CATGms/users/bt392/atlasses/gastrulation_multiome/data/processed/atac/archR/PeakCalls/peak_metadata.tsv.gz"
# args$peakset_name 'atlas'
# args$outdir <- file.path(io$archR.directory,"PeakCalls/manual")
## END TEST ##

dir.create(args$outdir, showWarnings=F, recursive=TRUE)

########################
## Load ArchR project ##
########################

source(here::here("atac/archR/load_archR_project.R"))
addArchRThreads(threads = args$threads)

###################
## Load peak set ##
###################

# Create GRanges object
peakset.gr <- fread(args$peak_set) %>% 
  #.[str_detect(peakset$GroupReplicate, paste(unique(ArchRProject.filt@cellColData$celltype), collapse='|'))] %>%  # if only keeping present celltypes
  makeGRangesFromDataFrame(., keep.extra.columns = T)

#####################
## Add peak matrix ##
#####################

ArchRProject = addFeatureMatrix( # Functionally the same as addPeakMatrix but for keeping both
  input = ArchRProject,
  features = peakset.gr,
  matrixName = args$peakset_name, ceiling = 4)

##################
## Save project ##
##################
saveArchRProject(ArchRProject)

# Save PeakSet
saveRDS(peakset.gr, file.path(args$outdir, paste0(args$peakset_name, ".rds")))
