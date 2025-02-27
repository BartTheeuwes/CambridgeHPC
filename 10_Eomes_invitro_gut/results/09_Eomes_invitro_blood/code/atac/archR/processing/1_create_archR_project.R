here::i_am("atac/archR/processing/1_create_archR_project.R")
source(here::here("settings.R"))

suppressPackageStartupMessages(library(ArchR))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--arrow_files',     type="character",  nargs='+',      help='Arrow files')
p$add_argument('--genome',          type="character", default="mm10",      help='Genome')
p$add_argument('--outdir',          type="character",                               help='Output directory')

args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# args$arrow_files <- file.path(io$basedir,sprintf("processed/atac/archR/%s.arrow",opts$samples))
# args$genome <- "mm10"
# args$threads <- 1
# args$outdir <- file.path(io$basedir,"processed/atac/archR")
## END TEST ##

#####################
## Define settings ##
#####################

addArchRGenome(args$genome)

############################
## create an ArchRProject ##
############################

ArchRProject <- ArchRProject(
  ArrowFiles = args$arrow_files, 
  outputDirectory = args$outdir,
  copyArrows = FALSE
)

##########
## Save ##
##########

saveArchRProject(ArchRProject)