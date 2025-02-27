suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(cisTopic))

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--samples',         type="character",                nargs='+',     help='Samples')
p$add_argument('--nfeatures',       type="integer",    default=1000,               help='Number of features')
p$add_argument('--ntopics',         type="integer",      nargs="+",    help='Number of topics')
p$add_argument('--indir',           type="character",                               help='Input directory')
p$add_argument('--outdir',          type="character",                               help='Output directory')
args <- p$parse_args(commandArgs(TRUE))

## START TEST
args$samples <- opts$samples
args$nfeatures <- 5e4
args$ntopics <- 50
args$colour_by <- c("celltype.mapped","sample","stage","log_nFrags_atac")
args$indir <- paste0(io$basedir,"/results/atac/cistopic/all_cells")
args$outdir <- paste0(io$basedir,"/results/atac/cistopic/all_cells")
## END TEST


#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
  source("/Users/ricard/gastrulation_multiome_10x/utils.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
  source("/homes/ricard/gastrulation_multiome_10x/utils.R")
} else {
  stop("Computer not recognised")
}

# I/O

# Options

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(io$metadata) %>%
  .[pass_atacQC==TRUE & doublet_call==FALSE & sample%in%args$samples] %>%
  .[,log_nFrags_atac:=log10(nFrags_atac)]

#####################
## Load input data ##
#####################

#########################
## Load cisTopic model ##
#########################

# io$input.file <- sprintf("%s/%s_cistopic_features%d_ntopics%d.txt.gz",args$indir, paste(args$samples,collapse="-"), args$nfeatures, args$ntopics)
cistopic.model <- readRDS("/Users/ricard/data/gastrulation_multiome_10x/results/atac/cistopic/E8.5/E8.5_rep1-E8.5_rep2_cistopic_bestModel_features10000_ntopics50.rds")
get

topics.mtx <- modelMatSelection(cisTopicObject, 'cell', opts$output) %>% t 
# foo <- t(cistopic.model@selected.model$topics) %*% cistopic.model@selected.model$document_expects

p_mat <- predictiveDistribution(cistopic.model)
