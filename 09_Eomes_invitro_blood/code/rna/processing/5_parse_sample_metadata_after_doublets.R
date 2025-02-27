here::i_am("rna/processing/5_parse_sample_metadata_after_doublets.R")

source(here::here("settings.R"))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--metadata',         type="character",   help='Metadata file to use as input')
p$add_argument('--doublet_files',    type="character", nargs="+",  help='Results of the doublet score detection algorithm')
p$add_argument('--outfile',          type="character",   help='Output file')
args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# args$metadata <- file.path(io$basedir,"results/rna/qc/sample_metadata_after_qc.txt.gz")
# args$doublet_files <- file.path(io$basedir,"results/rna/doublet_detection/doublets_AGTCAA_R7_L001_mm10_sorted_merged_rmdup_mtx2_1.25.txt.gz")
# args$outfile <- file.path(io$basedir,"results/rna/doublet_detection/sample_metadata_after_doublets.txt.gz")
## END TEST ##

##########################
## Load mapping results ##
##########################

doublet.dt <- args$doublet_files %>% map(~ fread(.)) %>% rbindlist

####################
## Merge and save ##
####################

to.save <- fread(args$metadata) %>% 
  merge(doublet.dt[,c("cell","hybrid_score","doublet_call")] %>% setnames("hybrid_score","doublet_score"), by="cell", all.x=TRUE)
fwrite(to.save, args$outfile, sep="\t", na="NA", quote=F)