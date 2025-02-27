suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(scds))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(argparse))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--sce',         type="character",                            help='SingleCellExperiment file')
p$add_argument('--metadata',    type="character",                            help='metadata file')
p$add_argument('--samples',                 type="character",   nargs='+',     help='Sample(s)')
p$add_argument('--hybrid_score_threshold',  type="double",      default=1.25,   help='Hybrid score threshold')
p$add_argument('--test',                    action = "store_true",             help='Testing mode')
p$add_argument('--outdir',                  type="character",                  help='Output directory')
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
} else {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
}

## START TEST ##
# args <- list()
args$sce <- io$rna.sce
args$metadata <- paste0(io$basedir,"/results/rna/mapping/sample_metadata_after_mapping.txt.gz") # .io$metadata
# args$samples <- c("E8.5_rep2")
# args$hybrid_score_threshold <- 1.0
# args$test <- TRUE
## END TEST ##

if (isTRUE(args$test)) print("Test mode activated...")

############################
## Update sample metadata ##
############################

sample_metadata <- fread(args$metadata) %>%
  .[pass_rnaQC==TRUE & sample%in%args$samples]
table(sample_metadata$sample)

###############
## Load data ##
###############

# Load SingleCellExperiment object
sce <- load_SingleCellExperiment(args$sce, cells=sample_metadata$cell, normalise = TRUE)
dim(sce)

# Add sample metadata as colData
colData(sce) <- sample_metadata %>% tibble::column_to_rownames("cell") %>% DataFrame

#############################
## Calculate doublet score ##
#############################

sce <- cxds_bcds_hybrid(sce, estNdbl=TRUE)

dt <- colData(sce) %>%
  .[,c("sample","cxds_score", "cxds_call", "bcds_score", "bcds_call", "hybrid_score", "hybrid_call")] %>%
  as.data.frame %>% tibble::rownames_to_column("cell") %>% as.data.table

# Call doublets
dt[,doublet_call:=hybrid_score>args$hybrid_score_threshold]
table(dt$doublet_call)

##########
## Save ##
##########

io$outfile <- sprintf("%s/%s_%s.txt.gz",args$outdir, paste(args$samples,collapse="-"),round(args$hybrid_score_threshold,2))
fwrite(dt, io$outfile, sep="\t", na="NA", quote=F)
