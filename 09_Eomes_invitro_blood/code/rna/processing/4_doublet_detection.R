here::i_am("rna/processing/4_doublet_detection.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(scds))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(scater))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--sce',         type="character",                            help='SingleCellExperiment file')
p$add_argument('--metadata',    type="character",                            help='metadata file')
p$add_argument('--samples',                 type="character",   nargs='+',     help='Sample(s)')
p$add_argument('--doublet_score_threshold',  type="double",      default=1.25,   help='Doublet score threshold')
p$add_argument('--outfile',                  type="character",                  help='Output file')
args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# args <- list()
# args$sce <- io$rna.sce
# args$metadata <- paste0(io$basedir,"/results/rna/qc/sample_metadata_after_qc.txt.gz")
# args$samples <- c("2_Eo_DEG_G9_day3_5_VC")
# args$doublet_score_threshold <- 1.25
## END TEST ##

##########################
## Load sample metadata ##
##########################

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
  .[,c("sample","cxds_score", "bcds_score", "hybrid_score")] %>%
  as.data.frame %>% tibble::rownames_to_column("cell") %>% as.data.table %>%
  .[,c("cxds_score","bcds_score","hybrid_score"):=list(round(cxds_score,2),round(bcds_score,2),round(hybrid_score,2))]

# Call doublets
dt[,doublet_call:=hybrid_score>args$doublet_score_threshold]
table(dt$doublet_call)

# Save
# io$outfile <- sprintf("%s/doublets_%s_%s.txt.gz",args$outdir, paste(args$samples,collapse="-"),round(args$doublet_score_threshold,2))
fwrite(dt, args$outfile, sep="\t", na="NA", quote=F)


