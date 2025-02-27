
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(argparse))

here::i_am("rna_atac/mofa/TFs/prepare_mofa.R")

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--metadata',        type="character",                               help='Cell metadata file')
p$add_argument('--stages',       type="character",  default="all",  nargs='+',  help='Stages to plot')
p$add_argument('--remove_ExE_cells', action="store_true",                                 help='Remove ExE cells?')
p$add_argument('--rna_sce',          type="character",   help='RNA SingleCellExperiment')
p$add_argument('--chromvar_se',          type="character",   help='ATAC chromVAR SummarizedExperiment')
p$add_argument('--outdir',          type="character",                               help='Output directory')
p$add_argument('--binarise',  action="store_true",  help='Binarise ATAC counts?')
p$add_argument('--vars_to_regress', type="character",                nargs='+',     help='Metadata columns to regress out')
p$add_argument('--test',  action="store_true",  help='Test mode?')

args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

source(here::here("settings.R"))
source(here::here("utils.R"))

## START TEST ##
args$metadata <- file.path(io$basedir,"results_new/atac/archR/celltype_assignment/sample_metadata_after_celltype_assignment.txt.gz")
args$rna_sce <- file.path(io$basedir, "processed_new/rna/SingleCellExperiment_CISBP.rds")
args$chromvar_se <- file.path(io$basedir, "results_new/atac/archR/chromvar_chip/chromVAR_chip_deviations_CISBP_archr.rds")
args$stages <- "E8.75"
args$remove_ExE_cells <- TRUE
args$outdir <- file.path(io$basedir, "results_new/rna_atac/mofa/TFs")
args$test <- TRUE
## END TEST ##

dir.create(args$outdir, showWarnings=F)

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(args$metadata) %>%
  .[nFeature_RNA>=2000,pass_rnaQC:=TRUE] %>% .[stage=="E8.7",stage:="E8.75"] %>%    # temporary
  .[nFrags_atac>=3500 & TSSEnrichment_atac>=9,pass_atacQC:=TRUE] %>%    # temporary
  .[pass_atacQC==TRUE & pass_rnaQC==TRUE & doublet_call==FALSE & stage%in%args$stages]

if (args$remove_ExE_cells) {
  print("Removing ExE cells...")
  sample_metadata <- sample_metadata %>%
    .[!celltype.mapped_mnn%in%c("Visceral_endoderm","ExE_endoderm","ExE_ectoderm","Parietal_endoderm")]
}

if (args$test) {
  print("Test mode activated, subsetting number of cells...")
  sample_metadata <- sample_metadata %>% head(n=250)
}

table(sample_metadata$stage)
table(sample_metadata$celltype.mapped_mnn)

####################
## Load ATAC data ##
####################

chromvar.se <- readRDS(args$chromvar_se)

##############################
## Load RNA expression data ##
##############################

rna.sce <- load_SingleCellExperiment(args$rna_sce, normalise = TRUE, cells = sample_metadata$cell)
dim(rna.sce)

# Add sample metadata to the colData of the SingleCellExperiment
colData(rna.sce) <- sample_metadata %>% as.data.frame %>% tibble::column_to_rownames("cell") %>%
  .[colnames(rna.sce),] %>% DataFrame()

##################
## Filter cells ##
##################

cells <- Reduce("intersect",list(colnames(rna.sce),sample_metadata$cell,colnames(chromvar.se)))
rna.sce <- rna.sce[,cells]
chromvar.se <- chromvar.se[,cells]
sample_metadata <- sample_metadata[cell%in%cells]

###########################
## Feature selection RNA ##
###########################

# Feature selection
rna.features <- which(apply(assay(rna.sce,"logcounts"),1,var)>0.01) %>% names

# Fetch matrix
rna.mtx <- assay(rna.sce[rna.features,],"logcounts") %>% as.matrix
dim(rna.mtx)

hist(rna.mtx)

############################
## Feature selection ATAC ##
############################

# Feature selection
atac.features <- which(apply(assay(chromvar.se,"z"),1,var,na.rm=T)>0.01) %>% names

# Transformation
chromvar.mtx <- assay(chromvar.se,"z")[atac.features,]

chromvar.mtx <- chromvar.mtx-min(chromvar.mtx,na.rm=T)
chromvar.mtx <- log(chromvar.mtx+0.5)

mean(is.na(chromvar.mtx))
hist(chromvar.mtx, breaks=50)

######################################
## Regress out technical covariates ##
######################################

###############
## Save data ##
###############

# Matrix::writeMM(rna.mtx, file.path(args$outdir,"rna.mtx"))
# Matrix::writeMM(atac_tfidf.mtx, file.path(args$outdir,"atac_tfidf.mtx"))
# write.table(rownames(rna.mtx), file.path(args$outdir,"rna_features.txt"), quote=F, row.names=F, col.names=F)
# write.table(rownames(atac_tfidf.mtx), file.path(args$outdir,"atac_features.txt"), quote=F, row.names=F, col.names=F)
# write.table(cells, file.path(args$outdir,"cells.txt"), quote=F, row.names=F, col.names=F)
# fwrite(sample_metadata, file.path(args$outdir,"sample_metadata.txt.gz"))

