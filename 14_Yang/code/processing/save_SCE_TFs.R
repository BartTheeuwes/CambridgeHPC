suppressMessages(library(argparse))

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--motif_annotation',    type="character",    help='group A')
p$add_argument('--TF_file',    type="character",    help='group A')
p$add_argument('--outfile',   type="character",    help='Output file')
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

## START TEST
args$TF_file <- file.path(io$basedir,"processed_new/atac/archR/Annotations/JASPAR_TFs.txt.gz")
args$motif_annotation <- "JASPAR"
args$outfile <- file.path(io$basedir,sprintf("processed_new/rna/SingleCellExperiment_%s.rds",args$motif_annotation))
## END TEST

#########################
## Load RNA expression ##
#########################

# Load SingleCellExperiment object
rna.sce <- load_SingleCellExperiment(io$rna.sce)

######################
## Load list of TFs ##
######################

TFs <- fread(args$TF_file)[["gene"]]

################
## Subset TFs ##
################

TFs <- intersect(TFs,toupper(rownames(rna.sce)))
rna_tf.sce <- rna.sce[stringr::str_to_title(TFs),]
rownames(rna_tf.sce) <- toupper(rownames(rna_tf.sce))

##########
## Save ##
##########

saveRDS(rna_tf.sce, args$outfile)
