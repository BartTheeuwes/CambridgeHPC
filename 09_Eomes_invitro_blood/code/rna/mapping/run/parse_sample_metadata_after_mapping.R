here::i_am("rna/mapping/run/parse_sample_metadata_after_mapping.R")

source(here::here("settings.R"))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--metadata',    type="character",  help='Metadata file to use as input')
p$add_argument('--mapping_mnn',    type="character",  nargs="+", help='Results of the MNN mapping')
p$add_argument('--outfile',          type="character",               help='Output file')
args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# args$query_samples <- opts$samples
# args$metadata <- file.path(io$basedir,"results/rna/qc/sample_metadata_after_qc.txt.gz")
# # args$mapping_dir <- file.path(io$basedir,"results/rna/mapping")
# args$mapping_mnn <- file.path(io$basedir,"results/rna/mapping(..)")
# args$mapping_seurat <- file.path(io$basedir,"results/rna/mapping/(..)")
# args$outfile <- file.path(io$basedir,"results/rna/mapping/sample_metadata_after_mapping.txt.gz")
## END TEST ##

###################
## Load metadata ##
###################

sample_metadata <- fread(args$metadata)

##########################
## Load mapping results ##
##########################

mapping_mnn.dt <- args$mapping_mnn %>% map(~ fread(.)) %>% rbindlist

stopifnot(mapping_mnn.dt$cell%in%sample_metadata$cell)

###########
## Merge ##
###########

colnames(mapping_mnn.dt) = c("cell","celltype.mapped_mnn","celltype.score_mnn","celltype_extended.mapped_mnn","celltype_extended.score_mnn","stage.mapped_mnn", "cellstage.score_mnn", "closest.cell_mnn")

to.save <- sample_metadata %>% 
  merge(mapping_mnn.dt,by="cell",all.x=TRUE) %>%
    .[,celltype_genotype:=sprintf("%s-%s",celltype.mapped_mnn,genotype)] %>%
    .[,celltype_extended_genotype:=sprintf("%s-%s",celltype_extended.mapped_mnn,genotype)]

#################
## Save output ##
#################

fwrite(to.save, args$outfile, sep="\t", na="NA", quote=F)
