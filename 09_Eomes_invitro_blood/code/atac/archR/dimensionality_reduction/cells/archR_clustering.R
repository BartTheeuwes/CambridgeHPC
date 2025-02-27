here::i_am("atac/archR/dimensionality_reduction/cells/archR_dimensionality_reduction.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

suppressPackageStartupMessages(library(ArchR))
suppressPackageStartupMessages(library(uwot))


######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--metadata',        type="character",                               help='Cell metadata file')
p$add_argument('--reduced_dim',          type="character",  default="Reduced dimension rds file",   help='Reduced dimension')
p$add_argument('--matrix',          type="character",  default="matrix that the reduced dimension was made from",   help='Reduced dimension matrix')
p$add_argument('--cluster_resolution',        type="double",     default=0.5,  nargs='+',     help='resolution for clustering')
p$add_argument('--outfile',          type="character",                               help='Output file')

args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# args <- list()
# args$metadata <- "/bi/group/reik/ricard/data/gastrulation_multiome_10x/results/atac/archR/qc/sample_metadata_after_qc.txt.gz"
# args$reduced_dim <- 15000
# args$outdir <- "/bi/group/reik/ricard/data/gastrulation_multiome_10x/results/atac/archR/dimensionality_reduction"
## END TEST ##

###############################################
## Load sample metadata & reduced dimensions ##
###############################################

sample_metadata <- fread(args$metadata) %>%
  .[pass_atacQC==TRUE & doublet_call==FALSE | pass_atacQC==TRUE & is.na(doublet_call)]

reduced_dim = readRDS(args$reduced_dim)

########################
## Load ArchR project ##
########################

source(here::here("atac/archR/load_archR_project.R"))

# Subset
ArchRProject.filt <- ArchRProject[sample_metadata$cell]

################
## Clustering ##
################

ArchRProject.filt@reducedDims$IterativeLSI = reduced_dim

ArchRProject.filt <- addClusters(
    input = ArchRProject.filt,
    reducedDims = "IterativeLSI",
    method = "Seurat",
    name = "Clusters",
    resolution = args$cluster_resolution
)

meta_updated = sample_metadata %>% 
    .[,Clusters:=ArchRProject.filt@cellColData$Clusters] %>%
    .[,paste0("Clusters_", args$matrix, "_genotype"):=paste0(Clusters,genotype)] %>%
    setnames('Clusters', paste0("Clusters_", args$matrix, args$cluster_resolution))

fwrite(meta_updated, args$outfile) 

# Save an extra file indicating that the code finished, only necessary for snakemake to know the order
completed_file = paste0(strsplit(args$outfile, '\\.') %>% map_chr(1), 
                    '_', args$matrix,
                    '_completed.txt.gz')
fwrite(data.table(), completed_file)
