suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(batchelor))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(edgeR))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--atlas_stages',    type="character",   nargs='+',  help='Atlas stage(s)')
p$add_argument('--query_samples',   type="character",   nargs='+',  help='Query batch(es)')
p$add_argument('--query_sce',       type="character",               help='SingleCellExperiment file for the query')
p$add_argument('--atlas_sce',       type="character",               help='SingleCellExperiment file for the atlas')
p$add_argument('--query_metadata',  type="character",               help='metadata file for the query')
p$add_argument('--atlas_metadata',  type="character",               help='metadata file for the atlas')
p$add_argument('--npcs',            type="integer",                 help='Number of principal components')
p$add_argument('--n_neighbours',    type="integer",                 help='Number of neighbours')
p$add_argument('--test',            action = "store_true",          help='Testing mode')
p$add_argument('--outdir',          type="character",               help='Output directory')

args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
  source("/Users/ricard/gastrulation_multiome_10x/rna/iterative_mapping/run/utils.R")
  # io$atlas.marker_genes <- "/Users/ricard/data/gastrulation10x/results/marker_genes/E8.5/marker_genes.txt.gz"
  io$script_load_data <- "/Users/ricard/gastrulation_multiome_10x/rna/iterative_mapping/run/load_data.R"
} else {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
  source("/homes/ricard/gastrulation_multiome_10x/rna/iterative_mapping/run/utils.R")
  # io$atlas.marker_genes <- "/hps/nobackup2/research/stegle/users/ricard/gastrulation10x/results/marker_genes/E8.5/marker_genes.txt.gz"
  io$script_load_data <- "/homes/ricard/gastrulation_multiome_10x/rna/iterative_mapping/run/load_data.R"
}
io$path2atlas <- io$atlas.basedir
io$path2query <- io$basedir
io$outdir <- paste0(io$basedir,"/results/rna/iterative_mapping")
io$tree <- paste0(io$atlas.basedir,"/results/phylogenetic_tree/celltypes/PAGA_distances.csv.gz")

## START TEST ##
args$atlas_stages <- c(
  # "E6.5",
  # "E6.75",
  # "E7.0",
  "E7.25",
  "E7.5",
  "E7.75"
  # "E8.0",
  # "E8.25",
  # "E8.5",
  # "mixed_gastrulation"
)
args$query_samples <- c("E7.5_rep1","E7.5_rep2") # opts$samples
args$query_sce <- io$sce
args$atlas_sce <- io$rna.atlas.sce
# args$query_metadata <- paste0(io$basedir,"/results/rna/qc/sample_metadata_after_qc.txt.gz")
args$query_metadata <- io$metadata
args$atlas_metadata <- io$rna.atlas.metadata
args$test <- FALSE
args$npcs <- 5
args$n_neighbours <- 25
## END TEST ##

if (isTRUE(args$test)) print("Test mode activated...")

###############
## Load data ##
###############

source(io$script_load_data)

########################################################
## Define distance matrix for hierarchical clustering ##
########################################################

opts$celltypes <- which(table(meta_atlas$celltype)>10) %>% names# %>% head(n=4)
# opts$celltypes <- unique(sample_metadata_atlas$celltype) 
# opts$celltypes <- c("Caudal_epiblast","NMP")

dist <- fread(io$tree) %>%
  as.data.frame %>% tibble::column_to_rownames("V1") %>% as.matrix %>%
  .[opts$celltypes,opts$celltypes] %>%
  as.dist

#######################
## Recursive mapping ##
#######################

sce_query$celltype_mapped <- paste(opts$celltypes,collapse="%")

while (any(grepl("%",sce_query$celltype_mapped))) {
  print(table(sce_query$celltype_mapped))
  mapping.dt <- recursive.fn(sce_query, sce_atlas, dist, npcs = args$npcs, k = args$n_neighbours, cosineNorm = FALSE)
  ids <- match(mapping.dt$cell,colnames(sce_query))
  sce_query$celltype_mapped[ids] <- mapping.dt$celltype_mapped
  sce_query$celltype_score[ids] <- mapping.dt$celltype_score
  sce_query$stage_mapped[ids] <- mapping.dt$stage_mapped
  sce_query$stage_score[ids] <- mapping.dt$stage_score
}

##########
## Save ##
##########

mapping.dt <- data.table(
  cell = colnames(sce_query), 
  celltype_mapped = sce_query$celltype_mapped,
  celltype_score = sce_query$celltype_score,
  stage_mapped = sce_query$stage_mapped,
  stage_score = sce_query$stage_score
)

fwrite(mapping.dt, sprintf("%s/%s_iterative_mnn.txt.gz",io$outdir,paste(args$query_samples,collapse="_")), sep="\t")

#############
## Explore ##
#############

foo <- mapping.dt %>% merge(asd[,c("cell","celltype.mapped","celltype.score")] %>% setnames(c("cell","celltype_old","score_old")))
fwrite(foo, sprintf("%s/%s_iterative_mnn_TEST.txt.gz",io$outdir,paste(args$query_samples,collapse="_")), sep="\t")
