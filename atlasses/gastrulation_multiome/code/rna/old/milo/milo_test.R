# here::i_am("rna/dimensionality_reduction/dimensionality_reduction_sce.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(miloR))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--sce',             type="character",                               help='SingleCellExperiment file')
p$add_argument('--metadata',        type="character",                               help='Cell metadata file')
p$add_argument('--stages',       type="character",  default="all",  nargs='+',  help='Stages to plot')
p$add_argument('--features',        type="integer",    default=1000,                help='Number of features')
p$add_argument('--npcs',            type="integer",    default=30,                  help='Number of PCs')
p$add_argument('--n_neighbors',     type="integer",    default=30,     help='(UMAP) Number of neighbours')
p$add_argument('--min_dist',        type="double",     default=0.3,     help='(UMAP) Minimum distance')
p$add_argument('--remove_ExE_cells',       type="character",  default="False",  help='Remove ExE cells? ("True"/"False")')
p$add_argument('--outdir',          type="character",                               help='Output file')
p$add_argument('--samples',         type="character",  default="all",              nargs='+',     help='Samples')
p$add_argument('--vars_to_regress', type="character",                nargs='+',     help='Metadata columns to regress out')
p$add_argument('--batch_correction',type="character",                               help='Metadata column to apply batch correction on')

args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
io$basedir <- file.path(io$basedir,"test")
args$sce <- file.path(io$basedir,"processed/rna/SingleCellExperiment.rds") # io$rna.sce
# args$metadata <- file.path(io$basedir,"results/rna/mapping/sample_metadata_after_mapping.txt.gz") # io$metadata
args$metadata <- file.path(io$basedir,"results/atac/archR/qc/sample_metadata_after_qc.txt.gz")
args$stages <- "all" # "E7.75"
args$samples <- c("E8.5_CRISPR_T_WT", "E8.5_CRISPR_T_KO")
args$features <- 2500
args$npcs <- 50
args$vars_to_regress <- NULL # c("nFeature_RNA","mitochondrial_percent_RNA")
args$batch_correction <- "sample"
args$remove_ExE_cells <- "True"
args$n_neighbors <- 25
args$min_dist <- 0.5
args$outdir <- paste0(io$basedir,"/results/rna/milo/test")
## END TEST ##

#####################
## Define settings ##
#####################

dir.create(args$outdir, showWarnings = F)

#####################
## Parse arguments ##
#####################

# Options
if (args$stages[1]=="all") {
  args$stages <- opts$stages
} else {
  stopifnot(args$stages%in%opts$stages)
}
if (args$samples[1]=="all") {
  args$samples <- opts$samples
} else {
  stopifnot(args$samples%in%opts$samples)
}

if (args$remove_ExE_cells=="True") {
  args$remove_ExE_cells <- TRUE
} else if (args$remove_ExE_cells=="False") {
  args$remove_ExE_cells <- FALSE 
} else {
  stop('remove_ExE_cells should be "True" or "False"')
}

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(args$metadata) %>%
  .[pass_rnaQC==TRUE & doublet_call==FALSE & stage%in%args$stages & sample%in%args$samples]

if (args$remove_ExE_cells) {
  print("Removing ExE cells...")
  sample_metadata <- sample_metadata %>%
    .[!celltype.mapped%in%c("Visceral_endoderm","ExE_endoderm","ExE_ectoderm","Parietal_endoderm")]
}

table(sample_metadata$stage)
table(sample_metadata$sample)
table(sample_metadata$celltype.mapped)

###################
## Sanity checks ##
###################

stopifnot(args$colour_by %in% colnames(sample_metadata))
# stopifnot(unique(sample_metadata$celltype.mapped) %in% names(opts$celltype.colors))

if (length(args$batch_correction)>0) {
  stopifnot(args$batch_correction%in%colnames(sample_metadata))
  if (length(unique(sample_metadata[[args$batch_correction]]))==1) {
    message(sprintf("There is a single level for %s, no batch correction applied",args$batch_correction))
    args$batch_correction <- NULL
  } else {
    library(batchelor)
  }
}

if (length(args$vars_to_regress)>0) {
  stopifnot(args$vars_to_regress%in%colnames(sample_metadata))
}


###############
## Load data ##
###############

# Load RNA expression data as SingleCellExperiment object
sce <- load_SingleCellExperiment(args$sce, cells=sample_metadata$cell, normalise = TRUE)

# Add sample metadata as colData
colData(sce) <- sample_metadata %>% tibble::column_to_rownames("cell") %>% DataFrame

#######################
## Feature selection ##
#######################

decomp <- modelGeneVar(sce)
decomp <- decomp[decomp$mean > 0.01,]
hvgs <- decomp[order(decomp$FDR),] %>% head(n=args$features) %>% rownames

# Subset SingleCellExperiment
sce_filt <- sce[hvgs,]

############################
## Regress out covariates ##
############################

if (length(args$vars_to_regress)>0) {
  print(sprintf("Regressing out variables: %s", paste(args$vars_to_regress,collapse=" ")))
  logcounts(sce_filt) <- RegressOutMatrix(
    mtx = logcounts(sce_filt),
    covariates = colData(sce_filt)[,args$vars_to_regress,drop=F]
  )
}

############################
## PCA + Batch correction ##
############################

if (length(args$batch_correction)>0) {
  suppressPackageStartupMessages(library(batchelor))
  print(sprintf("Applying MNN batch correction for variable: %s", args$batch_correction))
  pca <- multiBatchPCA(sce_filt, batch = colData(sce_filt)[[args$batch_correction]], d = args$npcs)
  pca.corrected <- reducedMNN(pca)$corrected
  colnames(pca.corrected) <- paste0("PC",1:ncol(pca.corrected))
  reducedDim(sce_filt, "PCA") <- pca.corrected[colnames(sce),]
} else {
  sce_filt <- runPCA(sce_filt, ncomponents = args$npcs, ntop=args$features)
}

##########
## UMAP ##
##########

set.seed(args$seed)
sce_filt <- runUMAP(sce_filt, dimred="PCA", n_neighbors = args$n_neighbors, min_dist = args$min_dist)

##########
## Milo ##
##########

# Create Milo object
milo.obj <- Milo(sce_filt)

# Build a k-nearest neighbour graph
milo.obj <- buildGraph(milo.obj, k = 30, d = 30, reduced.dim = "PCA")

# Define neighbourhoods on a graph
milo.obj <- makeNhoods(milo.obj, prop = 0.1, k = 30, d=30, refined = TRUE, reduced_dims = "PCA")

# Plot neighbourhood size
plotNhoodSizeHist(milo.obj)

# Counting cells in neighbourhoods
# Milo leverages the variation in cell numbers between replicates for the same experimental condition to test for differential abundance.
# This adds to the Milo object a n×m matrix, where n is the number of neighbourhoods and m is the number of experimental samples. Values indicate the number of cells from each sample counted in a neighbourhood
milo.obj <- countCells(milo.obj, meta.data = as.data.frame(colData(milo.obj)), sample="sample")
head(nhoodCounts(milo.obj))

# Defining experimental design
design_matrix.df <- data.frame(colData(milo.obj))[,c("sample", "genotype")] %>% unique

# Convert batch info from integer to factor
design_matrix.df$genotype <- as.factor(design_matrix.df$genotype)
rownames(design_matrix.df) <- design_matrix.df$sample

# Calculate within neighbourhood distances
milo.obj <- calcNhoodDistance(milo.obj, d=30, reduced.dim = "PCA")

# This calculates a Fold-change and corrected P-value for each neighbourhood, which indicates wheather there is significant differential abundance between developmental stages
stopifnot(rownames(design_matrix.df)==colnames(milo.obj@nhoodCounts))
da_results <- testNhoods(milo.obj, design = ~ sample, design.df = design_matrix.df)

# Overlay Milo results on the UMAP
milo.obj <- buildNhoodGraph(milo.obj)

# Plot single-cell UMAP
umap_pl <- plotReducedDim(milo.obj, dimred = "umap", colour_by="stage", text_by = "celltype", text_size = 3, point_size=0.5) +
  guides(fill="none")

## Plot neighbourhood graph
nh_graph_pl <- plotNhoodGraphDA(milo.obj, da_results, layout="umap",alpha=0.1) 
umap_pl + nh_graph_pl + plot_layout(guides="collect")

# While neighbourhoods tend to be homogeneous, we can define a threshold for celltype_fraction to exclude neighbourhoods that are a mix of cell types.
da_results$celltype <- ifelse(da_results$celltype_fraction < 0.7, "Mixed", da_results$celltype)

# visualize the distribution of DA Fold Changes in different cell types
plotDAbeeswarm(da_results, group.by = "celltype")
