here::i_am("rna/conversions/convert_SingleCellExperiment_to_anndata.R")

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

suppressPackageStartupMessages(library(reticulate))

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--python_path',   type="character",    help='Python path for reticulate')
p$add_argument('--metadata',   type="character",    help='Cell metadata')
p$add_argument('--sce',  type="character",              help='SingleCellExperiment input file') 
p$add_argument('--outfile',          type="character",                help='Anndata output file')
p$add_argument('--test_mode',    action="store_true",             help='Test mode? subset data')
args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# args <- list()
# args$python_path <- "/Users/argelagr/opt/anaconda3/envs/main/bin/python"
# args$metadata <- file.path(io$basedir,"results/atac/archR/celltype_assignment/sample_metadata_after_celltype_assignment.txt.gz")
# args$sce <- file.path(io$basedir,"processed/rna/SingleCellExperiment.rds")
# args$outfile <- file.path(io$basedir,"processed/rna/anndata.h5ad")
## END TEST ##

################
## Reticulate ##
################

reticulate::use_python(args$python_path, required = TRUE)

sc <- import("scanpy")

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(args$metadata) %>% 
  .[pass_rnaQC==TRUE & doublet_call==FALSE & !is.na(celltype.mapped_mnn)]
  # .[,c("cell", "sample", "stage", "nFeature_RNA", "mitochondrial_percent_RNA", "ribosomal_percent_RNA", "celltype.mapped")] %>%

if (args$test_mode) {
	print("Test mode activated, subsetting number of cells...")
	sample_metadata <- sample_metadata %>% head(n=100)
}

###############
## Load data ##
###############

# Load RNA expression data as SingleCellExperiment object
sce <- load_SingleCellExperiment(args$sce, cells=sample_metadata$cell, normalise = FALSE)

# Add sample metadata as colData
colData(sce) <- sample_metadata %>% tibble::column_to_rownames("cell") %>% DataFrame

#############################################
## Convert SingleCellExperiment to AnnData ##
#############################################

adata_sce <- sc$AnnData(
    X   = t(counts(sce)),
    obs = as.data.frame(colData(sce)),
    var = data.frame(gene=rownames(sce), row.names=rownames(sce))
)
# adata_sce$obsm$update(umap = reducedDim(sce, "umap"))

adata_sce

##########################
## Parse anndata object ##
##########################

# (TO-DO) Add stage colors
# Add cell type colors
# adata_sce$uns$update(celltype_colors = opts$celltype.colors[sort(unique(as.character(adata_sce$obs$celltype)))])
# adata_sce$uns["celltype_colors"]

##########
## Save ##
##########

adata_sce$write_h5ad(args$outfile, compression="gzip")
