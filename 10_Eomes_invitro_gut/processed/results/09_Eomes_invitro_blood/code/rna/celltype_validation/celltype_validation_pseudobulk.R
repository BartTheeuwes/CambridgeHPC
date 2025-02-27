suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(SingleCellExperiment))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--sce',             type="character",                               help='SingleCellExperiment file')
p$add_argument('--marker_genes_file', type="character", help='Marker genes file')
p$add_argument('--outdir',          type="character",                               help='Output file')

args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# args$sce <- "/Users/argelagr/data/gastrulation_histones/results/rna/pseudobulk/SingleCellExperiment_pseudobulk_celltype.mapped_mnn.rds"
# args$marker_genes_file <- "/Users/argelagr/data/gastrulation10x/results/marker_genes/all_stages/marker_genes.txt.gz"
# args$outdir <- "/Users/argelagr/data/gastrulation_histones/results/rna/celltype_validation/pseudobulk/celltype.mapped_mnn"
## END TEST ##

#####################
## Define settings ##
#####################

source(here::here("settings.R"))
source(here::here("utils.R"))

####################################
## Load pseudobulk RNA expression ##
####################################

sce <- readRDS(args$sce)

# Consider only cell types with sufficient number of cells
opts$min.cells <- 30
sce <- sce[,names(which(sce@metadata$n_cells>=opts$min.cells))]

#######################
## Load marker genes ##
#######################

marker_genes.dt <- fread(args$marker_genes_file) %>%
  .[gene%in%rownames(sce)]

##########
## Plot ##
##########

# args$celltypes <- opts$celltypes[1:3]

for (i in colnames(sce)) {
  
  # PCA of marker genes
  genes.to.plot <- unique(marker_genes.dt[celltype==i,gene])
  if (length(genes.to.plot)>75) genes.to.plot <- genes.to.plot %>% sample(size=75)
  sce_tmp <- sce[genes.to.plot,]
  
  annotation_df <- data.frame(row.names = colnames(sce_tmp), celltype=colnames(sce_tmp))
  
  pheatmap(
    mat = logcounts(sce_tmp),
    show_rownames=T, show_colnames=T, 
    cluster_rows=T, cluster_cols = T,
    # scale="row",
    annotation_colors = list("celltype"=opts$celltype.colors[colnames(sce_tmp)]),
    annotation_col = annotation_df,
    annotation_legend = FALSE,
    fontsize_row = 9, fontsize_col = 11,
    treeheight_row = 0, treeheight_col = 0,
    legend = FALSE,
    filename = file.path(args$outdir,sprintf("heatmap_%s_markers.pdf",i)),
    width = 8, height = 12
  )
  
}

# Create a completion token
file.create(file.path(args$outdir,"completed.txt"))
