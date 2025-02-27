suppressPackageStartupMessages(library(argparse))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--sce',             type="character",                               help='SingleCellExperiment file')
p$add_argument('--metadata',        type="character",                               help='Cell metadata file')
p$add_argument('--marker_genes_file', type="character", help='Marker genes file')
p$add_argument('--celltype_label', type="character", help='Cell type label')
p$add_argument('--celltypes', type="character", nargs="+", help='Cell type to consider')
p$add_argument('--umap', type="character", help='Marker genes file')
p$add_argument('--outdir',          type="character",                               help='Output file')

args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# args$sce <- "/bi/group/reik/ricard/data/gastrulation_multiome_10x/processed/rna_new/SingleCellExperiment.rds"
# args$metadata <- "/bi/group/reik/ricard/data/gastrulation_multiome_10x/results_new/rna/mapping/sample_metadata_after_mapping.txt.gz"
# args$marker_genes_file <- "/bi/group/reik/ricard/data/pijuansala2019_gastrulation10x/results/marker_genes/all_stages/marker_genes.txt.gz"
# args$celltype_label <- "celltype.mapped_mnn"
# args$celltypes <- "Epiblast"
# args$umap <- "/bi/group/reik/ricard/data/gastrulation_multiome_10x/results_new/rna/dimensionality_reduction/umap_features1000_pcs30_neigh25_dist0.3.txt.gz"
# args$outdir <- "/bi/group/reik/ricard/data/gastrulation_multiome_10x/results_new/rna/celltype_validation/single_cells"
## END TEST ##

#####################
## Define settings ##
#####################

source(here::here("settings.R"))
source(here::here("utils.R"))

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(args$metadata) %>%
  .[pass_rnaQC==TRUE & doublet_call==FALSE & !is.na(eval(as.name(args$celltype_label)))]

#########################
## Load RNA expression ##
#########################

# Load RNA expression data as SingleCellExperiment object
sce <- load_SingleCellExperiment(args$sce, cells=sample_metadata$cell, normalise = TRUE)

# Add sample metadata as colData
colData(sce) <- sample_metadata %>% tibble::column_to_rownames("cell") %>% DataFrame

#######################
## Load marker genes ##
#######################

marker_genes.dt <- fread(args$marker_genes_file) %>%
  .[gene%in%rownames(sce)]

##############################
## Load UMAP representation ##
##############################

umap.mtx <- fread(args$umap) %>% .[,c("cell","UMAP1","UMAP2")] %>% matrix.please
cells <- intersect(rownames(umap.mtx),colnames(sce))
sce <- sce[,cells]
umap.mtx <- umap.mtx[cells,]
reducedDim(sce,"UMAP") <- umap.mtx

# args$celltypes <- opts$celltypes[1:3]
for (i in args$celltypes) {
  
  # PCA of marker genes
  genes.to.use <- unique(marker_genes.dt[celltype==i,gene])
  sce_tmp <- sce[genes.to.use,]
  sce_tmp <- scater::runPCA(sce_tmp, ncomponents = 2, ntop=nrow(sce_tmp))
    
  p <- plotPCA(sce_tmp, colour_by=args$celltype_label) + 
    scale_color_manual(values=opts$celltype.colors) +
    theme(
      legend.position = "none"
    )
  
  pdf(file.path(args$outdir,sprintf("pca_%s_markers.pdf",i)))
  print(p)
  dev.off()
  
  # Scatterplot of marker genes
  genes.to.plot <- unique(marker_genes.dt[celltype==i & score>=0.85,gene])
  if (length(genes.to.plot)>16) genes.to.plot <- genes.to.plot %>% sample(size=16)
  # genes.to.plot <- c("Hbb-y","Hba-x","Hba-a2","Hba-a1")
  
  p_list <- genes.to.plot %>% map(function(j) {
    plotUMAP(sce, colour_by=j, point_size=0.5) + 
      labs(title=j) +
      scale_color_gradient(low = "gray80", high = "red") +
      ggplot_theme_NoAxes() + 
      theme(
        plot.title = element_text(hjust = 0.5, size=rel(1.25)),
        legend.position = "none"
        )
    }
  )
  
  pdf(file.path(args$outdir,sprintf("umap_%s_markers.pdf",i)))
  print(cowplot::plot_grid(plotlist=p_list))
  dev.off()
  
}
