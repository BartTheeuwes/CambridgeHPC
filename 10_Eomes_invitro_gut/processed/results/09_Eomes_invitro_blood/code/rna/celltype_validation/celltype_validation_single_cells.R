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
p$add_argument('--umap', type="character", help='umap to plot on')
p$add_argument('--outdir',          type="character",                               help='Output file')

args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# args$sce <- "/bi/group/reik/ricard/data/eomes_10x_multiome/processed/rna_new/SingleCellExperiment.rds"
# args$metadata <- "/bi/group/reik/ricard/data/eomes_10x_multiome/results_new/rna/mapping/sample_metadata_after_mapping.txt.gz"
# args$marker_genes_file <- "/bi/group/reik/ricard/data/pijuansala2019_gastrulation10x/results/marker_genes/all_stages/marker_genes.txt.gz"
# args$celltype_label <- "celltype.mapped_mnn"
# args$celltypes <- "Epiblast"
# args$umap <- "/bi/group/reik/ricard/data/eomes_10x_multiome/results_new/rna/dimensionality_reduction/umap_features1000_pcs30_neigh25_dist0.3.txt.gz"
# args$outdir <- "/bi/group/reik/ricard/data/eomes_10x_multiome/results_new/rna/celltype_validation/single_cells"
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
umap.df = as.data.frame(umap.mtx)

#########################
## Per cell type plots ##
#########################

for (i in args$celltypes) {
  
  # PCA of marker genes
  genes.to.use <- unique(marker_genes.dt[celltype==i & score>=0.85,gene])
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
  if (length(genes.to.plot)>16) genes.to.plot <- genes.to.plot %>% head(16)
  # genes.to.plot <- c("Hbb-y","Hba-x","Hba-a2","Hba-a1")
  
  p_list <- genes.to.plot %>% map(function(j) {
values = as.vector(logcounts(sce[j]))
minmax = (values- min(values)) /(max(values)-min(values))
ggplot(umap.df, aes(UMAP1, UMAP2, col=minmax, alpha = minmax)) + #, size = minmax
      ggrastr::rasterize(geom_point(size=minmax*2), dpi=150) + 
      viridis::scale_color_viridis() + 
      theme_void() +
      labs(title=j) +
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

########################
## Compare cell types ##
########################

marker_any = marker_genes.dt[score>0.85] %>% unique(by='gene')

# Determine cell type score per cell
score_celltypes = lapply(args$celltypes, function(i){
  # select marker genes for celltype  
  genes.to.use <- unique(marker_genes.dt[celltype==i & score>=0.85,gene])
  # Only keep top 20 genes
  if(length(genes.to.use)>30){genes.to.use = genes.to.use[1:30]}
  # subset sce
  sce_tmp <- sce[genes.to.use,]
  
  # Negative markers
  marker_neg = marker_genes.dt[celltype==i & score < 0.35 & gene %in% marker_any$gene, gene]  
  sce_tmp2 <- sce[marker_neg,] 
    
  # calculate celltypescore per cell
  tmp = data.table(cell = colnames(sce), 
                   celltype = colData(sce)$celltype, 
                   celltype_test = i, 
                   score = colMeans(t(scale(t(logcounts(sce_tmp))))),
                   score_neg = colMeans(t(scale(t(logcounts(sce_tmp2)))))) # score is mean of scaled logcounts
    return(tmp)
    }) %>% rbindlist(.) %>% .[,score_tot:=score-score_neg] %>% .[,zscore:=scale(score), by=celltype]  %>% # then finally zscore by celltype again
        .[,zscore_capped:=ifelse(zscore<=-3, -3, ifelse(zscore>=3, 3, zscore))] # cap z-scores

score_celltypes.plot = score_celltypes %>% dcast(cell ~ celltype_test, value.var='score_tot')# %>% matrix.please
umap.plot = merge(score_celltypes.plot, as.data.table(umap.df, keep.rownames=T), by.x='cell', by.y='rn')

  p_list <- args$celltypes %>% map(function(j) {
ggplot(umap.plot, aes_string('UMAP1', 'UMAP2', col=j, alpha=j, size=0.5)) + 
    ggrastr::rasterize(geom_point(), dpi=150) + 
    ggtitle(j) + 
    scale_color_gradient2(low='blue', mid='white', high='red') + 
    theme_void() + theme(legend.position='none', 
                        text=element_text(size=20)) 
      })

pdf(file.path(args$outdir,"umap_celltype_scores.pdf"), height=40, width = 40)
   print(cowplot::plot_grid(plotlist=p_list))
dev.off()

fwrite(score_celltypes.plot, file.path(args$outdir,"celltype_scores.csv"))