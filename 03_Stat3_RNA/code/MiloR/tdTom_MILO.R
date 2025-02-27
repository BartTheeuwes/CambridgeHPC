suppressPackageStartupMessages(library(miloR))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(argparse))

here::i_am("processing/1_create_seurat_rna.R")

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

p <- ArgumentParser(description='')
p$add_argument('--stages',             type="character",        nargs='+',                         help='stages to use')
args <- p$parse_args(commandArgs(TRUE))

args$sce <- io$rna.sce
#args$metadata <- io$metadata
args$metadata <- paste0(io$basedir,"/results/rna/mapping/sample_metadata_after_mapping.txt.gz")
args$atlas_metadata <- "/rds/project/rds-SDzz0CATGms/users/bt392/atlasses/extended/sample_metadata.txt.gz"
args$features <- 2500
args$npcs <- 50 # test, put at 50
args$nneighbours = 30 # put at 30
args$remove_ExE_cells = FALSE
args$outdir = paste0(io$basedir,"/results/rna/MiloR")

dir.create(args$outdir, recursive=TRUE, showWarnings = FALSE)

print(args$stages)

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(args$metadata) %>%
   .[pass_rnaQC==TRUE & doublet_call==FALSE]

if (args$remove_ExE_cells) {
  print("Removing ExE cells...")
  sample_metadata <- sample_metadata %>%
    .[!celltype.mapped_mnn%in%c("Visceral endoderm","ExE endoderm","ExE ectoderm","Parietal endoderm")]
}

###############
## Load data ##
###############

# Load RNA expression data as SingleCellExperiment object
sce <- load_SingleCellExperiment(args$sce, cells=sample_metadata$cell, normalise = TRUE)

# Add sample metadata as colData
colData(sce) <- sample_metadata %>% tibble::column_to_rownames("cell") %>% DataFrame

# Get Reduced Dims
# PCA
pca_file = sprintf("%s/results/rna/dimensionality_reduction/sce/pca_features%d_pcs%d.txt.gz", io$basedir, args$features, args$npcs)
pca = fread(pca_file) %>% tibble::column_to_rownames("cell") %>% as.matrix
# Umap
umap_file = sprintf("%s/results/rna/dimensionality_reduction/sce/umap_features%d_pcs%d_neigh30_dist0.5.txt.gz", io$basedir, args$features, args$npcs)
umap = fread(umap_file) %>% tibble::column_to_rownames("cell") %>% as.matrix
# Atlas umap
umap_atlas = fread(args$atlas_metadata)[,c('cell', 'umapX', 'umapY')]
umap_mapped = data.table(actual_cell = colnames(sce), cell = colData(sce)$closest.cell_mnn) %>% 
    .[,idx:=1:nrow(.)] %>%
    merge(umap_atlas, by='cell') %>% 
    .[order(idx), c('actual_cell', 'umapX', 'umapY')] %>% 
    tibble::column_to_rownames("actual_cell") %>% as.matrix

# Add Reduced Dims to SCE object
reducedDims(sce) <- list(pca=pca, 
                         umap=umap,
                         umap_atlas = umap_mapped)

tomato = logcounts(sce['tomato-td',])
test = cbind(colData(sce), t(as.data.frame(as.matrix(tomato)))) %>% as.data.frame()
test$tomato_pos = ifelse(test$tomato>0, TRUE, FALSE)
table(test$tdTom, test$tomato_pos)

# annotate any cell with tomato-td expression as being tdTom+
colData(sce)$tdTom_corr = ifelse(logcounts(sce['tomato-td',]) > 0 | colData(sce)$tdTom==TRUE, TRUE, FALSE)


########################
## Create MILO object ##
########################

Milo_sce <- Milo(sce)
Milo_sce = Milo_sce[,colData(Milo_sce)$stage==args$stages]

colData(Milo_sce)$sample_ko = paste0(colData(Milo_sce)$sample, '_',  colData(Milo_sce)$tdTom_corr)

# Build graph
Milo_sce <- buildGraph(Milo_sce, 
                       k = args$nneighbours, 
                       d = args$npcs, 
                       reduced.dim = "pca")

# Identify neighbourhoods from NN-cells
# lower prop for larger datasets
Milo_sce <- makeNhoods(Milo_sce, 
                       prop = 0.1, 
                       k = args$nneighbours, 
                       d = args$npcs, 
                       refined = TRUE, 
                       reduced_dims = "pca")

# count cells of different samples per neighbourhood
Milo_sce <- countCells(Milo_sce, meta.data = as.data.frame(colData(Milo_sce)), sample="sample_ko")

# Calculate Neighbourhood Distances (most time consuming step)
Milo_sce <- calcNhoodDistance(Milo_sce, d=args$npcs, reduced.dim = "pca")

# Build neighbourhood graph
Milo_sce <- buildNhoodGraph(Milo_sce)

# Save MILO object
outfile = sprintf("%s/processed/%s_Milo_features%d_pcs%d.rds",io$basedir, paste(args$stage, collapse ='_'), args$features, args$npcs)
saveRDS(Milo_sce, outfile)

# Embryo design table
embryo_design <- data.frame(colData(Milo_sce))[,c("sample_ko", "sample", 'tdTom_corr', 'tdTom')]

embryo_design <- distinct(embryo_design)
rownames(embryo_design) <- embryo_design$sample_ko

# Test differential abundance per hood
# design = what to test for, can correct for another variable as well
da_results <- testNhoods(Milo_sce, design = ~ tdTom_corr, design.df = embryo_design, reduced.dim="pca")
head(da_results)

# Find Neighbourhood groups
#da_results <- groupNhoods(Milo_sce, da_results, max.lfc.delta = 0.5, overlap = 20)

# Annotate hoods by celltype
da_results <- annotateNhoods(Milo_sce, da_results, coldata_col = "celltype.mapped_mnn")

# Annotate 'mixed' hoods when there is not one main celltype
da_results$celltype.mapped_mnn <- ifelse(da_results$celltype.mapped_mnn_fraction < 0.7, "Mixed", da_results$celltype.mapped_mnn)

##################
## Plot results ##
##################

options(repr.plot.width=15, repr.plot.height=6)

p1 = ggplot(da_results, aes(PValue)) + 
    geom_histogram(bins=50) + 
    theme_bw()

p2 = ggplot(da_results, aes(logFC, -log10(SpatialFDR))) + 
    geom_point() +
    geom_hline(yintercept = 1) + ## Mark significance threshold (10% FDR)
    theme_bw()

# Actually save this to PDF
outfile <- file.path(args$outdir,sprintf("%s_hood_significance.pdf", paste(args$stage, collapse ='_')))
pdf(outfile, width=10, height=5)
ggarrange(p1, p2)
dev.off()



## Plot single-cell UMAP
umap_pl = plotReducedDim(Milo_sce, dimred = "umap_atlas", colour_by="tdTom_corr", text_by = "celltype.mapped_mnn", 
                          text_size = 3, point_size=0.3, text_colour = "grey30",) +
    scale_color_manual(values=opts$tdTom.color, name = 'tdTom') + 
    theme_void() + 
    theme(legend.position='right') +
    guides(fill="none")

## Plot neighbourhood graph
nh_graph_pl <- plotNhoodGraphDA(Milo_sce, da_results, layout="umap_atlas",alpha=0.1) 

                     
outfile <- file.path(args$outdir,sprintf("%s_MILO_graph.pdf", paste(args$stage, collapse ='_')))
pdf(outfile, width=10, height=5)
    umap_pl + nh_graph_pl +
  plot_layout(guides="collect")
dev.off()

# Plot celltype fractions per hood
ggplot(da_results, aes(celltype.mapped_mnn_fraction)) + 
    geom_histogram(bins=50) + 
    theme_bw()

# plot DA beeswarm

                     
outfile <- file.path(args$outdir,sprintf("%s_beeswarm.pdf", paste(args$stage, collapse ='_')))
pdf(outfile, width=10, height=15)
 plotDAbeeswarm(da_results, group.by = "celltype.mapped_mnn")
dev.off()









