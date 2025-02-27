
####################
## Load libraries ##
####################
library(scAI)
library(dplyr)
library(cowplot)

#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/public_datasets/settings.R")
  source("/Users/ricard/gastrulation_multiome_10x/public_datasets/utils.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/public_datasets/settings.R")
  source("/homes/ricard/gastrulation_multiome_10x/public_datasets/utils.R")
} else {
  stop("Computer not recognised")
}

# Define I/O
io$outdir <- paste0(io$basedir,"/results/results/scAI")
io$umap_atac_coords<- paste0(io$outdir,'/umap_atac_coords.Rdata')
io$umap_atac <- paste0(io$outdir,'/umap_atac.pdf')
io$ranking_plot <- paste0(io$outdir,'/ranking_plot.pdf')
io$markers <- paste0(io$outdir,'/markers.Rdata')
io$obj <- paste0(io$outdir,'scAI_outs.Rdata')
io$obj_final <- paste0(io$outdir,'scAI_outs_final.Rdata')

opts$samples <- c(
  "E7.5_rep1",
  "E7.5_rep2",
  "E8.0_rep1",
  "E8.0_rep2",
  "E8.5_rep1",
  "E8.5_rep2"
)

opts$celltypes = c(
  "Epiblast",
  "Primitive_Streak",
  "Caudal_epiblast",
  "PGC",
  "Anterior_Primitive_Streak",
  "Notochord",
  "Def._endoderm",
  "Gut",
  "Nascent_mesoderm",
  "Mixed_mesoderm",
  "Intermediate_mesoderm",
  "Caudal_Mesoderm",
  "Paraxial_mesoderm",
  "Somitic_mesoderm",
  "Pharyngeal_mesoderm",
  "Cardiomyocytes",
  "Allantois",
  "ExE_mesoderm",
  "Mesenchyme",
  "Haematoendothelial_progenitors",
  "Endothelium",
  "Blood_progenitors_1",
  "Blood_progenitors_2",
  "Erythroid1",
  "Erythroid2",
  "Erythroid3",
  "NMP",
  "Rostral_neurectoderm",
  "Caudal_neurectoderm",
  "Neural_crest",
  "Forebrain_Midbrain_Hindbrain",
  "Spinal_cord",
  "Surface_ectoderm",
  "Visceral_endoderm",
  "ExE_endoderm",
  "ExE_ectoderm",
  "Parietal_endoderm"
)

#####################
## Load metadata ##
#####################

sample_metadata <- fread(io$metadata) %>%
  .[pass_rnaQC==TRUE & doublet_call==FALSE] %>%
  .[celltype.mapped%in%opts$celltypes & sample%in%opts$samples]


###############
## Load data ##
###############

# X is a list of two dgCMatrix named ATAC and RNA, with cells in col


###################
## scAI pipeline ##
###################

# preprocess
scAI_outs <- create_scAIobject(raw.data = X)
scAI_outs <- addpData(scAI_outs, pdata = sample_metadata)
scAI_outs <- preprocessing(scAI_outs, minFeatures = 200, minCells = 1,
                           libararyflag = F, logNormalize = T)


# run scAI (30 factors)
scAI_outs <- run_scAI(scAI_outs, K = 30, nrun = 5)
save(scAI_outs,file=io$obj)

# umap of raw ATAC and aggregated_ATAC
cell_coords.ori <- reducedDims(scAI_outs, data.use = scAI_outs@norm.data$ATAC, do.scale = F, method = "umap", return.object = F)
cell_coords.agg <- reducedDims(scAI_outs, data.use = scAI_outs@agg.data, do.scale = F, method = "umap", return.object = F)
save(cell_coords.ori,cell_coords.agg,file=io$umap_atac_coords)

pdf(io$umap_atac)
gg1 <- cellVisualization(scAI_outs, cell_coords.ori, color.by = "stage",show.legend = F, title = "scATAC-seq")
gg2 <- cellVisualization(scAI_outs, cell_coords.agg, color.by = "stage", ylabel = NULL, title = "Aggregated scATAC-seq")
cowplot::plot_grid(gg1, gg2)
dev.off()

# marker for factors

markers_RNA <- identifyFactorMarkers(scAI_outs, assay = 'RNA', n.top = 5)
markers_ATAC <- identifyFactorMarkers(scAI_outs, assay = 'ATAC', n.top = 5)
save(markers_RNA,markers_ATAC,file=io$markers)
# plot the factors
#£ you should select what you are interested

feature_genes = c('ZSWIM6','BASP1','TXNRD1','NR3C1','CKB','ABHD12','CDH16','NFKBIA','PER1','SCNN1A');
feature_loci <- c('5-17070007-17070367','5-17113334-17113956','5-17118743-17119143','5-17189512-17190061','5-17216210-17217826','5-17249185-17249758',
                  '20-25370461-25372119','16-66928981-66929756','12-6376107-6377003')
feature_loci_motif <- c('FOXA1','RELA/SP1','SP1','SMAD3','NR3C1/YY1','CEBPB/GATA3',
                        'CEBPB','SMAD3','CREB1/NR3C1')
pdf(io$ranking_plot)
featureRankingPlot(scAI_outs, assay = 'RNA', feature.show = feature_genes, top.p = 0.5, ylabel = "Gene score")
featureRankingPlot(scAI_outs, assay = 'ATAC', feature.show = feature_loci, feature.show.names = feature_loci_motif, top.p = 0.5, ylabel = "Locus score")
dev.off()



scAI_outs <- getEmbeddings(scAI_outs)
save(scAI_outs,file=io$obj_final)