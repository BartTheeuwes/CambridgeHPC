suppressPackageStartupMessages(library(SingleCellExperiment))

#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
  source("/Users/ricard/gastrulation_multiome_10x/utils.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
  source("/homes/ricard/gastrulation_multiome_10x/utils.R")
} else if (grepl("Workstation",Sys.info()['nodename'])){
  source("/home/lijingyu/gastrulation/gastrulation_multiome_10x/settings.R")
  source("/home/lijingyu/gastrulation/gastrulation_multiome_10x/utils.R")
} else {
  stop("Computer not recognised")
}
io$DORCs.mtx <-paste0(io$outdir,"/DORCs_score_matrix_per_cell.rds")
io$outdir <- paste0(io$basedir,"/results/rna/dimensionality_reduction")
io$umap <- paste0(io$outdir,'/E7.5_rep1-E7.5_rep2-E8.0_rep1-E8.0_rep2-E8.5_rep1-E8.5_rep2_umap_nfeatures50000_ndims50_neigh45_dist0.45.txt.gz')

###############
## Load data ##
###############

# Load metadata #
sample_metadata <- fread(io$metadata) %>%
  .[pass_rnaQC==TRUE & doublet_call==FALSE] %>%
  .[celltype.mapped%in%opts$celltypes & sample%in%opts$samples]

# Load SingleCellExperiment object
sce <- readRDS(io$rna.sce)[,as.character(sample_metadata$cell)]

# Remove genes that are not expressed
sce <- sce[rowMeans(counts(sce))>0,]

# Load DORC score
DORCs.mtx <-readRDS(io$DORCs.mtx) %>% t()%>%
  as.data.table(keep.rownames = T)%>%
  setnames("rn","cell")
umap.dt

# Load gene metadata
gene_metadata <- fread(io$gene_metadata) %>%
  .[symbol%in%rownames(sce)]

# Get UMAP coordinates 
umap.dt <-fread(io$umap)
DORC_umap=merge(umap.dt,DORCs.mtx)


###################################
## Plot dimensionality reduction ##
###################################

# plot celltype
umap.dt<-merge(umap.dt,sample_metadata)
ggplot(umap.dt, aes(x=umap1, y=umap2, color=celltype.mapped)) +
  geom_point(size=0.25) +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank()
  )

# plot gene expr and DORC score
genes.to.plot <- c("Sox17")

for (i in 1:length(genes.to.plot)) {
  gene <- genes.to.plot[i]
  print(sprintf("plot %s:",gene))
  ggplot(DORC_umap, aes(x=umap1, y=umap2, color=Slc40a1)) +
    scale_color_gradient(low = "gray80", high = "red") +
    geom_point(size=0.25) +
    theme_classic() +
    ggtitle(paste0('DORC score of ',gene))+
    theme(plot.title = element_text(hjust = 0.5))+
    theme(
      legend.title = element_blank(),
      legend.position = "none",
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank()
    )
  
  # Plot DORC score
  p1 <- ggplot(DORC_umap, aes(x=umap1, y=umap2, color=Slc40a1)) +
    scale_color_gradient(low = "gray80", high = "red") +
    geom_point(size=0.25) +
    theme_classic() +
    ggtitle(paste0('DORC score of ',gene))+
    theme(plot.title = element_text(hjust = 0.5))+
    theme(
      legend.title = element_blank(),
      legend.position = "none",
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank()
    )
  # plot epxr
  ## Create data.table to plot
  to.plot <- data.table(
    cell = colnames(sce),
    expr = counts(sce[gene,])[1,]
  ) %>% merge(umap.dt,by.x='cell',by.y='cell')
  p2 <-ggplot(to.plot, aes(x=umap1, y=umap2, color=expr)) +
    scale_color_gradient(low = "gray80", high = "red") +
    geom_point(size=0.25) +
    theme_classic() +
    ggtitle(paste0('expression of ',gene))+
    theme(plot.title = element_text(hjust = 0.5))+
    theme(
      legend.title = element_blank(),
      legend.position = "none",
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank()
    )
  p1|p2
  # pdf(sprintf("%s/%s.pdf",io$outdir,i), width=8, height=3.5, useDingbats = F)
  # ggsave("ggtest.png", width = 3.25, height = 3.25, dpi = 1200)
  jpeg(sprintf("%s/%s_umap.jpeg",io$outdir,gene), width = 600, height = 600)
  print(p)
  dev.off()
}
