library(RColorBrewer)
library(ggnewscale)
library(pheatmap)
library(ggpubr)
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
# Options
opts$motif_annotation <- "Motif_cisbp"
opts$celltypes = c(
  # "Mixed_mesoderm",
  "Haematoendothelial_progenitors",
  "Blood_progenitors_1",
  "Blood_progenitors_2",
  "Erythroid1",
  "Erythroid2",
  "Erythroid3"
  # "Cardiomyocytes",
  # 'Pharyngeal_mesoderm',
  # 'Mesenchyme'
)
opts$celltypes <- c(
 
)
opts$min.expr <- 0.1


# I/O

io$outdir <- paste0(io$basedir,"/results/rna_atac/rna_vs_acc/trajectories/blood_trajectory")
io$archR.pseudobulk.deviations.se <- sprintf("%s/pseudobulk/pseudobulk_DeviationMatrix_%s_summarized_experiment.rds",io$archR.directory,opts$motif_annotation)
io$dorc.pseudobulk <- paste0(io$basedir, '/results/rna_atac/DORCs/pseudobulk/DORCs_samples_distance10000.rds')
io$atac.peak.se <- paste0(io$basedir,'/processed/atac/archR/atac_SummarizedExperiment.rds')
io$trajectory <- paste0(io$basedir,"/results/rna/trajectories/blood_trajectory/blood_trajectory.txt.gz")

#####################
## Update metadata ##
#####################

sample_metadata <- fread(io$metadata) %>%
  .[pass_atacQC==TRUE & pass_rnaQC==TRUE& doublet_call==FALSE] %>%
  .[celltype.mapped%in%opts$celltypes & sample%in%opts$samples] %>%
  setnames("celltype.mapped","celltype")

########################
## Load ArchR Project ##
########################


# if (grepl("ricard",Sys.info()['nodename'])) {
#   source("/Users/ricard/gastrulation_multiome_10x/atac/archR/load_archR_project.R")
# } else if (grepl("ebi",Sys.info()['nodename'])) {
#   source("/homes/ricard/gastrulation_multiome_10x/atac/archR/load_archR_project.R")
# } else if(grepl('Workstation',Sys.info()['nodename'])){
#   # source("/home/lijingyu/gastrulation/gastrulation_multiome_10x/atac/archR/load_archR_project.R")
#   ArchRProject=readRDS("/home/lijingyu/gastrulation/data/gastrulation_multiome_10x/processed/atac/archR/Save-ArchR-Project.rds")
# }else{
#   stop("Computer not recognised")
# }


#####################
## Load trajectory ##
#####################

trajectory.dt <- fread(io$trajectory)
# Filter cells
trajectory.dt <- trajectory.dt[cell%in%sample_metadata$cell]
sample_metadata <- sample_metadata[cell%in%trajectory.dt$cell] %>% setkey(cell) %>% .[trajectory.dt$cell]
sample_metadata <- merge(sample_metadata,trajectory.dt)
stopifnot(!is.na(sample_metadata$sample))
stopifnot(sample_metadata$cell==trajectory.dt$cell)

#############################
## Load peak accessibility ##
#############################

if (grepl("ricard",Sys.info()['nodename'])) {
  atac.peak.se <- getMatrixFromProject(ArchRproject, binarize = TRUE, useMatrix = "PeakMatrix")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  atac.peak.se <- getMatrixFromProject(ArchRproject, binarize = TRUE, useMatrix = "PeakMatrix")
} else if(grepl('Workstation',Sys.info()['nodename'])) {
  atac.peak.se <- readRDS(io$atac.peak.se)
} else{
  stop("Computer not recognised")
}
# dosnt wkr???
# atac.peak.se[,sample_metadata$archR_cell]

# Define peak names
peak_names <- rowRanges(atac.peak.se) %>% as.data.table %>% .[,id:=sprintf("%s:%s-%s",seqnames,start,end)] %>% .$id
rownames(atac.peak.se) <- peak_names

# Rename cells
colnames(atac.peak.se) <- sample_metadata %>%
  .[cell%in%colnames(atac.peak.se)] %>%
  setkey(cell) %>% .[colnames(atac.peak.se)] %>% .$cell

# Filter cells
atac.peak.se <- atac.peak.se[,sample_metadata$cell]


# Create peak matrix 
# peak.mtx <- assay(atac.peak.se)
##############################################
## Load association between peaks and genes ##
##############################################


DORCs.granges <- readRDS(io$dorc.pseudobulk)


genes <- unique(DORCs.granges$gene) 
DORCs.mtx <- matrix(NA, nrow=length(genes), ncol=length(colnames(atac.peak.se)))
rownames(DORCs.mtx) <- genes; colnames(DORCs.mtx) <- colnames(atac.peak.se)

for (i in genes) {
  peaks <- DORCs.granges[DORCs.granges$gene==i]$peak %>% as.character
  DORCs.mtx[i,] <- colMeans(assay(atac.peak.se[peaks,]))
}
saveRDS(DORCs.mtx,file=paste0(io$outdir,"/DORCs.mtx_raw.rds"))


###############################
## Load SingleCellExperiment ##
###############################

sce <- load_SingleCellExperiment(
  file = io$rna.sce, 
  cells = sample_metadata$cell,
  normalise = TRUE, 
  remove_non_expressed_genes = FALSE
)
sce <- sce[rownames(DORCs.mtx)]

# add metadata

# colData <- merge(as.data.frame(sce@colData),sample_metadata)
# rownames(colData) <- colData$cell
# colData <- colData[colnames(sce),]
# sce@colData <-DataFrame(colData)


#############
## Denoise ##
#############

trajectory.mtx <- trajectory.dt[,c("cell","PC1")] %>% matrix.please

opts$knn <- 1
# RNA

rna.mtx <- smoother_aggregate_nearest_nb(mat=as.matrix(logcounts(sce)), D=pdist(trajectory.mtx), k=opts$knn)
colnames(rna.mtx) <- colnames(sce)
saveRDS(rna.mtx,paste0(io$outdir,"/rna_smoothed.rds"))

# DORC

DORCs.mtx <- smoother_aggregate_nearest_nb(mat=DORCs.mtx, D=pdist(trajectory.mtx), k=opts$knn)
colnames(DORCs.mtx) <- colnames(atac.peak.se)
# save
saveRDS(DORCs.mtx, paste0(io$outdir,"/DORCs_score.rds"))

# ATAC
# atac.peak.mtx <- smoother_aggregate_nearest_nb(mat=as.matrix(assay(atac.peakMatrix.se)), D=pdist(trajectory.mtx), k=opts$knn)
# colnames(atac.peak.mtx) <- colnames(atac.peakMatrix.se)


############################
## Compute residual score ##
############################
sce_norm <- minmax.normalisation(rna.mtx[rownames(DORCs.mtx),])
saveRDS(sce_norm,paste0(io$outdir,"/rna_score.rds"))
res_score <- DORCs.mtx - sce_norm
saveRDS(res_score,paste0(io$outdir,"/residual_score.rds"))

DORC.dt <- DORCs.mtx %>% t %>%
  as.data.table(keep.rownames = T) %>%
  setnames("rn","cell") %>%
  melt(id.vars=c("cell"), variable.name="gene", value.name="DORC_score")

rna.dt <- rna.mtx %>%
  as.data.table(keep.rownames = T) %>%
  setnames("rn","gene") %>%
  melt(id.vars="gene", variable.name="cell", value.name="expr")

DORC_rna_dt <- merge(
  rna.dt,
  DORC.dt %>% 
    merge(sample_metadata[,c("cell","celltype",'PC1','DC1')], by="cell"),
  by = c("cell","gene")
)
fwrite(DORC_rna_dt, sprintf("%s/DORC_rna_knn_5_1.txt.gz",io$outdir), quote=F, sep="\t", na="NA")


#plot RNA
# "Mrap",'Hemgn','Trim10'
gene <- 'Trim10'
#Create data.table to plot
to.plot <- data.table(
  cell = colnames(rna.mtx),
  expr = rna.mtx[gene,]
) %>% merge(sample_metadata, by="cell")

to.plot$celltype <- factor(to.plot$celltype,levels = c('Haematoendothelial_progenitors','Blood_progenitors_1','Blood_progenitors_2','Erythroid1','Erythroid2','Erythroid3'))
# Plot

p <- ggplot(to.plot, aes(x=celltype, y=expr, fill=celltype)) +
  geom_violin(scale = "width") +
  geom_boxplot(width=0.5, outlier.shape=NA) +
  scale_fill_manual(values=opts$celltype.colors) +
  labs(title=gene) +
  # facet_wrap(~class) +
  theme_classic() +
  labs(x="",y="RNA expression") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(colour="black",size=rel(1.0), angle=50, hjust=1),
    axis.text.y = element_text(colour="black",size=rel(1.0)),
    legend.position="none"
  )
# plot DORC
to.plot[, DORC :=DORCs.mtx[gene,colnames(sce)] ]
p1 <- ggplot(to.plot, aes(x=celltype, y=DORC, fill=celltype)) +
  geom_violin(scale = "width") +
  geom_boxplot(width=0.5, outlier.shape=NA) +
  scale_fill_manual(values=opts$celltype.colors) +
  labs(title=gene) +
  # facet_wrap(~class) +
  theme_classic() +
  labs(x="",y="DORC score") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(colour="black",size=rel(1.0), angle=50, hjust=1),
    axis.text.y = element_text(colour="black",size=rel(1.0)),
    legend.position="none"
  )
ggarrange(p1,p,ncol = 1,
          common.legend = T,
          legend = "right")
