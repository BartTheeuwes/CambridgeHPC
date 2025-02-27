
####################
## Load libraries ##
####################
library(Seurat)
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

# Define I/O
io$outdir <- paste0(io$basedir,"/results/rna_atac/rna_vs_acc/trajectories/blood_trajectory")
io$sce_res <-  paste0(io$outdir,"/residual_score.rds")
io$sce_rna <- paste0(io$outdir,"/blood_rna_sce.rds")
io$DORC_score <- paste0(io$outdir,"/DORCs_score_matrix_per_cell.rds")
opts$celltypes = c(
  # "Mixed_mesoderm",
  "Haematoendothelial_progenitors",
  "Blood_progenitors_1",
  "Blood_progenitors_2",
  "Erythroid1",
  "Erythroid2",
  "Erythroid3"
)

########################
## Load data needed ##
########################
res_score <- readRDS(io$sce_res)
sce_rna <- readRDS(io$sce_rna)
DORC_score <- readRDS(io$DORC_score)
#####################
## Load metadata ##
#####################

sample_metadata <- fread(io$metadata) %>%
  .[pass_rnaQC==TRUE&pass_rnaQC==TRUE] %>%
  .[celltype.mapped%in%opts$celltypes & sample%in%opts$samples &cell %in%colnames(sce_rna)] 

##################
## Based on rna ##
##################
meta.data <- as.data.frame(sample_metadata)
rownames(meta.data) <- meta.data$cell
seu_rna <- CreateSeuratObject(counts = counts(sce_rna),meta.data =meta.data[colnames(sce_rna),])
seu_rna <- FindVariableFeatures(seu_rna,nfeatures = 600)
VariableFeaturePlot(seu_rna)
seu_rna <- ScaleData(seu_rna)
seu_rna <- RunPCA(seu_rna)
ElbowPlot(seu_rna)

seu_rna <- FindNeighbors(seu_rna)
seu_rna <- RunUMAP(seu_rna,dims = 1:20)
seu_rna <- RunTSNE(seu_rna,dims = 1:20)


Idents(seu_rna) <- seu_rna$celltype.mapped
markers <- FindAllMarkers(seu_rna,only.pos = T,logfc.threshold = 0.1)
top10 <- markers %>% group_by(cluster) %>% top_n(10,wt=avg_logFC)


rownames(res_score) <- paste0('res_',rownames(res_score))
res_dt <- as.data.frame(t(res_score))

rownames(DORC_score) <- paste0('DORC_',rownames(DORC_score))
DORC_dt <- as.data.frame(t(DORC_score))
seu_rna@meta.data <- cbind(seu_rna@meta.data,res_dt,DORC_dt)
DimPlot(seu_rna,cols = opts$celltype.colors)


p1 <- FeaturePlot(seu_rna,'Foxf1')
p2 <- FeaturePlot(seu_rna,'DORC_Foxf1')
p3 <- DimPlot(seu_rna,cols = opts$celltype.colors)
p1|p2|p3

# res
seu_res <- CreateSeuratObject(counts = res_score,meta.data =meta.data[colnames(sce_rna),])
seu_res <- FindVariableFeatures(seu_res,nfeatures = 600)
VariableFeaturePlot(seu_res)
seu_res <- ScaleData(seu_res)
seu_res <- RunPCA(seu_res)
ElbowPlot(seu_res)

seu_res <- FindNeighbors(seu_res)
seu_res <- RunUMAP(seu_res,dims = 1:15)
seu_res <- RunTSNE(seu_res,dims = 1:15)


Idents(seu_res) <- seu_res$celltype.mapped
markers <- FindAllMarkers(seu_res,only.pos = T)
top10 <- markers %>% group_by(cluster) %>% top_n(10,wt=avg_logFC)


rownames(res_score) <- paste0('res_',rownames(res_score))
res_dt <- as.data.frame(t(res_score))

rownames(DORC_score) <- paste0('DORC_',rownames(DORC_score))
DORC_dt <- as.data.frame(t(DORC_score))
seu_res@meta.data <- cbind(seu_res@meta.data,res_dt,DORC_dt)
DimPlot(seu_res,cols = opts$celltype.colors,reduction = 'tsne')

FeaturePlot(seu_res,'Hba-x',reduction = 'tsne')




