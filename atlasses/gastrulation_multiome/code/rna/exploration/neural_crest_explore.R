here::i_am("rna/exploration/neural_crest_explore.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(scran))

#####################
## Define settings ##
#####################

io$basedir <- file.path(io$basedir,"test")
io$metadata <- file.path(io$basedir,"results/atac/archR/qc/sample_metadata_after_qc.txt.gz")
io$sce <- file.path(io$basedir,"processed/rna/SingleCellExperiment.rds") # io$rna.sce
io$outdir <- paste0(io$basedir,"/results/rna/dimensionality_reduction/test"); dir.create(io$outdir, showWarnings = F)

opts$celltypes <- c("Neural_crest") # c("E8.5_CRISPR_T_WT", "E8.5_CRISPR_T_KO")

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(io$metadata) %>%
  .[pass_rnaQC==TRUE & doublet_call==FALSE & stage%in%opts$stages & sample%in%opts$samples & celltype%in%opts$celltypes]

table(sample_metadata$stage)
table(sample_metadata$sample)
table(sample_metadata$celltype)

###############
## Load data ##
###############

# Load RNA expression data as SingleCellExperiment object
sce <- load_SingleCellExperiment(io$sce, cells=sample_metadata$cell, normalise = TRUE, remove_non_expressed_genes = T)

# Add sample metadata as colData
colData(sce) <- sample_metadata %>% tibble::column_to_rownames("cell") %>% DataFrame

#######################
## Feature selection ##
#######################

# decomp <- modelGeneVar(sce)
# decomp <- decomp[decomp$mean > 0.01,]
# hvgs <- decomp[order(decomp$FDR),] %>% head(n=1000) %>% rownames

hvgs <- grep("Hox",rownames(sce),value=T)

# Subset SingleCellExperiment
sce_filt <- sce[hvgs,]

sce_filt <- runPCA(sce_filt, ncomponents = 2)
reducedDim(sce,"PCA") <- reducedDim(sce_filt,"PCA") 

# sce_filt <- runUMAP(sce_filt, dimred="PCA", n_neighbors = 15, min_dist = 0.3)
# reducedDim(sce,"UMAP") <- reducedDim(sce_filt,"UMAP") 

rna.dt <- as.matrix(logcounts(sce_filt)) %>% 
  as.data.table(keep.rownames = "gene") %>% 
  melt(id.vars = "gene", variable.name = "cell", value.name = "expr")

##########
## Plot ##
##########

to.plot <- reducedDim(sce_filt,"PCA") %>% as.data.table %>% 
  .[,cell:=colnames(sce_filt)] %>%
  merge(sample_metadata, by="cell") 

ggplot(to.plot, aes_string(x="PC1", y="PC2", fill="stage")) +
  geom_point(size=2, shape=21, stroke=0.05) +
  # scale_fill_manual(values=opts$celltype.colors) +
  theme_classic() +
  theme(
    legend.position = "right"
  )



to.plot2 <- to.plot %>% merge(rna.dt, by="cell", allow.cartesian=T)

ggplot(to.plot2, aes_string(x="PC1", y="PC2", color="expr")) +
  facet_wrap(~gene) +
  geom_point(size=1) +
  scale_color_gradient(low = "gray80", high = "purple") +
  theme_classic() +
  theme(
    legend.position = "right"
  )

for (i in grep("Hox",rownames(sce_filt),value=T)) {
  plotPCA(sce, colour_by=i, ncomponents=c(1,2)) %>% print
}

plotPCA(sce, colour_by="Hox", ncomponents=c(1,2)) %>% print
