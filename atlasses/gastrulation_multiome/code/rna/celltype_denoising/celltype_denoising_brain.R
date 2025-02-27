#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
  source("/Users/ricard/gastrulation_multiome_10x/rna/celltype_denoising/utils.R")
} else {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
  source("/homes/ricard/gastrulation_multiome_10x/rna/celltype_denoising/utils.R")
}

# io$metadata <- paste0(io$basedir,"/results/rna/doublets/sample_metadata_after_doublets.txt.gz")
io$metadata.out <- paste0(io$basedir,"/results/rna/celltype_denoising/sample_metadata_after_celltype_denoising.txt.gz")
io$outdir <- paste0(io$basedir,"/results/rna/celltype_denoising")

opts$samples <- c(
  "E7.5_rep1",
  "E7.5_rep2",
  "E8.5_rep1",
  "E8.5_rep2"
)

opts$celltypes = c(
  "Rostral_neurectoderm",
  "Caudal_neurectoderm",
  "Neural_crest",
  "Forebrain_Midbrain_Hindbrain",
  "Spinal_cord"
)

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(io$metadata) %>%
  .[pass_rnaQC==TRUE & doublet_call==FALSE] %>%
  .[celltype.mapped%in%opts$celltypes & sample%in%opts$samples]
table(sample_metadata$celltype.mapped)

###############
## Load data ##
###############

sce <- load_SingleCellExperiment(io$sce, normalise = TRUE, cells = sample_metadata$cell)
dim(sce)

# Add sample metadata to the colData of the SingleCellExperiment
colData(sce) <- sample_metadata %>% as.data.frame %>% tibble::column_to_rownames("cell") %>%
  .[colnames(sce),] %>% DataFrame()

#######################
## Feature selection ##
#######################

decomp <- modelGeneVar(sce)
decomp <- decomp[decomp$mean > 0.01,]

hvgs <- decomp[order(decomp$FDR),] %>% head(n=1500) %>% rownames
# hvgs <- rownames(decomp)[decomp$p.value <= 0.05]

# Subset SingleCellExperiment
sce_filt <- sce[hvgs,]
dim(sce_filt)

##############################
## Dimensionality reduction ##
##############################

# PCA
sce_filt <- runPCA(sce_filt, ncomponents = 15, ntop = 1500)

# UMAP
set.seed(42)
sce_filt <- runUMAP(sce_filt, dimred="PCA", n_neighbors = 15, min_dist = 0.25)

##########
## Plot ##
##########

plotUMAP(sce_filt, colour_by="celltype.mapped") +
  scale_color_manual(values=opts$celltype.colors) + 
  theme(
    legend.position = "none"
  )

##############################################################
## Construct a similarity matrix using the UMAP coordinates ##
##############################################################

dimred.dt <- reducedDim(sce_filt,"PCA") %>% as.data.table %>%
  .[,cell:=colnames(sce_filt)] %>%
  merge(sample_metadata[,c("cell","celltype.mapped")],by="cell")

# Euclidean similarity 
S <- 1 / as.matrix(dist(reducedDim(sce_filt,"PCA")))
rownames(S) <- colnames(S) <- dimred.dt$cell
diag(S) <- 0

################
## Define kNN ##
################

opts$k <- 5

# get closest k cells
knn.matrix <- apply(S, 1, function(x) which.maxn(x,n=opts$k))
k.mapped  <- apply(knn.matrix, 2, function(x) dimred.dt$cell[x])

# get celltypes for the k nearest cells
celltypes <- apply(k.mapped, 2, function(x) dimred.dt$celltype.mapped[match(x, dimred.dt$cell)])

# For each cell calculate the most dominant cell type
celltype.denoised <- apply(celltypes, 2, function(x) getmode(x, 1:length(x)))

sce_filt$celltype.denoised <- celltype.denoised

###########################
## Plot before denoising ##
###########################

plotUMAP(sce_filt, colour_by="celltype.mapped") +
  scale_color_manual(values=opts$celltype.colors) + 
  theme(
    legend.position = "none"
  )

pdf(paste0(io$outdir,"/umap_before_denoising.pdf"), width=6, height=4)
print(p)
dev.off()


###########################
## Plot after denoising ##
###########################

plotUMAP(sce_filt, colour_by="celltype.denoised") +
  scale_color_manual(values=opts$celltype.colors) + 
  theme(
    legend.position = "none"
  )

pdf(sprintf("%s/umap_after_denoising_k%s.pdf",io$outdir,opts$k), width=6, height=4)
print(p)
dev.off()

############################
## Update sample metadata ##
############################

sample_metadata.updated <- fread(io$metadata) %>%
  merge(dimred.dt[,c("cell","celltype.denoised")],by=c("cell"), all.x=T) %>%
  .[is.na(celltype.denoised),celltype.denoised:=celltype.mapped]

# fwrite(sample_metadata.updated, io$metadata.out, sep="\t", na="NA", quote=F)
