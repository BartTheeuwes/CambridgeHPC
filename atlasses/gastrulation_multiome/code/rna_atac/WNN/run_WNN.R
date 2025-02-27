library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(scater))



# Load default settings
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

########################
## Load ArchR project ##
########################


if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/atac/archR/load_archR_project.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/atac/archR/load_archR_project.R")
} else if(grepl('Workstation',Sys.info()['nodename'])){
  # source("/home/lijingyu/gastrulation/gastrulation_multiome_10x/atac/archR/load_archR_project.R")
  ArchRProject=readRDS("/home/lijingyu/gastrulation/data/gastrulation_multiome_10x/processed/atac/archR/Save-ArchR-Project.rds")
}else{
  stop("Computer not recognised")
}
################
## Define I/O ##
################

io$atac.peaks.se <- paste0(io$basedir,'/processed/atac/archR/atac_SummarizedExperiment.rds')
# Options
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
  # "PGC",
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

opts$npeaks <- 1e4

########################
## Load cell metadata ##
########################

io$metadata <- paste0(io$basedir, "/sample_metadata.txt.gz")

sample_metadata <- fread(io$metadata) %>%
  .[(pass_atacQC==TRUE & pass_rnaQC==TRUE) & doublet_call==FALSE] %>%
  .[sample%in%opts$samples & celltype.mapped%in%opts$celltypes]

##########################
## Subset ArchR project ##
##########################

ArchRProject.filt <- ArchRProject[sample_metadata$cell]

######################
## Load peak matrix ##
#####################

# Get Peak Matrix from ArchR
# atac.peak.se <- getMatrixFromProject(ArchRProject.filt, useMatrix="PeakMatrix")
if (grepl("ricard",Sys.info()['nodename'])) {
  atac.peaks.se <- getMatrixFromProject(ArchRproject, binarize = TRUE, useMatrix = "PeakMatrix")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  atac.peaks.se <- getMatrixFromProject(ArchRproject, binarize = TRUE, useMatrix = "PeakMatrix")
} else if(grepl('Workstation',Sys.info()['nodename'])) {
  atac.peaks.se <- readRDS(io$atac.peaks.se)
}else{
  stop("Computer not recognised")
}
atac.peaks.se <- atac.peaks.se[,sample_metadata$cell]

# Define peak names
peak_names <- rowRanges(atac.peaks.se) %>% as.data.table %>% .[,id:=sprintf("%s_%s_%s",seqnames,start,end)] %>% .$id
rownames(atac.peaks.se) <- peak_names

#####################################
## Feature selection on ATAc peaks ##
#####################################

# Load peak variability estimates
peak.variability.dt <- fread(io$archR.peak.variability) %>%
  .[peak%in%peak_names]

# Define highly variable peaks
peaks <- peak.variability.dt %>% 
  setorder(-variance_pseudobulk) %>% 
  head(n=opts$npeaks) %>% 
  .$peak

# Subset
atac.peaks.se.filt <- atac.peaks.se[peaks,]
dim(atac.peaks.se.filt)

# Get matrix
peak.mtx <- assay(atac.peaks.se.filt, "PeakMatrix")

# TFIDF normalisation
tfidf.peak.mtx <- tfidf(peak.mtx)
dim(tfidf.peak.mtx)

# hist(tfidf.peak.mtx[1:10000,1:1000] %>% as.matrix)

##############################
## Load RNA expression data ##
##############################

sce <- load_SingleCellExperiment(io$rna.sce, normalise = TRUE, cells = sample_metadata$cell)
dim(sce)

# Add sample metadata to the colData of the SingleCellExperiment
colData(sce) <- sample_metadata %>% as.data.frame %>% tibble::column_to_rownames("cell") %>%
  .[colnames(sce),] %>% DataFrame()

# Feature selection
decomp <- modelGeneVar(sce)
decomp <- decomp[decomp$mean > 0.001,]
hvgs <- rownames(decomp)[decomp$p.value <= 0.10]
sce_filt <- sce[hvgs,]


seu <- CreateSeuratObject(sce_filt@assays@data$counts)
# seu@assays$RNA@data <- logcounts(sce_filt)
VariableFeatures(seu) <- rownames(seu)
seu <- ScaleData(seu)
seu <- RunPCA(seu)

atac.data <- atac.peaks.se.filt@assays@data$PeakMatrix
rownames(atac.data) <- rownames(atac.peaks.se.filt)
rownames(atac.peaks.se.filt)
seu[["ATAC"]] <- CreateChromatinAssay(
  counts = atac.data,
  sep = c("_",'_'),
)
DefaultAssay(seu) <- "ATAC"
seu <- RunTFIDF(seu)
seu <- RunSVD(seu)
seu <- FindMultiModalNeighbors(
  object = seu,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:50, 2:40),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)

# build a joint UMAP visualization
seu <- RunUMAP(
  object = seu,
  nn.name = "weighted.nn",
  assay = "RNA",
  verbose = TRUE
)


# add metadata
metadata <- as.data.frame(sce_filt@colData@listData)
rownames(metadata) <- colnames(sce_filt)
metadata <- cbind(metadata,seu@meta.data)
seu@meta.data <- metadata
Idents(seu) <- seu$celltype.mapped

# plot
DimPlot(seu, label = TRUE, repel = TRUE, reduction = "umap",group.by = 'celltype.mapped',cols = opts$celltype.colors ) + NoLegend()

VlnPlot(seu,features  = 'RNA.weight')
# seu@meta.data <- as.data.frame( sce_filt@colData@listData)
