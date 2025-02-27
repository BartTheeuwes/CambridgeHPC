here::i_am("rna_atac/rna_vs_acc/cells/trajectories/run_TFexpr_vs_peakAcc_trajectories.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--sce',        type="character",        help='SingleCellExperiment file')
p$add_argument('--metadata',        type="character",   help='Cell metadata file')
p$add_argument('--motif_annotation',  type="character", default="CISBP", help='Motif annotation')
p$add_argument('--motifmatcher',  type="character",              help='Motif annotation') 
p$add_argument('--trajectory',   type="character",    help='File with the trajectory')
p$add_argument('--trajectory_name',   type="character",    help='Name of the trajectory')
p$add_argument('--pca_rna',   type="character",    help='PCA for the RNA')
p$add_argument('--pca_atac',   type="character",    help='PCA for the RNA')
p$add_argument('--denoise', action="store_true", help='apply kNN denoising?')
p$add_argument('--knn',    type="integer", default=25,    help='Number of kNN')
p$add_argument('--outdir',   type="character",    help='Output directory')
p$add_argument('--test_mode',    action="store_true",             help='Test mode? subset data')
args <- p$parse_args(commandArgs(TRUE))

## START TEST
args <- list()
args$motif_annotation <- "CISBP"
args$motifmatcher <- sprintf("%s/Annotations/%s-Scores.rds",io$archR.directory,args$motif_annotation)
args$trajectory_name <- "nmp"
args$sce <- file.path(io$basedir,sprintf("results/rna/trajectories/%s/%s_SingleCellExperiment.rds",args$trajectory_name,args$trajectory_name))
args$metadata <- file.path(io$basedir,sprintf("results/rna/trajectories/%s/%s_sample_metadata.txt.gz",args$trajectory_name,args$trajectory_name))
args$trajectory <- file.path(io$basedir,sprintf("results/rna/trajectories/%s/%s_trajectory.txt.gz",args$trajectory_name,args$trajectory_name))
args$denoise <- TRUE
args$knn <- 25
args$pca_rna <- file.path(io$basedir,"results/rna/dimensionality_reduction/sce/pca_features2500_pcs50.txt.gz")
args$pca_atac <- file.path(io$basedir,"results/atac/archR/dimensionality_reduction/lsi_PeakMatrix_nfeatures25000_ndims50.txt.gz")
args$outdir <- file.path(io$basedir,sprintf("results/rna_atac/rna_vs_acc/trajectories/%s",args$trajectory_name))
args$test_mode <- TRUE
## END TEST

#####################
## Define settings ##
#####################

# I/O

# Options
if (isFALSE(args$denoise)) args$knn <- 0
opts$motif_annotation <- args$motif_annotation

#####################
## Load trajectory ##
#####################

trajectory.dt <- fread(args$trajectory)

###################
## Load metadata ##
###################

sample_metadata <- fread(args$metadata) %>%
  .[pass_atacQC==TRUE & pass_rnaQC==TRUE & doublet_call==FALSE]

cells <- intersect(trajectory.dt$cell,sample_metadata$cell)
trajectory.dt <- trajectory.dt[cell%in%cells] 
sample_metadata <- sample_metadata[cell%in%cells] %>% setkey(cell) %>% .[cells]

############################
## Load RNA and ATAC data ##
############################

# print("Loading data...")
# io$rna.sce <- args$sce
# source(here::here("rna_atac/load_rna_atac_single_cells.R"))

# rm(list=c("rna.sce"))
# rna_tf.sce <- rna_tf.sce[,cells]
# atac.peakMatrix.se <- atac.peakMatrix.se[,cells]

###################
## Load RNA data ##
###################

print("Loading RNA data...")

rna.sce <- load_SingleCellExperiment(
  file = args$sce, 
  cells = sample_metadata$cell, 
  normalise = TRUE, 
  remove_non_expressed_genes = FALSE
)

####################
## Load ATAC data ##
####################

source(here::here("atac/archR/load_archR_project.R"))

# Subset
ArchRProject.filt <- ArchRProject[sample_metadata$cell]

atac.peakMatrix.se <- getMatrixFromProject(ArchRProject.filt, binarize = FALSE, useMatrix = "PeakMatrix")

row_ranges.dt <- rowRanges(atac.peakMatrix.se) %>% as.data.table %>% 
  setnames("seqnames","chr") %>%
  .[,c("chr","start","end")] %>%
  .[,idx:=sprintf("%s:%s-%s",chr,start,end)]
rownames(atac.peakMatrix.se) <- row_ranges.dt$idx

# Load peak metadata
peak_metadata.dt <- fread(io$archR.peak.metadata) %>%
  .[,idx:=sprintf("%s:%s-%s",chr,start,end)]
atac.peakMatrix.se <- atac.peakMatrix.se[rownames(atac.peakMatrix.se) %in% peak_metadata.dt$idx]

###############################
## Load motifmatcher results ##
###############################

print("Load motifmatcher results...")

# source(here::here("load_motifmatchR.R"))
motifmatcher.se <- readRDS(args$motifmatcher)

# Subset peaks
stopifnot(rownames(atac.peakMatrix.se)%in%rownames(motifmatcher.se))
motifmatcher.se <- motifmatcher.se[rownames(atac.peakMatrix.se),]

################################
## Load motif2gene annotation ##
################################

source(here::here("atac/archR/load_motif_annotation.R"))
# motif2gene.dt <- fread(args$motif2gene)

################
## Filter TFs ##
################

motifs <- intersect(colnames(motifmatcher.se),motif2gene.dt$motif)
motifmatcher.se <- motifmatcher.se[,motifs]
motif2gene.dt <- motif2gene.dt[motif%in%motifs]


genes <- intersect(toupper(rownames(rna.sce)),motif2gene.dt$gene)
rna_tf.sce <- rna.sce[str_to_title(genes),]
rownames(rna_tf.sce) <- toupper(rownames(rna_tf.sce))
motif2gene.dt <- motif2gene.dt[gene%in%genes]

# Duplicated gene-motif pairs
# sum(duplicated(motif2gene.dt$motif))
# sum(duplicated(motif2gene.dt$gene))
# motif2gene.dt[,.N,by="motif"] %>% .[N>1]
# motif2gene.dt[,.N,by="gene"] %>% .[N>1]

#################
## Filter data ##
#################

print("Filtering data...")

# RNA 
rna_tf.sce <- rna_tf.sce[Matrix::rowSums(counts(rna_tf.sce))>100]

# ATAC
atac.peakMatrix.se <- atac.peakMatrix.se[Matrix::rowSums(assay(atac.peakMatrix.se))>100]

motifmatcher.se <- motifmatcher.se[rownames(atac.peakMatrix.se),]

#################
## Smooth data ##
#################

# stopifnot(colnames(rna_tf.sce) == trajectory.dt$cell)

if (args$denoise & args$knn>1) {
  
  print("kNN denoising...")
  
  # RNA
  pca.rna <- fread(args$pca_rna) %>% matrix.please %>% .[sample_metadata$cell,]
  rna_tf.mtx <- smoother_aggregate_nearest_nb(mat=as.matrix(logcounts(rna_tf.sce)), D=pdist(pca.rna), k=args$knn)
  colnames(rna_tf.mtx) <- colnames(rna_tf.sce)
  
  # ATAC peaks
  pca.atac <- fread(args$pca_atac) %>% matrix.please %>% .[sample_metadata$cell,]
  atac_peak.mtx <- smoother_aggregate_nearest_nb(mat=as.matrix(assay(atac.peakMatrix.se)), D=pdist(pca.atac), k=args$knn)
  colnames(atac_peak.mtx) <- colnames(atac.peakMatrix.se)
  
} else {
  
  rna_tf.mtx <- as.matrix(logcounts(rna_tf.sce))
  atac_peak.mtx <- as.matrix(assay(atac.peakMatrix.se))
  
}

# Make sure that cells have the same order to calculate the correlations
rna_tf.mtx <- rna_tf.mtx[,sample_metadata$cell]
atac_peak.mtx <- atac_peak.mtx[,sample_metadata$cell]
stopifnot(colnames(rna_tf.mtx)==colnames(atac_peak.mtx))

##############################################################
## Correlate ATAC peak accessibility with TF RNA expression ##
##############################################################

print("Correlating ATAC peak accessibility with TF RNA expression...")

TFs <- rownames(rna_tf.mtx)
if (args$test_mode) {
  print("Test mode activated...")
  TFs <- TFs %>% head(n=5)
}

# Prepare output data objects
cor.mtx <- matrix(as.numeric(NA), nrow=nrow(atac_peak.mtx), ncol=length(TFs))
pvalue.mtx <- matrix(as.numeric(NA), nrow=nrow(atac_peak.mtx), ncol=length(TFs))
rownames(cor.mtx) <- rownames(atac_peak.mtx); colnames(cor.mtx) <- TFs
dimnames(pvalue.mtx) <- dimnames(cor.mtx)

for (i in TFs) {
  print(i)

  motif_i <- motif2gene.dt[gene==i,motif]

  all_peaks_i <- rownames(motifmatcher.se)[which(assay(motifmatcher.se[,motif_i],"motifMatches")==1)]

  # calculate correlations
  corr_output <- psych::corr.test(rna_tf.mtx[i,], t(atac_peak.mtx[all_peaks_i,]), ci=FALSE)
  
  # Fill matrices
  cor.mtx[all_peaks_i,i] <- round(corr_output$r[1,],3)
  pvalue.mtx[all_peaks_i,i] <- round(corr_output$p[1,],5)
  
}


##########
## Save ##
##########

print("Saving output...")

to.save <- SummarizedExperiment(
  assays = SimpleList("cor" = dropNA(cor.mtx), "pvalue" = dropNA(pvalue.mtx)),
  rowData = rowData(atac.peakMatrix.se)
)

if (args$denoise) {
  saveRDS(to.save, file.path(args$outdir,sprintf("TFexpr_vs_peakAcc_%s_%s_denoise_knn%s.rds",args$motif_annotation,args$trajectory_name,args$knn)))
} else {
  saveRDS(to.save, file.path(args$outdir,sprintf("TFexpr_vs_peakAcc_%s_%s.rds",args$motif_annotation,args$trajectory_name)))
}

##########
## Test ##
##########

# i <- "GATA1"
# j <- colnames(corr_output$r)[which.max(corr_output$r)]
# 
# to.plot <- data.table(
#   rna = rna_tf.mtx[i,],
#   atac = atac_peak.mtx[j,],
#   cell = colnames(rna_tf.mtx)
# ) %>%
#   merge(sample_metadata[,c("cell","celltype.predicted")]) %>%
#   merge(trajectory.dt[,c("cell","PC1")])
# 
# ggscatter(to.plot, x="rna", y="atac", fill="celltype.predicted", size=1, shape=21,
#                add="reg.line", add.params = list(color="black", fill="lightgray"), conf.int=TRUE) +
#   stat_cor(method = "pearson") +
#   scale_fill_manual(values=opts$celltype.colors) +
#   labs(x=sprintf("%s expression",i), y="chromatin accessibility", title=j) +
#   guides(fill="none", color="none") +
#   theme(
#     plot.title = element_text(hjust = 0.5, size=rel(0.85)),
#     axis.text = element_text(size=rel(0.7))
#   )
# 
# ggscatter(to.plot[,.(rna=mean(rna), atac=mean(atac)), by="celltype.predicted"], x="rna", y="atac", fill="celltype.predicted", size=4, shape=21,
#           add="reg.line", add.params = list(color="black", fill="lightgray"), conf.int=TRUE) +
#   stat_cor(method = "pearson") +
#   scale_fill_manual(values=opts$celltype.colors) +
#   labs(x=sprintf("%s expression",i), y="chromatin accessibility", title=j) +
#   guides(fill="none", color="none") +
#   theme(
#     plot.title = element_text(hjust = 0.5, size=rel(0.85)),
#     axis.text = element_text(size=rel(0.7))
#   )
