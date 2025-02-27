here::i_am("atac/archR/dimensionality_reduction/cells/atac_dimensionality_reduction_cells.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

suppressPackageStartupMessages(library(ArchR))
suppressPackageStartupMessages(library(uwot))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--metadata',        type="character",                               help='Cell metadata file')
p$add_argument('--matrix',          type="character",   help='Matrix to use')
p$add_argument('--atac_matrix_file',          type="character",  help='Matrix file')
p$add_argument('--atac_feature_stats',          type="character",   help='Feature stats')
p$add_argument('--stages',       type="character",  default="all",  nargs='+',  help='Stages to plot')
p$add_argument('--samples',       type="character",  default="all",  nargs='+',  help='Samples to plot')
p$add_argument('--celltypes',       type="character",  default="all",  nargs='+',  help='Samples to plot')
p$add_argument('--remove_ExE_cells',       type="character",  default="False",  help='Remove ExE cells? ("True"/"False")')
p$add_argument('--binarise', action="store_true",                                 help='Binarise ATAC matrix?')
p$add_argument('--scale_dims', action="store_true",                                 help='Scale latent dimensions?')
p$add_argument('--nfeatures',       type="integer",    default=1000,               help='Number of features')
p$add_argument('--ndims',           type="integer",    default=30,                  help='Number of LSI dimensions')
p$add_argument('--batch_variable',  type="character",                               help='Metadata column to apply batch correction on')
p$add_argument('--batch_method',    type="character",  default="MNN",               help='Batch correctin method ("Harmony" or "MNN")')
p$add_argument('--n_neighbors',     type="integer",    default=30,   nargs='+',     help='(UMAP) Number of neighbours')
p$add_argument('--min_dist',        type="double",     default=0.3,  nargs='+',     help='(UMAP) Minimum distance')
p$add_argument('--colour_by',       type="character",  nargs='+',  help='Metadata columns to colour the UMAP by')
p$add_argument('--seed',            type="integer",    default=42,                  help='Random seed')
p$add_argument('--outdir',          type="character",                               help='Output directory')

args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
io$basedir <- file.path(io$basedir,"test")
args <- list()
args$matrix <- "PeakMatrix"
args$metadata <- file.path(io$basedir,"results/atac/archR/celltype_assignment/sample_metadata_after_celltype_assignment.txt.gz")
args$atac_feature_stats <- file.path(io$basedir,sprintf("results/atac/archR/feature_stats/%s_celltype_stats.txt.gz",args$matrix))
args$stages <- c("E8.5","E8.75") # "all"
args$samples <- c("E8.5_CRISPR_T_KO","E8.5_CRISPR_T_WT") # "all"
args$celltypes <- "Neural_crest"
args$nfeatures <- 15000
# args$matrix <- "PeakMatrix"
args$atac_matrix_file <- file.path(io$basedir,"processed/atac/archR/Matrices/PeakMatrix_summarized_experiment.rds")
args$ndims <- 50
args$scale_dims <- TRUE
args$outdir <- file.path(io$basedir,"results/atac/archR/dimensionality_reduction/test")
## END TEST ##

#####################
## Parse arguments ##
#####################

if (args$stages[1]=="all") {
  args$stages <- opts$stages
} else {
  stopifnot(args$stages%in%opts$stages)
}

if (args$celltypes[1]=="all") {
  args$celltypes <- opts$celltypes
} else {
  stopifnot(args$celltypes%in%opts$celltypes)
}

# I/O
dir.create(args$outdir, showWarnings=F, recursive=T)

# Options
opts$remove_dim_cor_seq_depth <- TRUE

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(args$metadata) %>%
  .[pass_atacQC==TRUE & doublet_call==FALSE & sample%in%args$samples & stage%in%args$stages & celltype%in%args$celltypes]

table(sample_metadata$sample,sample_metadata$celltype)

#######################
## Fetch ATAC matrix ##
#######################

print(sprintf("Fetching single-cell ATAC %s...",args$matrix))

if (file.exists(args$atac_matrix_file)) {
  atac.se <- readRDS(args$atac_matrix_file)[,sample_metadata$cell]
} else {
  print(sprintf("%s does not exist. Loading matrix from the ArchR project...",args$atac_matrix_file))
  source(here::here("atac/archR/load_archR_project.R"))
  stopifnot(args$matrix %in% getAvailableMatrices(ArchRProject))
  ArchRProject.filt <- ArchRProject[sample_metadata$cell]
  atac.se <- getMatrixFromProject(ArchRProject.filt, useMatrix=args$matrix, binarize = FALSE)[,sample_metadata$cell]
  dim(atac.se)

  # Define feature names
  if (grepl("peak",tolower(args$matrix),ignore.case=T)) {
    row.ranges.dt <- rowRanges(atac.se) %>% as.data.table %>% 
      setnames("seqnames","chr") %>%
      .[,c("chr","start","end")] %>%
      .[,idx:=sprintf("%s:%s-%s",chr,start,end)]
    rownames(atac.se) <- row.ranges.dt$idx
  } else if (grepl("gene",tolower(args$matrix),ignore.case=T)) {
    rownames(atac.se) <- rowData(atac.se)$name
  } else {
    stop("Matrix not recognised")
  }
}

#######################
## Feature selection ##
#######################

tmp <- apply(assay(atac.se),1,var) %>% sort
atac_features <- tmp[tmp>0.5] %>% names

#############################
## ATAC TFIDF normalisation ##
#############################

atac.mtx <- assay(atac.se[atac_features,])
atac.mtx[atac.mtx>1] <- 1
atac_tfidf.mtx <- tfidf(atac.mtx, method=1, scale.factor=1e4)

###########################
## Latent Semantic Index ##
###########################

svd <- irlba::irlba(atac_tfidf.mtx, args$ndims, args$ndims)
svdDiag <- matrix(0, nrow=args$ndims, ncol=args$ndims); diag(svdDiag) <- svd$d
lsi.mtx <- t(svdDiag %*% t(svd$v))
rownames(lsi.mtx) <- colnames(atac.se)
colnames(lsi.mtx) <- paste0("LSI",seq_len(ncol(lsi.mtx)))

###############
## Parse LSI ##
###############

# (Optional) Remove dimensions that correlate with nFrags
# dims.to.keep <- which(abs(cor(lsi.mtx,sample_metadata$nFrags_atac))<=0.75)
# lsi.mtx <- lsi.mtx[,dims.to.keep]

# (Optional) scale dimensions, Z scores
# lsi.mtx <- sweep(lsi.mtx - rowMeans(lsi.mtx), 1, matrixStats::rowSds(lsi.mtx),`/`)
# if (args$scale_dims) {
#   lsi.mtx <- sweep(lsi.mtx - colMeans(lsi.mtx), 2, matrixStats::colSds(lsi.mtx),`/`)
#   lsi.mtx[lsi.mtx>2] <- 2
#   lsi.mtx[lsi.mtx<(-2)] <- (-2)
# }

##########
## UMAP ##
##########

# # Run UMAP
# set.seed(args$seed)
# umap_embedding.mtx <- umap(lsi.mtx, n_neighbors=i, min_dist=j, metric="cosine", fast_sgd = TRUE) %>% round(2)
# rownames(umap_embedding.mtx) <- rownames(lsi.mtx)
# 
# # Fetch UMAP coordinates
# umap.dt <- umap_embedding.mtx %>%
#   as.data.table(keep.rownames = T) %>%
#   setnames(c("cell","umap1","umap2"))

##################
## Find regions ##
##################

lsi_weights.mtx <- t(svdDiag %*% t(svd$u))
rownames(lsi_weights.mtx) <- rownames(atac_tfidf.mtx)
colnames(lsi_weights.mtx) <- paste0("LSI",seq_len(ncol(lsi_weights.mtx)))

# lsi_weights.mtx <- lsi_weights.mtx[,dims.to.keep]

sort(lsi_weights.mtx[,2]) %>% head(n=15)
sort(lsi_weights.mtx[,2]) %>% tail(n=15)

##########
## Plot ##
##########

to.plot <- lsi.mtx[,c(1,2,3)] %>%
  as.data.table(keep.rownames = T) %>%
  setnames(c("cell","x","y","z")) %>% 
  merge(sample_metadata[,c("cell","sample","stage","genotype","celltype")],by="cell")

ggplot(to.plot, aes_string(x="y", y="z", fill="sample")) +
  geom_point(size=2, shape=21, stroke=0.05) +
  # scale_fill_manual(values=opts$stage.colors) +
  theme_classic()


to.plot <- lsi.mtx[,c(1,2,3)] %>%
  as.data.table(keep.rownames = T) %>%
  setnames(c("cell","x","y","z")) %>% 
  merge(sample_metadata[,c("cell","sample","stage","genotype","celltype")],by="cell")

ggplot(to.plot, aes_string(x="z", y="y", fill="stage")) +
  geom_point(size=2, shape=21, stroke=0.05) +
  # scale_fill_manual(values=opts$stage.colors) +
  theme_classic()

##########
## Plot ##
##########

# features.to.plot <- c(sort(lsi_weights.mtx[,2]) %>% head(n=10), sort(lsi_weights.mtx[,2]) %>% tail(n=10)) %>% names %>% unique
features.to.plot <- c("chr16:76225043-76225643","chr6:30636298-30636898")
foo <- atac_tfidf.mtx[features.to.plot,] %>% as.matrix %>% as.data.table(keep.rownames = T) %>% 
  setnames("rn","feature") %>%
  melt(id.vars="feature", variable.name="cell") %>%
  merge(to.plot,by="cell", allow.cartesian=T)

ggplot(foo, aes_string(x="y", y="z", color="value")) +
  facet_wrap(~feature) +
  geom_point(size=0.5) +
  scale_color_gradient(low = "gray80", high = "purple") +
  theme_classic() +
  theme(
    strip.text = element_text(size=rel(0.7)),
    legend.position = "right"
  )

##########
## TEST ##
##########

# to.plot <- umap.dt %>% merge(sample_metadata,by="cell")

# # to.plot[,foo:=F]
# # to.plot[pass_rnaQC==F & nFrags_atac_log10<=4,foo:=TRUE]
# # to.plot[,foo:=nFrags_atac<=10000]

# ggplot(to.plot, aes_string(x="umap1", y="umap2", fill="PromoterRatio_atac")) +
#   geom_point(size=1.5, shape=21, stroke=0.05) +
#   scale_fill_gradientn(colours = terrain.colors(10)) +
#   theme_classic() +
#   ggplot_theme_NoAxes()
