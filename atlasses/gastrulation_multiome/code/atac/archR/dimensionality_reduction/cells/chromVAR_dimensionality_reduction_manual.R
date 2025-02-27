suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(uwot))

here::i_am("atac/archR/dimensionality_reduction/cells/chromVAR_dimensionality_reduction_manual.R")

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--metadata',        type="character",                               help='Cell metadata file')
p$add_argument('--chromvar_se',          type="character",   help='ATAC chromVAR SummarizedExperiment')
p$add_argument('--stages',       type="character",  default="all",  nargs='+',  help='Stages to plot')
p$add_argument('--remove_ExE_cells', action="store_true",                                 help='Remove ExE cells?')
p$add_argument('--scale_dims', action="store_true",                                 help='Scale latent dimensions?')
p$add_argument('--ndims',           type="integer",    default=30,                  help='Number of LSI dimensions')
p$add_argument('--n_neighbors',     type="integer",    default=30,   nargs='+',     help='(UMAP) Number of neighbours')
p$add_argument('--min_dist',        type="double",     default=0.3,  nargs='+',     help='(UMAP) Minimum distance')
p$add_argument('--colour_by',       type="character",  nargs='+',  help='Metadata columns to colour the UMAP by')
p$add_argument('--seed',            type="integer",    default=42,                  help='Random seed')
p$add_argument('--outdir',          type="character",                               help='Output directory')

args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

source(here::here("settings.R"))
source(here::here("utils.R"))

## START TEST ##
args <- list()
args$metadata <- io$metadata #file.path(io$basedir,"results_new/atac/archR/qc/sample_metadata_after_qc.txt.gz")
args$stages <- "all" # 
args$remove_ExE_cells <- TRUE
args$chromvar_se <- file.path(io$basedir, "results_new/atac/archR/chromvar_chip/chromVAR_chip_deviations_CISBP_archr.rds")
# args$chromvar_se <- file.path(io$basedir, "results_new/atac/archR/chromvar/chromVAR_deviations_CISBP_archr.rds")
args$ndims <- 50
args$scale_dims <- FALSE
args$seed <- 42
args$n_neighbors <- 25
args$min_dist <- 0.50
args$colour_by <- c("stage","celltype")
args$outdir <- file.path(io$basedir,"results_new/atac/archR/dimensionality_reduction/test")
args$test <- FALSE
## END TEST ##

# I/O
dir.create(args$outdir, showWarnings=F)

# Options
if (args$stages[1]=="all") {
  args$stages <- opts$stages
} else {
  stopifnot(args$stages%in%opts$stages)
}

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(args$metadata) %>%
  # .[pass_atacQC==TRUE & doublet_call==FALSE & sample%in%args$samples] %>%
  .[nFrags_atac>=3500 & TSSEnrichment_atac>=9,pass_atacQC:=TRUE] %>% .[stage=="E8.7",stage:="E8.75"] %>%    # temporary
  .[pass_atacQC==TRUE & stage%in%args$stages] %>%
  .[,log_nFrags_atac:=log10(nFrags_atac)]

stopifnot(args$colour_by %in% colnames(sample_metadata))

if (args$remove_ExE_cells) {
  print("Removing ExE cells...")
  sample_metadata <- sample_metadata %>%
    .[!celltype%in%c("Visceral_endoderm","ExE_endoderm","ExE_ectoderm","Parietal_endoderm")]
}

if (args$test) {
  print("Test mode activated, subsetting number of cells...")
  sample_metadata <- sample_metadata %>% head(n=250)
}

table(sample_metadata$stage)
table(sample_metadata$celltype)

###########################
## Fetch chromVAR matrix ##
###########################

chromvar.se <- readRDS(args$chromvar_se)

chromvar.se <- chromvar.se[,colnames(chromvar.se)%in%sample_metadata$cell]

#######################
## Feature selection ##
#######################

# Feature selection
atac.features <- which(apply(assay(chromvar.se,"z"),1,var,na.rm=T)>0.01) %>% names

# Transformation
chromvar.mtx <- assay(chromvar.se,"z")[atac.features,] %>% as.matrix

# Mask outliers
chromvar.mtx[abs(chromvar.mtx)>=10] <- 10
hist(chromvar.mtx, breaks=50)

# Impute NAs
chromvar.mtx[is.na(chromvar.mtx)] <- 0

#########
## PCA ##
#########

pca <- irlba::prcomp_irlba(t(chromvar.mtx), n=args$ndims)$x
rownames(pca) <- colnames(chromvar.mtx)

###############
## Parse PCA ##
###############

# Remove outliers
# outlierQuantiles = c(0.02, 0.98)

# (Optional) scale dimensions, Z scores
# svd.mtx <- sweep(svd.mtx - rowMeans(svd.mtx), 1, matrixStats::rowSds(svd.mtx),`/`)
# if (args$scale_dims) {
#   svd.mtx <- sweep(svd.mtx - colMeans(svd.mtx), 2, matrixStats::colSds(svd.mtx),`/`)
#   svd.mtx[svd.mtx>2] <- 2
#   svd.mtx[svd.mtx<(-2)] <- (-2)
# }

##########
## UMAP ##
##########

# i <- args$n_neighbors[1]; j <- args$min_dist[1]
for (i in args$n_neighbors) {
  for (j in args$min_dist) {
    
    # Run UMAP
    set.seed(args$seed)
    umap_embedding.mtx <- umap(pca, n_neighbors=i, min_dist=j, metric="cosine") %>% round(2)
    rownames(umap_embedding.mtx) <- rownames(pca)
    
    # Fetch UMAP coordinates
    umap.dt <- umap_embedding.mtx %>%
      as.data.table(keep.rownames = T) %>%
      setnames(c("cell","umap1","umap2"))
    
    # Save UMAP coordinates
    # outfile <- sprintf("%s/umap_%s_nfeatures%d_ndims%d.txt.gz",args$outdir, args$matrix, args$nfeatures, args$ndims)
    # fwrite(umap.dt, outfile)

    # Plot
    to.plot <- umap.dt %>%
      merge(sample_metadata,by="cell")
    
    # k <- "celltype"
    for (k in args$colour_by) {

      # log10 large numeric values
      if (is.numeric(to.plot[[k]])) {
        if (max(to.plot[[k]],na.rm=T) - min(to.plot[[k]],na.rm=T) > 1000) {
          to.plot[[k]] <- log10(to.plot[[k]]+1)
          to.plot %>% setnames(k,paste0(k,"_log10")); k <- paste0(k,"_log10")
        }
      }
      
      p <- ggplot(to.plot, aes_string(x="umap1", y="umap2", fill=k)) +
        geom_point(size=1.5, shape=21, stroke=0.05) +
        # ggrastr::geom_point_rast(size=1.5, shape=21, stroke=0.05) +  # DOES NOT WORK IN THE CLUSTER
        theme_classic() +
        ggplot_theme_NoAxes()
      
      # Define colormap
      if (is.numeric(to.plot[[j]])) {
        p <- p + scale_fill_gradientn(colours = terrain.colors(10))
      }

      if (grepl("celltype",k)) {
        p <- p + scale_fill_manual(values=opts$celltype.colors) +
          theme(
            legend.position="none",
            legend.title=element_blank()
          )
      }
      if (grepl("stage",i)) {
        p <- p + scale_fill_manual(values=opts$stage.colors) +
          theme(
            legend.position="none",
            legend.title=element_blank()
          )
      }
      
      # Save UMAP plot
      # outfile <- sprintf("%s/umap_%s_nfeatures%d_ndims%d_neigh%d_dist%s_%s.pdf",args$outdir, args$matrix, args$nfeatures, args$ndims, i, j, k)
      outfile <- sprintf("%s/umap_chromvar_ndims%d_%s.pdf",args$outdir, args$ndims, k)
      pdf(outfile, width=7, height=5)
      print(p)
      dev.off()
    }
    
  }
}

