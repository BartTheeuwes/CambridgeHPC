suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(ArchR))
suppressPackageStartupMessages(library(uwot))

here::i_am("atac/archR/dimensionality_reduction/atac_dimensionality_reduction.R")

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--metadata',        type="character",                               help='Cell metadata file')
p$add_argument('--matrix',          type="character",  default="PeakMatrix",   help='Matrix to use')
p$add_argument('--stages',       type="character",  default="all",  nargs='+',  help='Stages to plot')
p$add_argument('--remove_ExE_cells', action="store_true",                                 help='Remove ExE cells?')
p$add_argument('--binarise', action="store_true",                                 help='Binarise ATAC matrix?')
p$add_argument('--scale_dims', action="store_true",                                 help='Scale latent dimensions?')
p$add_argument('--nfeatures',       type="integer",    default=1000,               help='Number of features')
p$add_argument('--ndims',           type="integer",    default=30,                  help='Number of LSI dimensions')
p$add_argument('--batch.variable',  type="character",                               help='Metadata column to apply batch correction on')
p$add_argument('--batch.method',    type="character",  default="MNN",               help='Batch correctin method ("Harmony" or "MNN")')
p$add_argument('--n_neighbors',     type="integer",    default=30,   nargs='+',     help='(UMAP) Number of neighbours')
p$add_argument('--min_dist',        type="double",     default=0.3,  nargs='+',     help='(UMAP) Minimum distance')
p$add_argument('--colour_by',       type="character",  nargs='+',  help='Metadata columns to colour the UMAP by')
p$add_argument('--seed',            type="integer",    default=42,                  help='Random seed')
p$add_argument('--outdir',          type="character",                               help='Output directory')

args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
args <- list()
args$metadata <- io$metadata #file.path(io$basedir,"results_new/atac/archR/qc/sample_metadata_after_qc.txt.gz")
args$stages <- "all" # "E8.75"
args$remove_ExE_cells <- TRUE
args$nfeatures <- 25000
args$matrix <- "PeakMatrix"
args$binarise <- FALSE
args$ndims <- 50
args$scale_dims <- FALSE
args$seed <- 42
args$n_neighbors <- 25
args$min_dist <- 0.50
args$colour_by <- c("stage","celltype.mapped_mnn","nFrags_atac")
args$outdir <- file.path(io$basedir,"results_new/atac/archR/dimensionality_reduction/test")
args$test <- FALSE
## END TEST ##


#####################
## Define settings ##
#####################

source(here::here("settings.R"))
source(here::here("utils.R"))

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
    .[!celltype.mapped_mnn%in%c("Visceral_endoderm","ExE_endoderm","ExE_ectoderm","Parietal_endoderm")]
}

if (args$test) {
  print("Test mode activated, subsetting number of cells...")
  sample_metadata <- sample_metadata %>% head(n=250)
}

table(sample_metadata$stage)
table(sample_metadata$celltype.mapped_mnn)

########################
## Load ArchR project ##
########################

source(here::here("atac/archR/load_archR_project.R"))

ArchRProject.filt <- ArchRProject[sample_metadata$cell]

###################
## Sanity checks ##
###################

stopifnot(args$matrix %in% getAvailableMatrices(ArchRProject))

if (length(args$batch.variable)>0) {
  stopifnot(args$batch.variable%in%colnames(sample_metadata))
  if (length(unique(sample_metadata[[args$batch.variable]]))==1) {
    message(sprintf("There is a single level for %s, no batch correction applied",args$batch.variable))
    args$batch.variable <- NULL
  } else {
    library(batchelor)
  }
}

#######################
## Fetch ATAC matrix ##
#######################

atac.se <- getMatrixFromProject(ArchRProject.filt, useMatrix=args$matrix, binarize = args$binarise)
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

#######################
## Feature selection ##
#######################

# Load feature stats
atac_featureStats.dt <- fread(file.path(io$basedir, sprintf("results_new/atac/archR/feature_stats/%s_stats.txt.gz",args$matrix)))

# Define highly variable features using the pseudobulk estimates
atac_features <- atac_featureStats.dt %>% 
  setorder(-var_pseudobulk) %>% 
  head(n=args$nfeatures) %>% .$feature

if (args$test) {
  print("Test mode activated, subsetting number of ATAC features...")
  atac_features <- atac_features %>% head(n=500)
}

#############################
## ATAC data normalisation ##
#############################

# TFIDF normalisation
atac_tfidf.mtx <- tfidf(assay(atac.se[atac_features,]), method=1, scale.factor=1e4)

###########################
## Latent Semantic Index ##
###########################

svd <- irlba::irlba(atac_tfidf.mtx, args$ndims, args$ndims)
svdDiag <- matrix(0, nrow=args$ndims, ncol=args$ndims)
diag(svdDiag) <- svd$d
svd.mtx <- t(svdDiag %*% t(svd$v))
rownames(svd.mtx) <- colnames(atac.se)
colnames(svd.mtx) <- paste0("LSI",seq_len(ncol(svd.mtx)))

rm(svd,svdDiag)

###############
## Parse LSI ##
###############

# Remove outliers
# outlierQuantiles = c(0.02, 0.98)

# (Optional) scale dimensions, Z scores
# svd.mtx <- sweep(svd.mtx - rowMeans(svd.mtx), 1, matrixStats::rowSds(svd.mtx),`/`)
if (args$scale_dims) {
  svd.mtx <- sweep(svd.mtx - colMeans(svd.mtx), 2, matrixStats::colSds(svd.mtx),`/`)
  svd.mtx[svd.mtx>2] <- 2
  svd.mtx[svd.mtx<(-2)] <- (-2)
}


######################
## Batch correction ##
######################

if (length(args$batch.variable)>0) {
  print(sprintf("Applying %s batch correction for variable: %s", args$batch.method, args$batch.variable))
  outfile <- sprintf("%s/lsi_%s_nfeatures%d_dims%d_%sbatchcorrection_by_%s.txt.gz",args$outdir, args$matrix, args$nfeatures, args$ndims, args$batch.method, paste(args$batch.variable,collapse="-"))
  
  # Harmony
  if (args$batch.method=="Harmony") {
    
    # (...)
    # library(harmony)
    # harmonyParams <- list(...)
    # harmonyParams$data_mat <- getReducedDims(
    #   ArchRProj = ArchRProj, 
    #   reducedDims = reducedDims, 
    #   dimsToUse = dimsToUse, 
    #   scaleDims = scaleDims, 
    #   corCutOff = corCutOff
    # )
    # harmonyParams$verbose <- verbose
    # harmonyParams$meta_data <- data.frame(getCellColData(
    #   ArchRProj = ArchRProj, 
    #   select = groupBy)[rownames(harmonyParams$data_mat), , drop = FALSE])
    # harmonyParams$do_pca <- FALSE
    # harmonyParams$vars_use <- groupBy
    # harmonyParams$plot_convergence <- FALSE
    
    lsi.dt <- getReducedDims(ArchRProject.filt, "IterativeLSI_Harmony") %>% round(3) %>% 
      as.data.table(keep.rownames = T) %>% setnames("rn","cell")

  } else if (args$batch.method=="MNN") {
    # library(batchelor)
    # Z <- reducedMNN(svd.mtx, batch=MOFAobject@samples_metadata$stage)$corrected
    stop("Not implemented")
    
  } else {
    stop("Batch correction method not recognised")
  }

} else {
  outfile <- sprintf("%s/lsi_%s_nfeatures%d_ndims%d.txt.gz",args$outdir, args$matrix, args$nfeatures, args$ndims)
  lsi.dt <- svd.mtx%>% round(3) %>% 
    as.data.table(keep.rownames = T) %>% setnames("rn","cell")
}

# Save LSI coordinates
fwrite(lsi.dt, outfile)

##########
## UMAP ##
##########

# i <- args$n_neighbors[1]
# j <- args$min_dist[1]
for (i in args$n_neighbors) {
  for (j in args$min_dist) {
    
    # Define the latent space to run UMAP on
    if (length(opts$batch.correction)>0) {
      if (args$batch.method=="Harmony") {
        dimred <- "IterativeLSI_Harmony"
      } else if  (args$batch.method=="MNN") {
        dimred <- "IterativeLSI_MNN"
      } else {
        stop("Batch correction method not recognised")
      }
    } else {
      dimred <- "IterativeLSI"
    }
    
    # Run UMAP
    set.seed(args$seed)
    umap_embedding.mtx <- umap(svd.mtx, n_neighbors=i, min_dist=j, metric="cosine") %>% round(2)
    rownames(umap_embedding.mtx) <- rownames(svd.mtx)
    
    # Fetch UMAP coordinates
    umap.dt <- umap_embedding.mtx %>%
      as.data.table(keep.rownames = T) %>%
      setnames(c("cell","umap1","umap2"))
    
    # Save UMAP coordinates
    # outfile <- sprintf("%s/umap_%s_nfeatures%d_ndims%d_neigh%d_dist%s.txt.gz",args$outdir, args$matrix, args$nfeatures, args$ndims, i, j)
    outfile <- sprintf("%s/umap_%s_nfeatures%d_ndims%d.txt.gz",args$outdir, args$matrix, args$nfeatures, args$ndims)
    fwrite(umap.dt, outfile)

    # Plot
    to.plot <- umap.dt %>%
      merge(sample_metadata,by="cell")
    
    # k <- "celltype.mapped_mnn"
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
      outfile <- sprintf("%s/umap_%s_nfeatures%d_ndims%d_%s.pdf",args$outdir, args$matrix, args$nfeatures, args$ndims, k)
      pdf(outfile, width=7, height=5)
      print(p)
      dev.off()
    }
    
  }
}

