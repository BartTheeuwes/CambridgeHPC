###################
## Sanity checks ##
###################


# Check that variables exist
# opts$motif_annotation
# io$pca.atac
# io$pca.rna
# ArchRproject.filt
stopifnot(!is.null(opts$motif_annotation))

###############################
## Load SingleCellExperiment ##
###############################

rna.sce <- load_SingleCellExperiment(
  file = io$rna.sce, 
  cells = sample_metadata$cell, 
  normalise = TRUE, 
  remove_non_expressed_genes = FALSE
)

#####################
## Load peakMatrix ##
#####################

if (file.exists(io$archR.peakMatrix.se)) {
  cat("Loading precomputed ATAC peak matrix...")
  atac.peakMatrix.se <- readRDS(io$archR.peakMatrix.se) %>% .[,colnames(.) %in%sample_metadata$cell]
} else {
  
  if (exi(ArchRProject))
  cat("Precomputed ATAC peak matrix not found, fetching from the ArchR object...")
  atac.peakMatrix.se <- getMatrixFromProject(ArchRProject, binarize = TRUE, useMatrix = "PeakMatrix") %>% .[,colnames(.) %in%sample_metadata$cell]
  
  row_ranges.dt <- rowRanges(atac.peakMatrix.se) %>% as.data.table %>% 
    setnames("seqnames","chr") %>%
    .[,c("chr","start","end")] %>%
    .[,idx:=sprintf("%s:%s-%s",chr,start,end)]
  rownames(atac.peakMatrix.se) <- row_ranges.dt$idx
  
  # Load peak metadata
  peak_metadata.dt <- fread(io$archR.peak.metadata) %>%
    .[,idx:=sprintf("%s:%s-%s",chr,start,end)]
  atac.peakMatrix.se <- atac.peakMatrix.se[rownames(atac.peakMatrix.se) %in% peak_metadata.dt$idx]
  
}

##########################
## Load chromVAR scores ##
##########################

# io$archr.chromvar.se <- sprintf("%s/chromVAR_deviations_summarized_experiment_%s_correlated_peaks_archr.rds",io$archr.chromvar.dir,opts$motif_annotation)
# opts$motif_annotation <- "Motif_cisbp_lenient"
# io$archr.chromvar.se <- sprintf("%s/results/atac/archR/chromvar/chromRNA_deviations_%s_archr.rds",io$basedir,opts$motif_annotation)
# if (file.exists(io$archr.chromvar.se)) {
#   cat("Loading precomputed chromVAR matrix...")
#   atac.chromvar.se <- readRDS(io$archr.chromvar.se) %>% .[,colnames(.) %in%sample_metadata$cell]
# } else {
#   stop("Precomputed chromVAR matrix not found")
# }

################################
## Load motif2gene annotation ##
################################

# source(here::here("atac/archR/load_motif_annotation.R"))
# 
# motif2gene.dt <- motif2gene.dt %>%
#   .[gene%in%rownames(atac.chromvar.se) & gene%in%toupper(rownames(rna.sce))] %>%
#   .[,N:=length(unique(motif)),by="gene"] %>% .[N==1] %>% .[,N:=NULL]
# 
# rna_tf.sce <- rna.sce[str_to_title(motif2gene.dt$gene),]
# rownames(rna_tf.sce) <- toupper(rownames(rna_tf.sce))
# 
# atac.chromvar.se <- atac.chromvar.se[motif2gene.dt$gene,]

#############################
## Create long data.tables ##
#############################

# chromvar.dt <- atac.chromvar.mtx %>% t %>% 
#   as.data.table(keep.rownames = T) %>%
#   setnames("rn","cell") %>% 
#   melt(id.vars=c("cell"), variable.name="gene", value.name="chromvar_zscore")

# rna.dt <- logcounts(rna.sce) %>%
#   as.data.table(keep.rownames = T) %>%
#   setnames("rn","gene") %>%
#   melt(id.vars="gene", variable.name="cell", value.name="expr")

# rna_tf.dt <- logcounts(rna_tf.sce) %>%
#   as.data.table(keep.rownames = T) %>%
#   setnames("rn","gene") %>%
#   melt(id.vars="gene", variable.name="cell", value.name="expr")
