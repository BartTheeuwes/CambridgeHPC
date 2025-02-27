suppressPackageStartupMessages(library(argparse))

here::i_am("rna_atac/rna_vs_acc/pseudobulk/TFexpr_vs_peakAcc/run_cor_TFexpr_vs_peakAcc.R")

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--sce',  type="character",              help='RNA SingleCellExperiment (pseudobulk)') 
p$add_argument('--atac_peak_matrix',  type="character",              help='ATAC Peak matrix (pseudobulk)') 
# p$add_argument('--motif2gene',  type="character",              help='Motif annotation') 
p$add_argument('--motifmatcher',  type="character",              help='Motif annotation') 
p$add_argument('--motif_annotation',  type="character",              help='Motif annotation') 
p$add_argument('--outdir',          type="character",                help='Output directory')
p$add_argument('--test',  action="store_true",  help='Test mode?')
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

source(here::here("settings.R"))
source(here::here("utils.R"))

## START TEST ##
# args <- list()
# args$sce <- file.path(io$basedir,"results_new/rna/pseudobulk/SingleCellExperiment_pseudobulk_celltype.mapped_mnn.rds") # io$rna.pseudobulk.sce
# args$atac_peak_matrix <- file.path(io$basedir,"results_new/atac/archR/pseudobulk/celltype.mapped_mnn/pseudobulk_PeakMatrix_summarized_experiment.rds")
# args$outdir <- file.path(io$basedir,"results_new/rna_atac/rna_vs_acc/pseudobulk/TFexpr_vs_peakAcc")
# args$motif_annotation <- "JASPAR"
# args$motifmatcher <- sprintf("%s/Annotations/%s-Scores.rds",io$archR.directory,args$motif_annotation)
# args$motif2gene <- sprintf("%s/Annotations/%s_TFs.txt.gz",io$archR.directory,args$motif_annotation)
# args$test <- TRUE
## END TEST ##

# I/O
dir.create(args$outdir, showWarnings = F)

# Options
opts$motif_annotation <- args$motif_annotation

##################################
## Load pseudobulk RNA and ATAC ##
##################################

# Load SingleCellExperiment
rna_pseudobulk.sce <- readRDS(args$sce)

# Load ATAC SummarizedExperiment
atac_pseudobulk_PeakMatrix.se <- readRDS(args$atac_peak_matrix)

# Rename features
# if (any(grepl("f+",rownames(atac_pseudobulk_PeakMatrix.se)))) {
#   peak_names <- rowData(atac_pseudobulk_PeakMatrix.se) %>% as.data.table %>% .[,idx:=sprintf("%s:%s-%s",seqnames,start,end)] %>% .$id
#   rownames(atac_pseudobulk_PeakMatrix.se) <- peak_names
# }

# Make sure that samples are consistent
celltypes <- intersect(colnames(rna_pseudobulk.sce),colnames(atac_pseudobulk_PeakMatrix.se))
rna_pseudobulk.sce <- rna_pseudobulk.sce[,celltypes]
atac_pseudobulk_PeakMatrix.se <- atac_pseudobulk_PeakMatrix.se[,celltypes]

###############################
## Load motifmatcher results ##
###############################

# source(here::here("load_motifmatchR.R"))
motifmatcher.se <- readRDS(args$motifmatcher)

# Subset peaks
stopifnot(sort(rownames(motifmatcher.se))==sort(rownames(atac_pseudobulk_PeakMatrix.se)))
motifmatcher.se <- motifmatcher.se[rownames(atac_pseudobulk_PeakMatrix.se),]

###########################
## Load motif annotation ##
###########################

source(here::here("atac/archR/load_motif_annotation.R"))
# motif2gene.dt <- fread(args$motif2gene)

################
## Filter TFs ##
################

motifs <- intersect(colnames(motifmatcher.se),motif2gene.dt$motif)
motifmatcher.se <- motifmatcher.se[,motifs]
motif2gene.dt <- motif2gene.dt[motif%in%motifs]

genes <- intersect(toupper(rownames(rna_pseudobulk.sce)),motif2gene.dt$gene)
rna_tf_pseudobulk.sce <- rna_pseudobulk.sce[str_to_title(genes),]
rownames(rna_tf_pseudobulk.sce) <- toupper(rownames(rna_tf_pseudobulk.sce))
motif2gene.dt <- motif2gene.dt[gene%in%genes]

# Duplicated gene-motif pairs
# sum(duplicated(motif2gene.dt$motif))
# sum(duplicated(motif2gene.dt$gene))
# motif2gene.dt[,.N,by="motif"] %>% .[N>1]
# motif2gene.dt[,.N,by="gene"] %>% .[N>1]

#########################################################
## Correlate peak accessibility with TF RNA expression ##
#########################################################

# Sanity checks
stopifnot(colnames(rna_tf_pseudobulk.sce)==colnames(atac_pseudobulk_PeakMatrix.se))

# TFs <- intersect(colnames(motifmatcher.se),rownames(rna_pseudobulk.sce))# %>% head(n=5)
# TFs <- c("GATA1","TAL1")
TFs <- genes

if (args$test) {
  TFs <- TFs %>% head(n=5)
}

# Prepare output data objects
cor.mtx <- matrix(as.numeric(NA), nrow=nrow(atac_pseudobulk_PeakMatrix.se), ncol=length(TFs))
pvalue.mtx <- matrix(as.numeric(NA), nrow=nrow(atac_pseudobulk_PeakMatrix.se), ncol=length(TFs))
rownames(cor.mtx) <- rownames(atac_pseudobulk_PeakMatrix.se); colnames(cor.mtx) <- TFs
dimnames(pvalue.mtx) <- dimnames(cor.mtx)

# i <- "GATA1"
for (i in TFs) {
  motif_i <- motif2gene.dt[gene==i,motif]
  
  print(i)
  all_peaks_i <- rownames(motifmatcher.se)[which(assay(motifmatcher.se[,motif_i],"motifMatches")==1)]
  
  # calculate correlations
  corr_output <- psych::corr.test(
    x = t(logcounts(rna_tf_pseudobulk.sce[i,])), 
    y = t(assay(atac_pseudobulk_PeakMatrix.se[all_peaks_i,])), 
    ci = FALSE
  )
  
  # Fill matrices
  cor.mtx[all_peaks_i,i] <- round(corr_output$r[1,],3)
  pvalue.mtx[all_peaks_i,i] <- round(corr_output$p[1,],10)
}

##########
## Save ##
##########

# Save SummarizedExperiment object
to.save <- SummarizedExperiment(
  assays = SimpleList("cor" = dropNA(cor.mtx), "pvalue" = dropNA(pvalue.mtx)),
  rowData = rowData(atac_pseudobulk_PeakMatrix.se)
)
saveRDS(to.save, file.path(args$outdir,sprintf("%s_cor_TFexpr_vs_peakAcc.rds",args$motif_annotation)))

##########
## TEST ##
##########

# target_peaks_i <- which(!is.na(dropNA2matrix(assay(to.save[,i],"pvalue"))[,1]))

