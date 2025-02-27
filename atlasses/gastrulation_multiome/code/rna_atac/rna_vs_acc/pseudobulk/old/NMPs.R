
#####################
## Define settings ##
#####################

# Load default settings
if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
  source("/Users/ricard/gastrulation_multiome_10x/utils.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
  source("/homes/ricard/gastrulation_multiome_10x/utils.R")
} else {
  stop("Computer not recognised")
}

# Options
opts$motif_annotation <- "Motif_cisbp"

opts$celltypes = c(
  "Epiblast",
  "Primitive_Streak",
  "Caudal_epiblast",
  # "Caudal_Mesoderm",
  "Somitic_mesoderm",
  "NMP",
  "Spinal_cord"
)

# I/O
# io$outdir <- paste0(io$basedir,"/results/rna_atac/rna_vs_acc/pseudobulk/TFexpr_vs_peakAcc")
io$archR.pseudobulk.deviations.se <- sprintf("%s/pseudobulk/pseudobulk_DeviationMatrix_%s_summarized_experiment.rds",io$archR.directory,opts$motif_annotation)

#######################################
## Load pseudobulk RNA and ATAC data ##
#######################################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/rna_atac/rna_vs_acc/pseudobulk/load_rna_atac_pseudobulk.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/rna_atac/rna_vs_acc/pseudobulk/load_rna_atac_pseudobulk.R")
} else {
  stop("Computer not recognised")
}

rna.sce <- rna.sce[,opts$celltypes]
rna.sce.tf <- rna.sce.tf[,opts$celltypes]
atac.chromvar.se <- atac.chromvar.se[,opts$celltypes]
atac.peakMatrix.se <- atac.peakMatrix.se[,opts$celltypes]

###############################
## Load motifmatcher results ##
###############################

source("/Users/ricard/gastrulation_multiome_10x/load_motifmatchR.R")

#######################
## Load TF2peak data ##
#######################

io$file <- paste0(io$basedir,"/results/rna_atac/rna_vs_acc/pseudobulk/TFexpr_vs_peakAcc/cor_TFexpr_vs_peakAcc_SummarizedExperiment.rds")
tf2peak_cor.se <- readRDS(io$file)

io$file <- paste0(io$basedir,"/results/rna_atac/rna_vs_acc/pseudobulk/TFexpr_vs_peakAcc/cor_TFexpr_vs_peakAcc.txt.gz")
tf2peak_cor.dt <- fread(io$file) %>%
  .[!is.na(cor)] %>%
  .[,cor_sign:=c("-","+")[(cor>0)+1]]

#######################
## Load marker genes ##
#######################

marker_genes.dt <- fread(io$rna.atlas.marker_genes) %>%
  .[celltype%in%opts$celltypes]
table(marker_genes.dt$celltype)

#######################
## Load TF2gene data ##
#######################

io$file <- paste0(io$basedir,"/results/rna_atac/rna_vs_acc/pseudobulk/TFexpr_vs_Geneexpr/cor_TFexpr_vs_Geneexpr_SummarizedExperiment_mincor0.25.rds")
tf2gene_cor.se <- readRDS(io$file)


##########
## Test ##
##########

# - Identify TFs for each marker gene
# - Differenetial accessibility between Spinal cord and Caudal mesoderm
# - Differential expression between Spinal cord and caudal mesoderm -> Identify genes with T/Sox2 motif
# rna.diff <- fread(paste0(io$rna.differential,"/Caudal_epiblast_vs_NMP.txt.gz")) %>%
rna.diff <- fread(paste0(io$rna.differential,"/NMP_vs_Spinal_cord.txt.gz")) %>%
# rna.diff <- fread(paste0(io$rna.differential,"/Somitic_mesoderm_vs_Spinal_cord.txt.gz")) %>%
  .[sig==T]
# - Characterisation of the TF network that determines NMP identity

# atac.diff <- fread(paste0(io$archR.peak.differential.dir,"/PeakMatrix_Somitic_mesoderm_vs_Spinal_cord.txt.gz")) %>%
atac.diff <- fread(paste0(io$archR.peak.differential.dir,"/PeakMatrix_Caudal_epiblast_vs_NMP.txt.gz")) %>%
  .[abs(MeanDiff)>0.10 & FDR<0.01] %>% setorder(FDR)

marker_genes.dt[celltype=="NMP"]
i <- "Cdx4"
stopifnot(i%in%colnames(tf2gene_cor.se))
which(assay(tf2gene_cor.se[,i],"cor")[,1] > 0.5)

T.genes <- which(abs(assay(tf2gene_cor.se["T",],"cor")[1,]) > 0.5) %>% names
Sox2.genes <- which(abs(assay(tf2gene_cor.se["SOX2",],"cor")[1,]) > 0.5) %>% names
intersect(T.genes,Sox2.genes)

foo <- which(abs(assay(tf2gene_cor.se[,rna.diff[logFC>0,gene]],"cor")[1,]) > 0.5) %>% names

T.peaks <- which(assay(tf2peak_cor.se[,"T"],"cor")[,1] > 0.5) %>% names

tf2peak_cor.dt[TF=="T"] %>% View

logcounts(rna.sce.tf[,"NMP"])[,1] %>% sort
