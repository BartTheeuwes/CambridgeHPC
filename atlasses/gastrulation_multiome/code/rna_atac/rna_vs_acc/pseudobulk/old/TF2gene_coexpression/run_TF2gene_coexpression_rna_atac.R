
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
opts$celltypes = c(
  "Epiblast",
  "Primitive_Streak",
  "Caudal_epiblast",
  "PGC",
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

# opts$motif_annotation <- "Motif_JASPAR2020_human"
opts$motif_annotation <- "Motif_cisbp"


# I/O
io$outdir <- paste0(io$basedir,"/results/rna_atac/rna_vs_acc/pseudobulk/TFexpr_vs_Geneexpr"); dir.create(io$outdir, showWarnings = F)
io$motifmatcher.se <- sprintf("%s/Annotations/%s-Matches-In-Peaks.rds",io$archR.directory,opts$motif_annotation)

##################################
## Load pseudobulk RNA and ATAC ##
##################################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/rna_atac/rna_vs_acc/pseudobulk/load_rna_atac_pseudobulk.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/rna_atac/rna_vs_acc/pseudobulk/load_rna_atac_pseudobulk.R")
} else {
  stop("Computer not recognised")
}

peak_metadata.dt %>% setkey("idx")

# Filter genes with no variability
rna.sce <- rna.sce[apply(logcounts(rna.sce),1,var)>0.01,]

###############################
## Load gene2peak associations ##
###############################

# peak2gene.dt <- fread(io$archR.peak2gene.all) %>% 
#   .[,peak:=sprintf("chr%s:%s-%s",chr,peak.start,peak.end)]

peak2gene_nearest.dt <- fread(io$archR.peak2gene.nearest) %>% 
  .[,peak:=sprintf("chr%s:%s-%s",chr,peak.start,peak.end)] %>%
  setkey(peak)

###############################
## Load motifmatcher results ##
###############################

source("/Users/ricard/gastrulation_multiome_10x/load_motif_annotation.R")

# Subset peaks
motifmatcher.se <- motifmatcher.se[rownames(atac.peakMatrix.se),]

###############################################################################
## Load correlation results between TF RNA expression and peak accessibility ##
###############################################################################


io$file <- paste0(io$basedir,"/results/rna_atac/rna_vs_acc/pseudobulk/TFexpr_vs_peakAcc/cor_TFexpr_vs_peakAcc_SummarizedExperiment.rds")
tf2peak_cor.se <- readRDS(io$file)

TFs <- colnames(tf2peak_cor.se)

############################################################################################
## Calculate correlation results between TF RNA expression and target gene RNA expression ##
############################################################################################

cor.mtx <- matrix(as.numeric(NA), ncol=nrow(rna.sce), nrow=length(TFs))
pvalue.mtx <- matrix(as.numeric(NA), ncol=nrow(rna.sce), nrow=length(TFs))
colnames(cor.mtx) <- rownames(rna.sce); rownames(cor.mtx) <- TFs
dimnames(pvalue.mtx) <- dimnames(cor.mtx)

opts$min.cor.TF2peak <- 0.00

for (i in TFs) {
  # i <- "FOXA2"
  
  peaks <- names(which(abs(dropNA2matrix(assay(tf2peak_cor.se[,i], "cor")))[,1]>=opts$min.cor.TF2peak))
  genes <- unique(peak2gene_nearest.dt[peaks,gene])
  genes <- genes[genes%in%rownames(rna.sce)]
  
  # calculate correlations
  corr_output <- psych::corr.test(t(logcounts(rna.sce.tf[i,])), t(logcounts(rna.sce[genes,])), ci=FALSE)

  # Fill matrices
  cor.mtx[i,genes] <- round(corr_output$r[1,],3)
  pvalue.mtx[i,genes] <- round(corr_output$p[1,],5)
}

# Save SummarizedExperiment object
to.save <- SummarizedExperiment(
  assays = SimpleList("cor" = dropNA(cor.mtx), "pvalue" = dropNA(pvalue.mtx))
)
saveRDS(to.save, sprintf("%s/cor_TFexpr_vs_Geneexpr_SummarizedExperiment_mincor%s.rds",io$outdir,opts$min.cor.TF2peak))

