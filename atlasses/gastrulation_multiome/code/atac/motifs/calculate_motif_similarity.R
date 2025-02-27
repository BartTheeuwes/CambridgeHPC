#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
  source("/Users/ricard/gastrulation_multiome_10x/atac/motifs/utils.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
  source("/homes/ricard/gastrulation_multiome_10x/atac/motifs/utils.R")
} else {
  stop("Computer not recognised")
}

# I/O
io$motifs <- paste0(io$basedir,"/results/atac/archR/motif_seqlogo/PWMatrixList.rds")
io$outdir <- paste0(io$basedir,"/results/atac/archR/motif_seqlogo")


##############
## Load PWM ##
##############

motifs <- readRDS(io$motifs)

##################################
## Calculate motif similarities ##
##################################

motif_similarity.mtx <- matrix(nrow=length(motifs), ncol=length(motifs))
rownames(motif_similarity.mtx) <- colnames(motif_similarity.mtx) <- names(motifs)
diag(motif_similarity.mtx) <- 1

      
for (i in 1:length(names(motifs))) {
  motif1 <- names(motifs)[i]
  print(motif1)
  for (j in i:length(names(motifs))) {
    if (i!=j) {
      motif2 <- names(motifs)[j]
      m1 <- 0.25*exp(as.matrix(motifs[[motif1]]))
      m2 <- 0.25*exp(as.matrix(motifs[[motif2]]))
      tryCatch(tmp <- motifSimilarity(m1, m2, trim = 0.4, self.sim = FALSE), error = function(e) { NA })
      motif_similarity.mtx[i,j] <- motif_similarity.mtx[j,i] <- tmp
    }
  }
}

# save
saveRDS(motif_similarity.mtx, paste0(io$outdir,"/motif_similarity.rds"))

# Load pre-computed
motif_similarity.mtx <- readRDS(paste0(io$outdir,"/motif_similarity.rds"))

##############
## Validate ##
##############

hist(motif_similarity.mtx)

#########################################################
## Load RNA expression and chromVAR scores, pseudobulk ##
#########################################################

opts$motif_annotation <- "Motif_cisbp"
io$archR.pseudobulk.deviations.se <- sprintf("%s/pseudobulk/pseudobulk_DeviationMatrix_%s_summarized_experiment.rds",io$archR.directory,opts$motif_annotation)

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/rna_atac/load_rna_atac_pseudobulk.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/rna_atac/load_rna_atac_pseudobulk.R")
} else {
  stop("Computer not recognised")
}

##########
## Plot ##
##########

# Select variable TFs
tmp <- apply(logcounts(rna.sce.tf),1,var)
tfs.to.plot <- tmp[tmp>0.5] %>% names
stopifnot(tfs.to.plot %in% rownames(motif_similarity.mtx))

to.plot <- motif_similarity.mtx[tfs.to.plot,tfs.to.plot]
to.plot[to.plot<0] <- 0

pheatmap::pheatmap(
  mat = to.plot, 
  cluster_rows  = T, 
  cluster_cols = T, 
  show_rownames = F, 
  show_colnames = F, 
  treeheight_row = 0,
  treeheight_col = 0,
  fontsize = 3,
  filename = paste0(io$outdir,"/motif_similarity_heatmap_all.pdf"),
  width = NA, height = NA
)


# HOX genes
tfs.to.plot <- grep("^HOX",rownames(motif_similarity.mtx),value=T)
to.plot <- motif_similarity.mtx[tfs.to.plot,tfs.to.plot]
to.plot[to.plot<0] <- 0

pheatmap::pheatmap(
  mat = to.plot, 
  cluster_rows  = T, 
  cluster_cols = T, 
  show_rownames = T, 
  show_colnames = T, 
  treeheight_row = 0,
  treeheight_col = 0,
  fontsize = 6,
  filename = paste0(io$outdir,"/motif_similarity_heatmap_HOX.pdf"),
  width = NA, height = NA
)


# FOX TFs
tfs.to.plot <- grep("^FOX",rownames(motif_similarity.mtx),value=T)
to.plot <- motif_similarity.mtx[tfs.to.plot,tfs.to.plot]
to.plot[to.plot<0] <- 0

pheatmap::pheatmap(
  mat = to.plot, 
  cluster_rows  = T, 
  cluster_cols = T, 
  show_rownames = T, 
  show_colnames = T, 
  treeheight_row = 0,
  treeheight_col = 0,
  fontsize = 6,
  filename = paste0(io$outdir,"/motif_similarity_heatmap_FOX.pdf"),
  width = NA, height = NA
)

# GATA TFs
tfs.to.plot <- grep("^GATA",rownames(motif_similarity.mtx),value=T)
to.plot <- motif_similarity.mtx[tfs.to.plot,tfs.to.plot]
to.plot[to.plot<0] <- 0

pheatmap::pheatmap(
  mat = to.plot, 
  cluster_rows  = T, 
  cluster_cols = T, 
  show_rownames = T, 
  show_colnames = T, 
  treeheight_row = 0,
  treeheight_col = 0,
  fontsize = 8,
  filename = paste0(io$outdir,"/motif_similarity_heatmap_GATA.pdf"),
  width = NA, height = NA
)
