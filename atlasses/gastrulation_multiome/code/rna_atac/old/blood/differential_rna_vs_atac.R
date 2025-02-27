
#####################
## Define settings ##
#####################

# Load default settings
if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
} else {
  stop("Computer not recognised")
}

# I/O
io$rna.differential <- paste0(io$basedir,"/results/rna/differential")
io$atac.differential <- paste0(io$basedir,"/results/atac/archR/differential/PeakMatrix")
io$outdir <- paste0(io$basedir,"/results/rna_atac/blood")

# Options
opts$samples <- c(
  "E7.5_rep1",
  "E7.5_rep2",
  "E8.5_rep1",
  "E8.5_rep2"
)

opts$celltypes <- c(
  # "ExE_mesoderm",
  # "Mesenchyme",
  # "Allantois",
  "Haematoendothelial_progenitors",
  "Blood_progenitors_1",
  "Blood_progenitors_2",
  "Erythroid1",
  "Erythroid2",
  "Erythroid3",
  "Endothelium"
)

opts$atac.matrix <- "PeakMatrix"
opts$min.Log2FC.rna <- 1
opts$min.diff.atac <- 0.1
opts$min.FDR <- 0.1

##############################################
## Load differential RNA expression results ##
##############################################

rna_diff.dt <- opts$celltypes %>% map(function(i) { opts$celltypes %>% map(function(j) {
  file <- sprintf("%s/%s_vs_%s.txt.gz", io$rna.differential,i,j)
  if (file.exists(file)) {
    fread(file, select = c(1,2,4)) %>% 
      .[,c("celltypeA","celltypeB"):=list(as.factor(i),as.factor(j))] %>%
    return
  }
}) %>% rbindlist }) %>% rbindlist %>%
  .[,celltypeA:=factor(celltypeA,levels=opts$celltypes)] %>%
  .[,celltypeB:=factor(celltypeB,levels=opts$celltypes)] %>%
  .[,sig:=abs(logFC)>=opts$min.Log2FC.rna & padj_fdr<=opts$min.FDR] %>%
  .[is.na(sig),sig:=FALSE]

####################################
## Load differential ATAC results ##
####################################

atac_diff.dt <- opts$celltypes %>% map(function(i) { opts$celltypes %>% map(function(j) {
  file <- sprintf("%s/%s_%s_vs_%s.txt.gz", io$atac.differential,opts$atac.matrix,i,j)
  if (file.exists(file)) {
    fread(file, select = c(2,6,7)) %>% 
      .[,c("celltypeA","celltypeB"):=list(as.factor(i),as.factor(j))] %>%
      return
  }
}) %>% rbindlist }) %>% rbindlist %>%
  .[,celltypeA:=factor(celltypeA,levels=opts$celltypes)] %>%
  .[,celltypeB:=factor(celltypeB,levels=opts$celltypes)] %>%
  .[,sig:=abs(MeanDiff)>=opts$min.diff.atac & FDR<=opts$min.FDR]

#########
## XXX ##
#########

# Calculate number of differential changes for each stage transition
rna_diff.sum <- rna_diff.dt[,mean(sig), by = c("celltypeA","celltypeB")]
atac_diff.sum <- atac_diff.dt[,mean(sig), by = c("celltypeA","celltypeB")]

# split by sign
rna_diff.sum <- rna_diff.dt %>%
  .[,sign:=c("Up","Down")[as.numeric(logFC>0)+1]] %>%
  .[!is.na(sign),mean(sig), by = c("celltypeA","celltypeB","sign")] %>%
  dcast(celltypeA+celltypeB~sign,value.var="V1") %>%
  .[,rna.diff:=Up-Down]

atac_diff.sum <- atac_diff.dt %>%
  .[,sign:=c("Up","Down")[as.numeric(MeanDiff>0)+1]] %>%
  .[!is.na(sign),mean(sig), by = c("celltypeA","celltypeB","sign")] %>%
  dcast(celltypeA+celltypeB~sign,value.var="V1") %>%
  .[,atac.diff:=Up-Down]

# CAN'T BECAUSE OF THE ORDER cellypeA celltypeB
options(scipen=999)

diff <- merge(
  rna_diff.sum[,c("celltypeA","celltypeB","rna.diff")], 
  atac_diff.sum[,c("celltypeA","celltypeB","atac.diff")], 
  by=c("celltypeA","celltypeB")
)
#############
## heatmap ##
#############

to.plot <- rna_diff.sum %>%
# to.plot <- atac_diff.sum %>%
  dcast(celltypeA~celltypeB, value.var="V1", drop=FALSE) %>%
  matrix.please 


# Fill NAs
for (i in rownames(to.plot)) {
  for (j in colnames(to.plot)) {
    if (is.na(to.plot[i,j])) to.plot[i,j] <- to.plot[j,i] 
  }
}
diag(to.plot) <- 0


# Plot heatmap
# pdf(sprintf("%s/heatmap_differential.pdf",io$outdir), width=8, height=6)
pheatmap::pheatmap(to.plot, cluster_rows  = F, cluster_cols = F, fontsize = 8)
# dev.off()



## diff

to.plot.diff <- to.plot.rna - to.plot.atac
pheatmap::pheatmap(to.plot.diff, cluster_rows  = F, cluster_cols = F, fontsize = 8)
