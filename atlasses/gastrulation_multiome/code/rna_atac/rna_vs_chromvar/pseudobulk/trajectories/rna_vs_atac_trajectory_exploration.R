library(pheatmap)
library(igraph)

#####################
## Define settings ##
#####################

# load default setings
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
  # "Epiblast",
  # "Primitive_Streak",
  # "Caudal_epiblast",
  # "PGC",
  # "Anterior_Primitive_Streak",
  # "Notochord",
  # "Def._endoderm",
  # "Gut",
  # "Nascent_mesoderm",
  # "Mixed_mesoderm",
  # "Intermediate_mesoderm",
  # "Caudal_Mesoderm",
  # "Paraxial_mesoderm",
  # "Somitic_mesoderm",
  # "Pharyngeal_mesoderm",
  # "Cardiomyocytes",
  # "Allantois",
  # "ExE_mesoderm",
  # "Mesenchyme",
  "Haematoendothelial_progenitors",
  # "Endothelium",
  "Blood_progenitors_1",
  "Blood_progenitors_2",
  "Erythroid1",
  "Erythroid2",
  "Erythroid3"
  # "NMP",
  # "Rostral_neurectoderm",
  # "Caudal_neurectoderm",
  # "Neural_crest",
  # "Forebrain_Midbrain_Hindbrain",
  # "Spinal_cord",
  # "Surface_ectoderm",
  # "Visceral_endoderm",
  # "ExE_endoderm",
  # "ExE_ectoderm",
  # "Parietal_endoderm"
)

opts$motif_annotation <- "Motif_cisbp" # Motif_JASPAR2020_human

# I/O
# io$archR.pseudobulk.deviations.se <- sprintf("%s/pseudobulk/pseudobulk_DeviationMatrix_%s_summarized_experiment.rds",io$archR.directory,opts$motif_annotation)
io$archR.pseudobulk.deviations.se <- sprintf("%s/results/atac/archR/chromvar/pseudobulk/chromVAR_deviations_summarized_experiment_%s_pseudobulk_correlated_peaks.rds",io$basedir,opts$motif_annotation)
io$outdir <- paste0(io$basedir,"/results/rna_atac/rna_vs_chromvar/pseudobulk/trajectories/test")

#####################
## Load pseudotime ##
#####################

blood_trajectory.dt <- fread("/Users/ricard/data/gastrulation_multiome_10x/results/rna/trajectories/blood_trajectory/blood_trajectory.txt.gz")

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(io$metadata) %>%
  .[pass_rnaQC==TRUE & doublet_call==FALSE] %>%
  .[celltype.mapped%in%opts$celltypes] %>%
  .[,celltype.mapped:=factor(celltype.mapped,levels=opts$celltypes)]

sample_metadata <- sample_metadata[cell%in%blood_trajectory.dt$cell]

########################################
## Load single-cell RNA and ATAC data ##
########################################

sce <- load_SingleCellExperiment(io$rna.sce, normalise = TRUE, cells = sample_metadata$cell)
dim(sce)

# Select variable genes
decomp <- modelGeneVar(sce)
decomp <- decomp[decomp$mean > 0.1,]
hvgs <- rownames(decomp)[decomp$p.value <= 0.10]

hvgs <- hvgs[grep("^mt-",hvgs,invert=T)]
hvgs <- hvgs[grep("^Rps-",hvgs,invert=T)]
hvgs <- hvgs[grep("Rik$",hvgs,invert=T)]

# Subset SingleCellExperiment
sce_filt <- sce[hvgs,]
dim(sce_filt)

# Denoise
opts$knn <- 25
# io$pca.rna <- paste0(io$basedir,"/results/rna/dimensionality_reduction/all_cells/E7.5_rep1-E7.5_rep2-E8.0_rep1-E8.0_rep2-E8.5_rep1-E8.5_rep2_pca_features2500_pcs30_batchcorrectionbysample.txt.gz")
# pca.rna <- fread(io$pca.rna) %>% matrix.please %>% .[sample_metadata$cell,]
pca.rna <- blood_trajectory.dt[,c("cell","PC1")] %>% matrix.please
rna.mtx <- smoother_aggregate_nearest_nb(mat=as.matrix(logcounts(sce_filt)), D=pdist(pca.rna), k=opts$knn)
colnames(rna.mtx) <- colnames(sce_filt)

#######################################
## Load pseudobulk RNA and ATAC data ##
#######################################

# if (grepl("ricard",Sys.info()['nodename'])) {
#   source("/Users/ricard/gastrulation_multiome_10x/rna_atac/rna_vs_acc/pseudobulk/load_rna_atac_pseudobulk.R")
# } else if (grepl("ebi",Sys.info()['nodename'])) {
#   source("/homes/ricard/gastrulation_multiome_10x/rna_atac/rna_vs_acc/pseudobulk/load_rna_atac_pseudobulk.R")
# } else {
#   stop("Computer not recognised")
# }

sce.pseudobulk <- readRDS(io$rna.pseudobulk.sce)[,opts$celltypes]

###############################
## Load motifmatcher results ##
###############################

motifmatcher.se <- readRDS(sprintf("%s/Annotations/Motif_cisbp-Matches-In-Peaks.rds",io$archR.directory))
colnames(motifmatcher.se) <- colnames(motifmatcher.se) %>% toupper %>% stringr::str_split(.,"_") %>% map_chr(1)
motifmatcher.se <- motifmatcher.se[,!duplicated(colnames(motifmatcher.se))]
tmp <- rowRanges(motifmatcher.se)
rownames(motifmatcher.se) <- sprintf("%s:%s-%s",seqnames(tmp), start(tmp), end(tmp))
# motifmatcher.se <- motifmatcher.se[unique(cor_dt$peak),]

######################################
## Load TF2peak correlation results ##
######################################

io$file <- paste0(io$basedir,"/results/rna_atac/rna_vs_acc/pseudobulk/TFexpr_vs_peakAcc/cor_TFexpr_vs_peakAcc_SummarizedExperiment.rds")
tf2peak_cor.se <- readRDS(io$file)

#############################
## Load peak2gene linkages ##
#############################

peak2gene.dt <- fread(io$archR.peak2gene.all) %>% 
  .[,peak:=sprintf("chr%s:%s-%s",chr,peak.start,peak.end)]

########################################################
## Correlations between pseudotime and RNA expression ##
########################################################

cor.gene.trajectory <- psych::corr.test(t(rna.mtx), pca.rna)
hist(cor.gene.trajectory$p)

# rna.diff <- fread("/Users/ricard/data/gastrulation_multiome_10x/results/rna/differential/Erythroid3_vs_Haematoendothelial_progenitors.txt.gz") %>%
#   .[,gene:=toupper(gene)] %>%
#   .[abs(logFC)>1.6] %>%
#   .[gene%in%colnames(tf2peak_cor.se)]

###############################################
## Heatmap of peaks vs celltypes for each TF ##
###############################################

genes.to.plot <- rna.diff$gene %>% head(n=3)

# i <- "GATA1"
for (i in genes.to.plot) {
  
  peaks.to.plot <- assay(tf2peak_cor.se[,i])[,1][assay(tf2peak_cor.se[,i])[,1]!=0]
  peaks.to.plot <- sort(peaks.to.plot)
  
  to.plot.rna <- logcounts(rna.sce[i,])
  to.plot.acc <- assay(atac.peakMatrix.se[names(peaks.to.plot),])
  
  # Filter peaks
  to.plot.acc <- to.plot.acc[apply(to.plot.acc,1,var)>0.001,]
  peaks.to.plot <- peaks.to.plot[names(peaks.to.plot)%in%rownames(to.plot.acc)]
  
  # Prepare side information
  annotation_col.df <- data.frame(celltype=opts$celltypes); rownames(annotation_col.df) <- opts$celltypes
  mycolors <- list("celltype"=opts$celltype.colors[opts$celltypes])
  
  annotation_row.df <- data.frame(cor=peaks.to.plot)
  
  
  # Plot heamtap
  p.rna <- pheatmap(
    mat = to.plot.rna, 
    cluster_cols = F, cluster_rows = F,
    color = colorRampPalette(c("gray80", "purple"))(100),
    show_colnames = FALSE, show_rownames = FALSE,
    annotation_col = annotation_col.df,
    annotation_colors = mycolors,
    legend = F, annotation_legend = F,
    silent = TRUE
  )
  
  p.atac_clustered <- pheatmap(
    mat = to.plot.acc, 
    cluster_cols = F, cluster_rows = T,
    annotation_row = annotation_row.df,
    show_colnames = FALSE, show_rownames = FALSE,
    legend = F, annotation_legend = F,
    annotation_colors =  list("cor"=colorRampPalette(c("red", "blue"))(10)),
    treeheight_row = 10,
    silent = TRUE
  )
  
  # p.atac <- pheatmap(
  #   mat = to.plot.acc, 
  #   cluster_cols = F, cluster_rows = F,
  #   annotation_row = annotation_row.df,
  #   show_colnames = FALSE, show_rownames = FALSE,
  #   legend = F, annotation_legend = F,
  #   annotation_colors =  list("cor"=colorRampPalette(c("red", "blue"))(10)),
  #   treeheight_row = 0,
  #   silent = TRUE
  # )
  
  p <- cowplot::plot_grid(p.rna$gtable, p.atac_clustered$gtable, nrow = 2, rel_heights = c(1/10,9/10))
  pdf(sprintf("%s/%s_test_clustered.pdf",io$outdir,i), width=5, height=10)
  print(p)
  dev.off()
  
  # p <- cowplot::plot_grid(p.rna$gtable, p.atac$gtable, nrow = 2, rel_heights = c(1/10,9/10))
  # pdf(sprintf("%s/%s_test.pdf",io$outdir,i), width=5, height=10)
  # print(p)
  # dev.off()
}

########################################################
## Calculate correlation between TFs and target genes ##
########################################################

# TO-DO: INCLUDE CORRELATIONS BETWEEN TFS

TFs <- rownames(atac.chromvar.se)[rownames(atac.chromvar.se)%in%toupper(rownames(sce))]

# Select TFs that are highly variable along the trajectory
# TFs <- TFs[apply(logcounts(sce)[stringr::str_to_title(TFs),],1,var)>2]
# target_genes <- rownames(sce)[apply(logcounts(sce),1,var)>2]; target_genes <- target_genes[!target_genes%in%TFs]
# TFs <- TFs[stringr::str_to_title(TFs)%in%hvgs]
# target_genes <- hvgs[!toupper(hvgs)%in%TFs]
corr_genes <- names(which(cor.gene.trajectory$p[,1]<0.001))
TFs <- toupper(corr_genes)[toupper(corr_genes) %in% rownames(atac.chromvar.se)]
target_genes <- corr_genes[!toupper(corr_genes)%in%TFs]

# sce.TFs <- sce[stringr::str_to_title(TFs),]; rownames(sce.TFs) <- toupper(rownames(sce.TFs))
# sce.target_genes <- sce[target_genes,]
# atac.chromvar.se.filt <- atac.chromvar.se[TFs,]
rna.mtx.TFs <- rna.mtx[stringr::str_to_title(TFs),]; rownames(rna.mtx.TFs) <- toupper(rownames(rna.mtx.TFs))
rna.mtx.target_genes <- rna.mtx[rownames(rna.mtx) %in% target_genes,]

sce.pseudobulk.TFs <- sce.pseudobulk[stringr::str_to_title(TFs),]; rownames(sce.pseudobulk.TFs) <- toupper(rownames(sce.pseudobulk.TFs)) 
sce.pseudobulk.target_genes <- sce.pseudobulk[target_genes,]

peak2gene.dt.filt <- peak2gene.dt[gene%in%target_genes & dist<1e4]
peaks <- unique(peak2gene.dt.filt$peak)

motifmatcher.se.filt <- motifmatcher.se[,TFs]
motifmatcher.se.filt <- motifmatcher.se.filt[peaks,]
# motifmatcher.se.filt <- motifmatcher.se.filt[toupper(rowRanges(motifmatcher.se.filt)$nearestGene)%in%target_genes,]

corr.mtx <- psych::corr.test(t(rna.mtx.TFs),t(rna.mtx.target_genes))
# corr.mtx <- psych::corr.test(t(logcounts(sce.TFs)),t(logcounts(sce.target_genes)))

rownames(corr.mtx$r) <- toupper(rownames(corr.mtx$r))
rownames(corr.mtx$p) <- toupper(rownames(corr.mtx$p))

####################################################################################
## Filter peaks by correlation between RNA expression and chromatin accessibility ##
####################################################################################

# TO-DO: RE-DO CORRELATION USING ONLY THE TRAJECTORY

io$file <- paste0(io$basedir,"/results/rna_atac/rna_vs_acc/pseudobulk/TFexpr_vs_peakAcc/cor_TFexpr_vs_peakAcc.txt.gz")
tf2peak_cor.dt <- fread(io$file) %>%
  .[!is.na(cor) & peak%in%peaks] %>%
  .[,cor_sign:=c("-","+")[(cor>0)+1]]

io$file <- paste0(io$basedir,"/results/rna_atac/rna_vs_acc/pseudobulk/TFexpr_vs_peakAcc/cor_TFexpr_vs_peakAcc_SummarizedExperiment.rds")
tf2peak_cor.se <- readRDS(io$file)[peaks,]

# tf2peak_cor.dt[TF=="GATA1"]

corr.mtx.filt <- corr.mtx
for (i in TFs) {
  
  # Filter peaks by RNA expr & ATAC cor
  target_peaks_i <- which(dropNA2matrix(assay(tf2peak_cor.se[,i],"pvalue"))[,1]<=0.15 & abs(dropNA2matrix(assay(tf2peak_cor.se[,i],"cor"))[,1])>0.15) %>% names
  
  # Filter peaks by motif presence (no RNAexpr & ATACcor)
  # target_peaks_i <- names(which(assay(motifmatcher.se.filt[,i])[,1]))
  
  target_genes_i <- peak2gene.dt.filt[peak%in%target_peaks_i,gene] %>% unique
  corr.mtx.filt$p[i,!colnames(corr.mtx.filt$p)%in%target_genes_i] <- NA
  corr.mtx.filt$r[i,!colnames(corr.mtx.filt$p)%in%target_genes_i] <- NA
}

mean(is.na(corr.mtx.filt$p))
mean(is.na(corr.mtx.filt$r))

##########################################
## Create network of TF to target genes ##
##########################################

# tmp <- corr.mtx.filt$p
# tmp[tmp>0.50] <- NA
# tmp[tmp<1e-100] <- 1e-100
# tmp <- minmax.normalisation(-log(tmp))
# tmp <- tmp[rowSums(!is.na(tmp))>=5,]
# tmp <- tmp[,colSums(!is.na(tmp))>=1]

tmp <- corr.mtx.filt$r
tmp[abs(tmp)<0.25] <- NA
tmp <- tmp[rowSums(!is.na(tmp))>=5,]
tmp <- tmp[,colSums(!is.na(tmp))>=1]

# Prepare data
node_list.dt <- data.table(node_id=1:nrow(tmp), node_name=rownames(tmp))
target_list.dt <- data.table(target_id=1:ncol(tmp), target_name=colnames(tmp))

edge_list.dt <- as.data.table(tmp,keep.rownames = T) %>% 
  setnames("rn","from") %>%
  melt(id.vars=c("from"), variable.name="to", value.name="weight") %>%
  .[!is.na(weight)]


node_list_metadata.dt <- data.table(
  label = c(node_list.dt$node_name, target_list.dt$target_name),
  class = c(rep("TF",nrow(node_list.dt)), rep("gene",nrow(target_list.dt)))
)


net <- graph_from_data_frame(d = edge_list.dt, vertices = node_list_metadata.dt)

##################
## Plot network ##
##################

# Modify nodes
V(net)$group <- factor(V(net)$class, levels=c("TF","gene"))
V(net)$size <- c(2.75,6)[factor(V(net)$class)]
V(net)$label.color <- "black"
# V(net)$color <- c("tomato", "gold")[factor(V(net)$class)]
V(net)$color <- "gray60"
V(net)$label.cex <- c(0.15,0.5)[factor(V(net)$class)]

# Modify edges
E(net)$arrow.size <- 0.025
# E(net)$edge.color <- "gray80"
# E(net)$edge.width <- 0.25
# E(net)$width <- E(net)$weight/2

E(net)[which(E(net)$weight<0)]$color <- "darkblue"  # Colour negative correlation edges as blue
E(net)[which(E(net)$weight>0)]$color <- "darkred"   # Colour positive correlation edges as red
E(net)$weight <- abs(E(net)$weight)   # Convert edge weights to absolute values
# net <- delete_edges(net, E(net)[which(E(net)$weight<0.8)]) # Remove edges below absolute Pearson correlation 0.8
# net <- delete.vertices(net, degree(net)==0)  # Remove any vertices remaining that have no edges

# set.seed(42)
# plot(net, layout = layout_with_fr)

pdf(sprintf("%s/blood_network_test.pdf",io$outdir), width=10, height=9, bg="white")
# png(sprintf("%s/blood_network_test.png",io$outdir), width=900, height=800)
set.seed(42)
plot(net, layout = layout_with_fr)
dev.off()

# pal <- grDevices::colorRamp(c("gray60", "purple"))( (1:100)/100 )
# pal <- cbind(pal, seq(100, 255, length.out = 100))

for (i in opts$celltypes) {
  # rna.expr.i <- c(logcounts(sce.pseudobulk)[TFs,i], logcounts(sce.target_genes)[tolower(target_genes),i])[names(V(net))]
  rna.expr.i <- c(logcounts(sce.pseudobulk.TFs)[TFs,i], logcounts(sce.pseudobulk.target_genes)[target_genes,i])[names(V(net))] %>% minmax.normalisation
  V(net)$color <- colourvalues::colour_values(rna.expr.i, palette = "viridis")
  
  pdf(sprintf("%s/%s_blood_network_test.pdf",io$outdir,i), width=10, height=9, bg="white")
  # png(sprintf("%s/%s_blood_network_test.png",io$outdir,i), width=900, height=800)
  set.seed(42)
  plot(net, layout = layout_with_fr, vertex.label=NA, main=i)
  dev.off()
}


################
## Create GIF ##
################

library(magick)

## list file names and read in
imgs.files <- c(
  sprintf("%s/blood_network_test.pdf",io$outdir),
  sprintf("%s/%s_blood_network_test.pdf",io$outdir,opts$celltypes)
)
img_list <- lapply(imgs.files, image_read)

# join the images together
img_joined <- image_join(img_list)

# animate
img_animated <- image_animate(img_joined, fps = 1)

# save
image_write(image = img_animated, path = sprintf("%s/test.gif",io$outdir))


