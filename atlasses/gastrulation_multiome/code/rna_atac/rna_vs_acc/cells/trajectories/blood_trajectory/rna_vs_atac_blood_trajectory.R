library(RColorBrewer)
library(ggnewscale)
library(pheatmap)

#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
  source("/Users/ricard/gastrulation_multiome_10x/utils.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
  source("/homes/ricard/gastrulation_multiome_10x/utils.R")
} else {
  stop("Computer not recognised")
}

########################
## Load ArchR Project ##
########################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/atac/archR/load_archR_project.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/atac/archR/load_archR_project.R")
} else {
  stop("Computer not recognised")
}

#####################
## Define settings ##
#####################

# I/O
io$outdir <- paste0(io$basedir,"/results/rna_atac/rna_vs_acc/trajectories/blood_trajectory")
# io$pca.rna <- paste0(io$basedir,"/results/rna/dimensionality_reduction/all_cells/E7.5_rep1-E7.5_rep2-E8.0_rep1-E8.0_rep2-E8.5_rep1-E8.5_rep2_pca_features2500_pcs30_batchcorrectionbysample.txt.gz")
# io$pca.atac <- paste0(io$basedir,"/results/atac/archR/dimensionality_reduction/PeakMatrix/all_cells/E7.5_rep1-E7.5_rep2-E8.0_rep1-E8.0_rep2-E8.5_rep1-E8.5_rep2_lsi_features50000_ndims30.txt.gz")
io$trajectory <- paste0(io$basedir,"/results/rna/trajectories/blood_trajectory/blood_trajectory.txt.gz")

# Options

opts$celltypes = c(
  # "Mixed_mesoderm",
  "Haematoendothelial_progenitors",
  "Blood_progenitors_1",
  "Blood_progenitors_2",
  "Erythroid1",
  "Erythroid2",
  "Erythroid3"
)

opts$min.expr <- 0.1
opts$motif_annotation <- "Motif_cisbp"

#####################
## Load trajectory ##
#####################

trajectory.dt <- fread(io$trajectory)

###################
## Load metadata ##
###################

sample_metadata <- fread(io$metadata) %>%
  .[pass_atacQC==TRUE & pass_rnaQC==TRUE] %>%
  .[celltype.predicted%in%opts$celltypes & cell%in%trajectory.dt$cell] %>%
  setnames("celltype.predicted","celltype")

# Filter cells
trajectory.dt <- trajectory.dt[cell%in%sample_metadata$cell]
sample_metadata <- sample_metadata[cell%in%trajectory.dt$cell] %>% setkey(cell) %>% .[trajectory.dt$cell]
stopifnot(!is.na(sample_metadata$sample))
stopifnot(sample_metadata$cell==trajectory.dt$cell)

##################
## Subset ArchR ##
##################

ArchRProject <- ArchRProject[sample_metadata$cell,]

####################################################
## Load genes that correlate along the trajectory ##
####################################################

# io$cor_gene_trajectory <- sprintf("%s/results/rna/trajectories/blood_trajectory/correlation_gene_vs_pseudotime_blood.txt.gz",io$basedir)
# cor_expr_vs_pseudotime.dt <- fread(io$cor_gene_trajectory)

#####################################################
## Load RNA expression and ATAC peak accessibility ##
#####################################################

# io$rna_atac.file <- sprintf("%s/chromvar_rna.txt.gz",args$outdir)
# if (file.exists(io$rna_atac.file)) {
#   
#   sprintf("Found precomputed values, loading %s...",io$rna_chromvar.file)
#   chromvar_rna_dt <- fread(io$rna_chromvar.file)
#   
# } else {
#   
if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/rna_atac/rna_vs_acc/load_rna_atac_single_cells.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/rna_atac/rna_vs_acc/load_rna_atac_single_cells.R")
} else {
  stop("Computer not recognised")
}
# }

# Filter genes
rna.sce <- rna.sce[!grepl("Rik|Gm|mt-|Rps|Rpl",rownames(rna.sce)),]
rna.sce <- rna.sce[apply(logcounts(rna.sce),1,var)>0.25,]
# rna.mtx <- rna.mtx[apply(rna.mtx,1,var)>0.5,]
# rna.mtx <- rna.mtx[!grepl("Rik|Gm|mt-|Rps|Rpl",rownames(rna.mtx)),]

# Filter peaks
atac.peakMatrix.se <- atac.peakMatrix.se[Matrix::rowSums(assay(atac.peakMatrix.se))>50,]
# atac.peak.mtx <- atac.peak.mtx[rowSums(atac.peak.mtx)>25,]

##############################################
## Load association between peaks and genes ##
##############################################

peak2gene.dt <- fread(io$archR.peak2gene.all) %>%
  .[gene%in%rownames(rna.sce)] %>%
  .[,peak:=sub("_",":",peak)] %>%
  .[,peak:=sub("_","-",peak)] %>%
  .[,peak:=paste0("chr",peak)] %>%
  .[peak%in%rownames(atac.peakMatrix.se)]

opts$max.dist <- 5e4
peak2gene.dt <- peak2gene.dt[dist<opts$max.dist]

# Load DORCs
# io$dorcs <- paste0(io$basedir,"/results/rna_atac/DORCs/DORCs_correlation_estimates_E7.5_rep1-E7.5_rep2-E8.5_rep1-E8.5_rep2_25000.tsv.gz")
# DORCs.dt <- fread(io$dorcs)

# Subset genes and peaks with DORCs
# atac.peak.mtx <- atac.peak.mtx[rownames(atac.peak.mtx) %in% unique(DORCs.dt$peak)]
# rna.mtx <- rna.mtx[rownames(rna.mtx) %in% unique(DORCs.dt$gene)]

# io$dorcs <- paste0(io$basedir,"/results/rna_atac/DORCs_v2/DORCs_samplesE7.5_rep1-E7.5_rep2-E8.0_rep1-E8.0_rep2-E8.5_rep1-E8.5_rep2_distance10000_noExE.rds")
# DORCs <- readRDS(io$dorcs)

#############
## Denoise ##
#############
stop("WRONG")

opts$knn <- 50

trajectory.mtx <- trajectory.dt[,c("cell","PC1")] %>% matrix.please

# RNA
rna.mtx <- smoother_aggregate_nearest_nb(mat=as.matrix(logcounts(rna.sce)), D=pdist(trajectory.mtx), k=opts$knn)
colnames(rna.mtx) <- colnames(rna.sce)

# RNA TF
rna_tf.mtx <- smoother_aggregate_nearest_nb(mat=as.matrix(logcounts(rna.sce.tf)), D=pdist(trajectory.mtx), k=opts$knn)
colnames(rna_tf.mtx) <- colnames(rna.sce.tf)

# ATAC
atac.peak.mtx <- smoother_aggregate_nearest_nb(mat=as.matrix(assay(atac.peakMatrix.se)), D=pdist(trajectory.mtx), k=opts$knn)
colnames(atac.peak.mtx) <- colnames(atac.peakMatrix.se)

###############################################
## Correlate gene expression with pseudotime ##
###############################################

cor.output.rna <- cor(trajectory.mtx[,1],t(rna.mtx))[1,]

##################################################
## Correlate peak accessibility with pseudotime ##
##################################################

cor.output.atac <- cor(trajectory.mtx[,1],t(atac.peak.mtx))[1,]
cor.output.atac <- cor.output.atac[abs(cor.output.atac)>0.35]

#############
## Explore ##
#############

# genes.to.explore <- cor.output[abs(cor.output)>0.25] %>% abs %>% sort %>% tail(n=50) %>% names
genes.to.explore <- cor.output[abs(cor.output)>0.75] %>% abs %>% sort %>% names

# annotation_col.df <- trajectory.dt[,c("cell","PC1")]
# mycolors <- colorRampPalette(brewer.pal(9, "Reds"))(length(unique(annotation_col.df$PC1)))
# names(mycolors) <- unique(annotation_col.df$V1)

annotation_col.df <- sample_metadata[,c("cell","celltype")] %>% tibble::column_to_rownames("cell")
mycolors <- list(celltype = opts$celltype.colors[unique(sample_metadata$celltype)])

genes.to.explore <- i <- "Hemgn"
for (i in genes.to.explore) {
  # Select peaks
  peak_idx <- peak2gene.dt[gene==i,peak]
  peak_idx <- peak_idx[peak_idx%in%rownames(atac.peak.mtx)]
  peak_idx <- peak_idx[peak_idx%in%names(cor.output.atac)]
  
  
  if (length(peak_idx)>0) {
    
    #############
    ## Heatmap ##
    #############
    
    # RNA
    to.plot.rna <- rna.mtx[i,trajectory.dt$cell,drop=F]
    
    # ATAC
    # to.plot.acc <- atac.peakMatrix.se[peak_idx,] %>% assay %>% as.matrix %>% t %>% as.data.table(keep.rownames = T) %>%
    #   setnames("rn","cell") %>% 
    #   .[cell%in%colnames(to.plot.rna)] %>%
    #   matrix.please %>% t
    to.plot.acc <- atac.peak.mtx[peak_idx,,drop=F]
    
    # Sort cells according to the pseudotime values
    # to.plot.rna <- to.plot.rna[,trajectory.dt$cell,drop=F]
    # rownames(to.plot.rna) <- i
    to.plot.acc <- to.plot.acc[,trajectory.dt$cell]
    stopifnot(all(rownames(annotation_col.df)==colnames(to.plot.rna)))
    
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
    p.atac <- pheatmap(
      mat = to.plot.acc, 
      cluster_cols = F, cluster_rows = T,
      color = colorRampPalette(c("white", "black"))(100),
      show_colnames = FALSE, show_rownames = FALSE,
      legend = F,
      treeheight_row = 0,
      silent = TRUE
    )
    
    p.heatmap <- cowplot::plot_grid(p.rna$gtable, p.atac$gtable, nrow = 2, rel_heights = c(1/5,4/5))
    
    ##############
    ## Lineplot ##
    ##############
    
    rna_dt <- data.table(cell=colnames(rna.mtx), feature=i, value=rna.mtx[i,]) %>%
      .[,modality:="RNA"]
    
    # acc_dt <- atac.peakMatrix.se[peak_idx,] %>% assay %>% as.matrix %>% t %>% as.data.table(keep.rownames = T) %>%
    acc_dt <- atac.peak.mtx[peak_idx,,drop=F] %>% t %>% as.data.table(keep.rownames = T) %>%
      setnames("rn","cell") %>% 
      melt(id.vars=c("cell"), variable.name="feature", value.name="value") %>%
      .[,modality:="ATAC"]
    
    rna_acc.dt <- rbind(rna_dt, acc_dt) %>%
      merge(sample_metadata[,c("cell","celltype")], by="cell") %>%
      merge(trajectory.dt,by="cell")
    
    to.plot <- rna_acc.dt %>% 
      .[,value:=value/max(value),by=c("feature")]
    
    p.lineplot <- ggplot(to.plot, aes(x=PC1, y=value, group=feature)) +
      stat_smooth(aes(color=modality), fill=NA, method="loess", alpha=0.75, span=0.5) +
      scale_color_manual(values=c("ATAC"="black", "RNA"="purple")) +
      new_scale_color() +
      geom_rug(aes(color=celltype), sides="b") +
      scale_color_manual(values=opts$celltype.colors) +
      guides(color=F) +
      labs(x="Pseudotime", y=i) +
      theme_classic() +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        legend.position="top"
      )
    
    #############
    ## Combine ##
    #############
    
    p <- cowplot::plot_grid(p.heatmap, p.lineplot, ncol = 2, rel_widths = c(1/2,1/2))
    
    pdf(sprintf("%s/individual_genes/%s_rna_atac_vs_pseudotime_blood.pdf",io$outdir,i), width=9, height=4)
    print(p)
    dev.off()
    
  }
  
} 



################
## parse data ##
################

# chromvar_dt <- atac.deviation.mtx %>% 
chromvar_dt <- atac.deviation.mtx.smoothed %>% 
  .[motif2gene.dt$motif,] %>%
  as.matrix %>% t %>% as.data.table(keep.rownames = T) %>%
  setnames("rn","cell") %>% 
  melt(id.vars=c("cell"), variable.name="motif", value.name="chromvar_zscore")

# rna_dt <- sce %>%
#   .[motif2gene.dt$gene,] %>%
#   logcounts %>% as.matrix %>% 
rna_dt <- rna.mtx.smoothed[motif2gene.dt$gene,] %>%
  as.data.table(keep.rownames = T) %>%
  setnames("rn","gene") %>%
  melt(id.vars="gene", variable.name="cell", value.name="expr")



##########
## Plot ##
##########

genes.to.plot <- chromvar_rna_dt %>%
  .[,.(expr=mean(expr),chromvar_zscore=mean(chromvar_zscore)),by=c("gene","motif")] %>% 
  .[expr>opts$min.expr,gene]

# genes.to.plot <- unique(chromvar_rna_dt$gene) %>% head(n=3)

for (i in genes.to.plot) {
  
  to.plot <- chromvar_rna_dt[gene==i] %>% 
    # .[chromvar_zscore>7,chromvar_zscore:=7] %>%
    melt(id.vars=c("cell","V1","celltype"), measure.vars=c("chromvar_zscore","expr"), variable.name="modality")# %>%
    # .[,value_scaled:=value/max(value),by="modality"]
  # .[,value_scaled:=(value-min(value))/(max(value)-min(value)), by="modality"]
    
  p1 <- ggplot(to.plot, aes(x=V1, y=value)) +
    geom_point(aes(fill=celltype), size=1.25, shape=21, stroke=0.1) +
    stat_smooth(method="loess", color="black", alpha=0.75, span=0.5) +
    geom_rug(aes(color=celltype), sides="b") +
    facet_wrap(~modality, nrow=2, scales="free_y") +
    scale_color_manual(values=opts$celltype.colors) +
    scale_fill_manual(values=opts$celltype.colors) +
    guides(fill=F, color=F) +
    # scale_fill_manual(values=opts$celltype.colors) +
    # scale_fill_brewer(palette="Dark2") +
    labs(x="Pseudotime", y=i) +
    theme_classic() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.title = element_blank(),
      legend.position="top"
    )
  
  p2 <- ggboxplot(to.plot, x="celltype", y="value", fill="celltype", outlier.shape=NA) +
    scale_fill_manual(values=opts$celltype.colors) +
    labs(x="", y="") +
    # stat_compare_means(comparisons = list( c("Somitic_mesoderm", "NMP"), c("NMP", "Spinal_cord") ), label="p.signif", hide.ns=T) +
    facet_wrap(~modality, nrow=2, scales="free_y") +
    # coord_cartesian(ylim=c(0,5)) +
    theme_classic() +
    guides(x = guide_axis(angle = 90)) +
    theme(
      # axis.text.x = element_text(color="black", size=rel(0.75)),
      # legend.title = element_blank(),
      legend.position = "none",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  
  p <- cowplot::plot_grid(plotlist=list(p1,p2), nrow = 1, rel_widths = c(1/2,1/2))
  
  pdf(sprintf("%s/individual_genes/%s_rna_chromvar_vs_pseudotime_blood.pdf",io$outdir,i), width=8, height=6)
  print(p)
  dev.off()
}




##############################################
## Plot RNA expression of TFs vs pseudotime ##
##############################################

rna_tf.mtx <- smoother_aggregate_nearest_nb(mat=as.matrix(logcounts(rna.sce.tf)), D=pdist(trajectory.mtx), k=500)
colnames(rna_tf.mtx) <- colnames(rna.sce.tf)

cor.output.rna <- cor(trajectory.mtx[,1],t(rna_tf.mtx))[1,]
cor.output.rna <- cor.output.rna[abs(cor.output.rna)>0.75 & !is.na(cor.output.rna)]
blood.tfs <- intersect(names(which(apply(rna_tf.mtx.filt,1,mean)>0.25)),names(cor.output.rna))

rna_tf.mtx.filt <- rna_tf.mtx[blood.tfs,trajectory.dt$cell]

annotation_col.df <- sample_metadata[,c("cell","celltype")] %>% tibble::column_to_rownames("cell")
mycolors <- list(celltype = opts$celltype.colors[unique(sample_metadata$celltype)])
stopifnot(all(rownames(annotation_col.df)==colnames(rna_tf.mtx.filt)))

pheatmap(
  mat = rna_tf.mtx.filt, 
  cluster_cols = F, cluster_rows = T,
  # color = colorRampPalette(c("gray80", "purple"))(100),
  show_colnames = FALSE, show_rownames = TRUE,
  annotation_col = annotation_col.df,
  annotation_colors = mycolors,
  scale = "row",
  legend = F, annotation_legend = F
)

#################################
## Find genes that correlate with blood ##
#################################

# RNA TF
rna.mtx <- smoother_aggregate_nearest_nb(mat=as.matrix(logcounts(rna.sce)), D=pdist(trajectory.mtx), k=250)
colnames(rna.mtx) <- colnames(rna.sce)

cor.output.rna <- cor(trajectory.mtx[,1],t(rna.mtx))[1,]
cor.output.rna <- cor.output.rna[abs(cor.output.rna)>0.5 & !is.na(cor.output.rna)]
blood.genes <- intersect(names(which(apply(rna.mtx.filt,1,mean)>0.25)),names(cor.output.rna))
rna.mtx.filt <- rna.mtx[blood.genes,trajectory.dt$cell]

pheatmap(
  mat = rna.mtx.filt, 
  cluster_cols = F, cluster_rows = T,
  show_colnames = FALSE, show_rownames = TRUE,
  annotation_col = annotation_col.df,
  annotation_colors = mycolors,
  scale = "row",
  legend = F, annotation_legend = F
)


###############################
## Load motifmatcher results ##
###############################

source("/Users/ricard/gastrulation_multiome_10x/load_motifmatchR.R")

###################################################
## Load TF2peak correlation results (pseudobulk) ##
###################################################

io$tf2peak_cor.dt <- paste0(io$basedir,"/results/rna_atac/rna_vs_acc/pseudobulk/TFexpr_vs_peakAcc/cor_TFexpr_vs_peakAcc.txt.gz")
tf2peak_cor.dt <- fread(io$tf2peak_cor.dt) %>%
  .[!is.na(cor)] %>%
  .[,cor_sign:=c("-","+")[(cor>0)+1]]

io$tf2peak_cor.se <- paste0(io$basedir,"/results/rna_atac/rna_vs_acc/pseudobulk/TFexpr_vs_peakAcc/cor_TFexpr_vs_peakAcc_SummarizedExperiment.rds")
tf2peak_cor.se <- readRDS(io$tf2peak_cor.se)

#############################
## Load peak2gene linkages ##
#############################

peak2gene.dt <- fread(io$archR.peak2gene.all) %>% 
  .[,peak:=sprintf("chr%s:%s-%s",chr,peak.start,peak.end)]

#########################################################
## Correlate RNA expression of TFs versus target genes ##
#########################################################

cor_TFexpr_vs_GeneEXPR.se <- readRDS("/Users/ricard/data/gastrulation_multiome_10x/results/rna_atac/rna_vs_acc/pseudobulk/TFexpr_vs_Geneexpr/cor_TFexpr_vs_Geneexpr_SummarizedExperiment_mincor0.25.rds")
cor_TFexpr_vs_GeneEXPR.se <- cor_TFexpr_vs_GeneEXPR.se[rownames(cor_TFexpr_vs_GeneEXPR.se)%in%blood.tfs,colnames(cor_TFexpr_vs_GeneEXPR.se) %in% blood.genes]
# cor.mtx <- dropNA2matrix(assay(to.save,"cor"))
# pval.mtx <- dropNA2matrix(assay(to.save,"pvalue"))

# Compute correlations
# corr_TFexpr_vs_Geneexpr.mtx <- psych::corr.test(t(rna_tf.mtx.filt,t(rna.mtx.filt), ci = F))

##############################
## Load pseudobulk RNA data ##
##############################

sce.pseudobulk <- readRDS(io$rna.pseudobulk.sce)[,opts$celltypes]
sce.pseudobulk.TF <- sce.pseudobulk[toupper(rownames(sce.pseudobulk))%in%blood.tfs]
rownames(sce.pseudobulk.TF) <- toupper(rownames(sce.pseudobulk.TF))

#####################
## Create network ##
#####################

suppressMessages(library(GGally))
suppressMessages(library(igraph))
suppressMessages(library(network))
suppressMessages(library(sna))
suppressMessages(library(intergraph))

tmp <- dropNA2matrix(assay(cor_TFexpr_vs_GeneEXPR.se,"cor"))
tmp[abs(tmp)<0.50] <- NA
tmp <- tmp[,colSums(!is.na(tmp))>=1]
tmp <- tmp[rowSums(!is.na(tmp))>=10,]

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

# Define groups
V(net)$group <- factor(V(net)$class, levels=c("TF","gene"))


##################
## Plot network ##
##################

set.seed(42)
layout <- layout.fruchterman.reingold(net)
# layout <- layout.sphere(net)

celltypes.to.plot <- opts$celltypesi


pal <- grDevices::colorRamp(c("gray60", "purple"))( (1:100)/100 )
pal <- cbind(pal, seq(100, 255, length.out = 100))

rna.mtx.TFs.scaled <- apply(logcounts(sce.pseudobulk.TF[rownames(tmp),]),1,minmax.normalisation) %>% t
rna.mtx.target_genes.scaled <- apply(logcounts(sce.pseudobulk[colnames(tmp),]),1,minmax.normalisation) %>% t

for (i in celltypes.to.plot) {
  
  # Define node colors
  foo <- rna.mtx.TFs.scaled[,i]; names(foo) <- rownames(rna.mtx.TFs.scaled)
  rna.expr.i <- c(foo, rna.mtx.target_genes.scaled[,i])[names(V(net))]# %>% minmax.normalisation
  
  # Discretise a bit
  rna.expr.i[rna.expr.i<0.33] <- 0
  rna.expr.i[rna.expr.i>0.66] <- 1
  
  # V(net)$color <- colourvalues::colour_values(rna.expr.i, palette = "viridis")
  V(net)$color <- colourvalues::colour_values(rna.expr.i, palette = pal)
  
  # Define edge colors
  E(net)[which(E(net)$weight<0)]$color <- "darkblue"  # Colour negative correlation edges as blue
  E(net)[which(E(net)$weight>0)]$color <- "darkred"   # Colour positive correlation edges as red
  
  # Define edge width
  edge_weight <- abs(E(net)$weight)   # Convert edge weights to absolute values
  for (j in seq_len(length(E(net)))) {
    edge_weight[[j]] <- edge_weight[[j]]*(rna.expr.i[tail_of(net, j)] * rna.expr.i[head_of(net, j)])
  }
  edge_weight[edge_weight==0] <- 0.01
  # edge_weight <- minmax.normalisation(edge_weight)
  
  alphas <- rna.expr.i
  # alphas[rna.expr.i<0.5] <- 0.25
  # alphas[rna.expr.i>0.75] <- 1
  
  p <- ggnet2(
    mode = layout, # "fruchtermanreingold",
    net = net,
    color = V(net)$color,
    edge.color = E(net)$color,
    edge.size = edge_weight,
    node.size = c(2.5,10)[factor(V(net)$class)],
    node.alpha = alphas,
    label = TRUE,
    label.color = "black",
    label.size = c(1.5,3)[factor(V(net)$class)],
    arrow.size = 0.1,
    legend.position = "none"
  ) 
  
  outfile <- sprintf("%s/%s_%s_stringent.png",args$outdir,paste(args$genes,collapse="-"),i)
  png(outfile, width = 1200, height = 550, bg = "white")
  print(p)
  dev.off()
}
