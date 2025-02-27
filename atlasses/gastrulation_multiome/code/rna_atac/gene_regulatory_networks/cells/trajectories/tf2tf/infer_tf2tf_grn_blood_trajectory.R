source(here::here("settings.R"))
source(here::here("utils.R"))

suppressMessages(library(GGally))
suppressMessages(library(igraph))
suppressMessages(library(network))
suppressMessages(library(sna))
suppressMessages(library(intergraph))

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--distance',  type="integer",            default=1e5,      help='Maximum distance for a linkage between a peak and a gene')
p$add_argument('--min_chip_threshold',  type="double",            default=0.25,      help='')
p$add_argument('--cor_threshold',  type="double",            default=0.30,      help='')
p$add_argument('--outdir',       type="character",                help='Output file')
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

## START TEST
args$distance <- 5e4
args$trajectory <- paste0(io$basedir,"/results/rna/trajectories/blood_trajectory/blood_trajectory.txt.gz")
args$trajectory_name <- "blood"
args$min_chip_threshold <- 0.25
args$cor_threshold <- 0.30
args$outdir <- paste0(io$basedir,"/results/rna_atac/gene_regulatory_networks/trajectories/blood")
## END TEST

# I/O
io$virtual_chip.dir <- file.path(io$basedir,"results_new/rna_atac/virtual_chipseq/CISBP")
io$virtual_chip.mtx <- file.path(io$basedir,"results_new/rna_atac/virtual_chipseq/CISBP/virtual_chip_matrix.rds")
# io$rna.sce <- file.path(io$basedir,"results_new/rna/trajectories/blood/blood_SingleCellExperiment.rds")

#####################
## Load trajectory ##
#####################

trajectory.dt <- fread(args$trajectory)

###################
## Load metadata ##
###################

sample_metadata <- fread(io$metadata) %>%
  .[pass_atacQC==TRUE & pass_rnaQC==TRUE & doublet_call==FALSE]

##################
## Filter cells ##
##################

cells <- intersect(trajectory.dt$cell,sample_metadata$cell)
trajectory.dt <- trajectory.dt[cell%in%cells] 
sample_metadata <- sample_metadata[cell%in%cells] %>% setkey(cell) %>% .[cells]

# opts$celltypes.subset <- opts$celltypes[opts$celltypes%in%unique(sample_metadata$celltype.predicted)]
opts$celltypes_trajectory <- c("Haematoendothelial_progenitors", "Blood_progenitors_1", "Blood_progenitors_2", "Erythroid1", "Erythroid2", "Erythroid3")

#######################################
## Load virtual ChIP-seq annotation ###
#######################################

io$virtual_chip.mtx <- file.path(io$basedir,"results_new/rna_atac/virtual_chipseq/CISBP/motifmatchr_virtual_chip.rds")
virtual_chip.mtx <- readRDS(io$virtual_chip.mtx) %>% assay()

virtual_chip.mtx[virtual_chip.mtx<=args$min_chip_threshold] <- 0

#####################################
## Load single-cell RNA expression ##
#####################################

rna.sce <- load_SingleCellExperiment(
  file = io$rna.sce, 
  cells = sample_metadata$cell, 
  normalise = TRUE, 
  remove_non_expressed_genes = FALSE
)

rna_tf.sce <- rna.sce[rownames(rna.sce)%in%str_to_title(colnames(virtual_chip.mtx)),]
rownames(rna_tf.sce) <- toupper(rownames(rna_tf.sce))

# rm(rna.sce)

####################################
## Load pseudobulk RNA expression ##
####################################

io$rna.pseudobulk.sce <- file.path(io$basedir,"results_new/rna/pseudobulk/SingleCellExperiment_CISBP_pseudobulk_celltype.mapped_mnn.rds")
rna_tf_pseudobulk.sce <- readRDS(io$rna.pseudobulk.sce)[,opts$celltypes_trajectory]
# rna_tf_pseudobulk.sce <- rna_pseudobulk.sce[str_to_title(rownames(rna_tf.sce))]
# rownames(rna_tf_pseudobulk.sce) <- toupper(rownames(rna_tf_pseudobulk.sce))

#############################
## Load peak2gene linkages ##
#############################

peak2gene.dt <- fread(io$archR.peak2gene.all) %>%
  .[gene%in%str_to_title(rownames(rna_tf.sce)) & dist<args$distance] %>%
  .[,gene:=toupper(gene)] %>%
  .[,peak:=sprintf("chr%s:%s-%s",chr,peak.start,peak.end)]

# peak2gene.dt <- peak2gene.dt[peak%in%peaks]

#################
## Smooth data ##
#################

pca.rna <- fread(io$pca.rna) %>% matrix.please %>% .[sample_metadata$cell,]
logcounts(rna_tf.sce) <- smoother_aggregate_nearest_nb(mat=as.matrix(logcounts(rna_tf.sce)), D=pdist(pca.rna), k=50)

#######################
## Feature selection ##
#######################

rna.var <- apply(logcounts(rna_tf.sce),1,var) %>% sort 
TFs.rna <- rna.var[rna.var>0.05] %>% names
# rna_tf_filt.sce <- rna_tf.sce[TFs.rna,]

# Using pseudobulk data
# tmp <- apply(logcounts(rna_tf_pseudobulk.sce),1,var) %>% sort %>% tail(n=50)
# rna_tf.sce <- rna_tf.sce[names(tmp),]
# foo <- data.table(celltype=colnames(rna_tf_pseudobulk.sce), logcounts(rna_tf_pseudobulk.sce)["GATA1",])

############################################################
## Select TFs with TF motif variability in the trajectory ##
############################################################

# chromvar_chip.mtx <- readRDS("/Users/argelagr/data/gastrulation_multiome_10x/results_new/atac/archR/chromvar_chip/pseudobulk/chromVAR_deviations_CISBP_archr_chip.rds")[,opts$celltypes_trajectory] %>% assay("z")
# chromvar_chip.var <- apply(chromvar_chip.mtx,1,var,na.rm=T) %>% sort
# TFs.atac <- chromvar_chip.var[chromvar_chip.var>=200] %>% names

chromvar_chip.mtx <- readRDS(file.path(io$basedir,"results_new/atac/archR/chromvar_chip/chromVAR_chip_deviations_CISBP_archr.rds"))[,sample_metadata$cell] %>% assay("z")
chromvar_chip.var <- apply(chromvar_chip.mtx,1,var,na.rm=T) %>% sort
TFs.atac <- chromvar_chip.var[chromvar_chip.var>=5] %>% names

#################
## Filter data ##
#################

TFs <- Reduce("intersect",list(TFs.rna,colnames(virtual_chip.mtx),TFs.atac))
# TFs <- c("TAL1","RUNX1","GATA1","KLF1","JUN")

rna_tf_filt.sce <- rna_tf.sce[TFs,]
rna_tf_pseudobulk_filt.sce <- rna_tf_pseudobulk.sce[TFs,]
virtual_chip_filt.mtx <- virtual_chip.mtx[,TFs]

peaks <- intersect(rownames(virtual_chip_filt.mtx), unique(peak2gene.dt[gene%in%TFs,peak]))
virtual_chip_filt.mtx <- virtual_chip_filt.mtx[peaks,]

######################################################
## Calculate RNA expression correlation between TFs ##
######################################################

corr.mtx <- psych::corr.test(t(logcounts(rna_tf_filt.sce)),t(logcounts(rna_tf_filt.sce)), ci = F)
diag(corr.mtx$r) <- NA

##############################
## Link TFs to target genes ##
##############################

# opts$pvalue.threshold <- 0.25

corr.mtx.filt <- corr.mtx
for (i in TFs) {
    
    # Select target peaks
    target_peaks_i <- names(which(virtual_chip_filt.mtx[,i]>=args$min_chip_threshold))
    
    if (length(target_peaks_i)>=1) {
      
      # Select target genes
      target_genes_i <- intersect(
        x = names(which(abs(corr.mtx$r[i,])>=args$cor_threshold)),  # all correlations
        # x = names(which(corr.mtx$r[i,]>=args$cor_threshold)),   # only positive correlations
        y = unique(peak2gene.dt[peak%in%target_peaks_i,gene])
      )
      
      print(sprintf("%s: %s target peaks and %s target genes",i,length(target_peaks_i),length(target_genes_i)))
      
      # Update correlation matrices
      corr.mtx.filt$p[i,!colnames(corr.mtx.filt$p)%in%target_genes_i] <- NA
      corr.mtx.filt$r[i,!colnames(corr.mtx.filt$p)%in%target_genes_i] <- NA
      
    } else {
      
      corr.mtx.filt$p[i,] <- NA
      corr.mtx.filt$r[i,] <- NA
  }
}

####################
## Create network ##
####################

tmp <- corr.mtx.filt$r
diag(tmp) <- NA

# Filter edges
# tmp[abs(tmp)<0.30] <- NA

# Filter nodes
# tmp <- tmp[rowSums(!is.na(tmp))>=1,]
# tmp <- tmp[,colSums(!is.na(tmp))>=1]
# hist(rowSums(!is.na(tmp)))
# hist(colSums(!is.na(tmp)))
# dim(tmp)

# Save
# saveRDS(tmp, sprintf("%s/corr_matrix.rds",args$outdir))

# Prepare data
node_list.dt <- data.table(node_id=1:nrow(tmp), node_name=rownames(tmp))
target_list.dt <- data.table(target_id=1:ncol(tmp), target_name=colnames(tmp))

edge_list.dt <- as.data.table(tmp,keep.rownames = T) %>% 
  setnames("rn","from") %>%
  melt(id.vars=c("from"), variable.name="to", value.name="weight") %>%
  .[!is.na(weight)]

# node_list_metadata.dt <- data.table(
#   label = c(node_list.dt$node_name, target_list.dt$target_name),
#   class = c(rep("TF",nrow(node_list.dt)), rep("gene",nrow(target_list.dt)))
# )


###########################
## Create network layout ##
###########################

# Create network
net <- graph_from_data_frame(d = edge_list.dt)
nw <- asNetwork(net)

# Define layout
# set.seed(42)
# layout <- layout_with_fr(net, weights = E(net)$weight)c
# layout <- layout.sphere(net)
# layout <- layout.mds(net)
node_coords.mtx <- sna::gplot.layout.fruchtermanreingold(nw, layout.par = NULL)

# Define node colours
node_colors <- opts$celltype.colors[colnames(rna_tf_pseudobulk_filt.sce)[apply(logcounts(rna_tf_pseudobulk_filt.sce[names(V(net)),]),1,which.max)]] %>% unname

# Define edge colors
E(net)[which(E(net)$weight<0)]$color <- "darkblue"  # Colour negative correlation edges as blue
E(net)[which(E(net)$weight>0)]$color <- "darkred"   # Colour positive correlation edges as red


# Plot
ggnet2(
  net = nw,
  mode = node_coords.mtx,
  color = node_colors,
  edge.color = E(net)$color,
  edge.size = 0.2,
  node.size = 5,
  # node.alpha = alphas,
  label = TRUE,
  label.color = "black",
  label.size = 2,
  arrow.gap = 0.025,
  arrow.size = 5,
  legend.position = "none"
) 


#################################
## Plot number of edges per TF ##
#################################

to.plot <- edge_list.dt %>%
  .[,sign:=sign(weight)] %>%
  .[,.N,by=c("from","sign")] %>% setnames("from","TF")

p <- ggdensity(to.plot, x="N", fill="gray70") +
  geom_vline(xintercept=mean(to.plot$N), linetype="dashed") +
  labs(x="TF connectivity") +
  theme_classic() +
  theme(
    axis.ticks.y = element_blank(),
    axis.text = element_text(size=rel(0.75), color="black")
  )

to.plot <- edge_list.dt[,.N,by="from"] %>% 
  setnames("from","TF") %>% 
  setorder(-N) %>% head(n=25)

to.plot[,TF:=factor(TF,levels=to.plot$TF)]

p <- ggplot(to.plot, aes_string(x="TF", y="N"), fill="gray70") +
  geom_point(size=2) +
  geom_segment(aes_string(xend="TF"), size=0.5, yend=0) +
  coord_flip() +
  labs(x="", y="Number of connections with other TFs") +
  theme_classic() +
  theme(
    axis.ticks.y = element_blank(),
    axis.text = element_text(size=rel(0.75), color="black")
  )

# pdf(sprintf("%s/%s/%s_%s_cobinding.pdf",io$outdir,i,i,j), width = 5, height = 4)
print(p)
# dev.off()

###########################################################################
## Plot network, colouring each gene by celltype-specific RNA expression ##
###########################################################################

# Load PAGA graph
source(here::here("load_paga_graph.R"))

# Define colors
pal <- grDevices::colorRamp(c("gray60", "purple"))( (1:100)/100 )
pal <- cbind(pal, seq(100, 255, length.out = 100))

rna_scaled.mtx <- apply(logcounts(rna_tf_pseudobulk_filt.sce),1,minmax.normalisation) %>% t

for (i in opts$celltypes_trajectory) {
  
  #####################
  ## Plot PAGA graph ##
  #####################
  
  colors <- rep("gray70",length(opts$celltypes)); names(colors) <- opts$celltypes
  colors[[i]] <- opts$celltype.colors[[i]]
  
  alphas <- rep(0.5,length(opts$celltypes)); names(alphas) <- opts$celltypes
  alphas[[i]] <- 1.0
  
  sizes <- rep(6,length(opts$celltypes)); names(sizes) <- opts$celltypes
  sizes[[i]] <- 13
  
  p1 <- ggnet2(
    net = net.paga,
    mode = c("x", "y"),
    color = colors,
    node.alpha = alphas,
    node.size = sizes,
    edge.size = 0.15,
    edge.color = "grey",
    label = TRUE,
    label.size = 3.5,
    legend.position = "none"
  )
  
  ##################
  ## Plot network ##
  ##################
  
  # Define node colors
  foo <- rna_scaled.mtx[,i]; names(foo) <- rownames(rna_scaled.mtx)
  rna.expr.i <- c(foo, rna_scaled.mtx[,i])[names(V(net))]# %>% minmax.normalisation
  
  # Discretise a bit
  # rna.expr.i[rna.expr.i<0.25] <- 0
  # rna.expr.i[rna.expr.i>0.75] <- 1
  
  # V(net)$color <- colourvalues::colour_values(rna.expr.i, palette = "viridis")
  V(net)$color <- colourvalues::colour_values(rna.expr.i, palette = pal)
  
  # Define edge colors
  E(net)[which(E(net)$weight<0)]$color <- "darkblue"  # Colour negative correlation edges as blue
  E(net)[which(E(net)$weight>0)]$color <- "darkred"   # Colour positive correlation edges as red
  
  # Define edge width
  edge_weight <- abs(E(net)$weight)   # Convert edge weights to absolute values
  for (j in seq_len(length(E(net)))) {
    # edge_weight[[j]] <- edge_weight[[j]]*(rna.expr.i[tail_of(net, j)] * rna.expr.i[head_of(net, j)])
    edge_weight[[j]] <- edge_weight[[j]]*(rna.expr.i[tail_of(net, j)])
  }
  edge_weight[edge_weight==0] <- 0.03
  edge_weight[edge_weight>0.75] <- 0.75
  
  alphas <- rna.expr.i
  # alphas[rna.expr.i<0.5] <- 0.25
  # alphas[rna.expr.i>0.75] <- 1
  
  p2 <- ggnet2(
    net = nw,
    mode = node_coords.mtx,
    color = V(net)$color,
    edge.color = E(net)$color,
    edge.size = edge_weight,
    node.size = 7,
    node.alpha = alphas,
    label = TRUE,
    label.color = "black",
    label.size = 3.5,
    arrow.gap = 0.025,
    arrow.size = 5,
    legend.position = "none"
  ) 
  
  ###################
  ## Combine plots ##
  ###################
  
  p <- cowplot::plot_grid(plotlist=list(p1,p2), rel_widths = c(1.5/5,3.5/5), nrow=1)
    
  outfile <- sprintf("%s/%d_%s.png",args$outdir,match(i,opts$celltypes_trajectory),i)
  png(outfile, width = 1200, height = 550, bg = "white")
  print(p)
  dev.off()
}



#######################################
## Plot network along the trajectory ##
#######################################

ntimepoints <- 25

trajectory.dt[,pseudotime_group:=as.numeric(cut(PC1,ntimepoints))]

rna_scaled.mtx <- apply(logcounts(rna_tf_filt.sce),1,minmax.normalisation) %>% t

for (i in unique(trajectory.dt$pseudotime_group)) {
  
  cells <- trajectory.dt[pseudotime_group==i,cell]
  rna.expr.i <- apply(rna_scaled.mtx[,cells],1,mean)[names(V(net))]

  #####################
  ## Plot pseudotime ##
  #####################

  gene.to.plot <- "GATA1"
  to.plot <- data.table(cell=colnames(rna_scaled.mtx), expr=rna_scaled.mtx[gene.to.plot,]) %>% 
    merge(sample_metadata[,c("cell","celltype.predicted")]) %>%
    merge(trajectory.dt[,c("cell","PC1")])
  
  tmp <- trajectory.dt[pseudotime_group==i,mean(PC1)]
  to.plot[,foo:=(PC1>(tmp-1) & PC1<(tmp+1))]
  
  p1 <- ggplot(to.plot, aes(x=PC1, y=expr, alpha=foo, size=foo)) +
    # ggrastr::geom_point_rast(aes(fill=celltype), size=1.25, shape=21, stroke=0.1) +
    geom_point(aes(fill=celltype.predicted), shape=21, stroke=0.1) +
    geom_rug(aes(color=celltype.predicted), sides="b") +
    stat_smooth(method="loess", color="black", alpha=0.75, span=0.5, size=0.5) +
    scale_color_manual(values=opts$celltype.colors) +
    scale_fill_manual(values=opts$celltype.colors) +
    scale_alpha_manual(values=c(0.05,1)) +
    scale_size_manual(values=c(1,2)) +
    guides(fill=F, color=F, size=F, alpha=F) +
    labs(x="Pseudotime", y=sprintf("%s expression",gene.to.plot)) +
    theme_classic() +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      legend.title = element_blank(),
      legend.position="top"
    )
  
  
  ##################
  ## Plot network ##
  ##################
  
  # V(net)$color <- colourvalues::colour_values(rna.expr.i, palette = "viridis")
  V(net)$color <- colourvalues::colour_values(rna.expr.i, palette = pal)
  
  # Define edge colors
  E(net)[which(E(net)$weight<0)]$color <- "darkblue"  # Colour negative correlation edges as blue
  E(net)[which(E(net)$weight>0)]$color <- "darkred"   # Colour positive correlation edges as red
  
  # Define edge width
  edge_weight <- abs(E(net)$weight)   # Convert edge weights to absolute values
  for (j in seq_len(length(E(net)))) {
    # edge_weight[[j]] <- edge_weight[[j]]*(rna.expr.i[tail_of(net, j)] * rna.expr.i[head_of(net, j)])
    edge_weight[[j]] <- edge_weight[[j]]*(rna.expr.i[tail_of(net, j)])
  }
  edge_weight[edge_weight==0] <- 0.03
  # edge_weight <- minmax.normalisation(edge_weight)
  edge_weight[edge_weight>0.75] <- 0.75
  
  alphas <- rna.expr.i
  # alphas[rna.expr.i<0.5] <- 0.25
  # alphas[rna.expr.i>0.75] <- 1
  
  p2 <- ggnet2(
    net = nw,
    mode = node_coords.mtx,
    color = V(net)$color,
    edge.color = E(net)$color,
    edge.size = edge_weight,
    node.size = 7,
    node.alpha = alphas,
    label = TRUE,
    label.color = "black",
    label.size = 3.5,
    arrow.gap = 0.025,
    arrow.size = 5,
    legend.position = "none"
  ) 
  
  ###################
  ## Combine plots ##
  ###################
  
  p <- cowplot::plot_grid(plotlist=list(p1,p2), nrow = 1, rel_widths = c(1/4,3/4))
  
  outfile <- sprintf("%s/gif/%d_blood_pseudotime_GRN.png",args$outdir,i)
  png(outfile, width = 1100, height = 350, bg = "white")
  print(p)
  dev.off()
}


################
## Create GIF ##
################

library(magick)

sprintf("%s/gif/%d_blood_pseudotime_GRN.png",args$outdir,unique(trajectory.dt$pseudotime_group)) %>%
  map(image_read) %>%
  image_join %>%
  image_animate(fps=5) %>%
  image_write(quality=100, path = sprintf("%s/gif/test.gif",args$outdir))


#############
## Explore ##
#############

sort(igraph::eigen_centrality(net)$vector)

total_edges <- igraph::degree(net) %>% sort

TFs.to.plot <- tail(total_edges,n=25) %>% names

A <- get.adjacency(net, attr="weight", sparse=F)[TFs.to.plot,]

to.plot <- data.table(
  tf = rownames(A),
  # degree_centrality = sort(igraph::degree(net)),
  nnegative = rowSums(A<0),
  npositive = rowSums(A>0)
) %>% .[,total:=nnegative+npositive] %>% melt(id.vars=c("tf","total")) %>% setorder(-total)

ggbarplot(to.plot, x="tf", y="value", fill="variable") +
  coord_flip()

hist(A["CLOCK",])

