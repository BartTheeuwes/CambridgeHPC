
source(here::here("settings.R"))
source(here::here("utils.R"))

suppressMessages(library(GGally))
suppressMessages(library(igraph))
suppressMessages(library(network))
suppressMessages(library(sna))
suppressMessages(library(intergraph))

#####################
## Define settings ##
#####################

io$basedir <- file.path(io$basedir,"test")
# args$rna_cells.sce <- io$rna.sce
io$rna_metacells.sce <- file.path(io$basedir, 'results/rna/metacells/trajectories/nmp/SingleCellExperiment_metacells.rds')
io$rna_pseudobulk.sce <- file.path(io$basedir,"results/rna/pseudobulk/celltype/SingleCellExperiment_pseudobulk.rds")
io$trajectory_file <- "/Users/argelagr/data/gastrulation_multiome_10x/results/rna/trajectories/nmp/nmp_trajectory.txt.gz"
io$trajectory <- "nmp"
io$grn_coef <- file.path(io$basedir,'results/rna_atac/gene_regulatory_networks/metacells/trajectories/nmp/global_chip_GRN_coef.txt.gz')
io$outdir <-  file.path(io$basedir,"results/rna_atac/gene_regulatory_networks/metacells/trajectories/nmp/fig"); dir.create(io$outdir, showWarnings = F)

# Options
opts$min_coef <- 0.25
opts$min_tf_score <- 0.75

if (io$trajectory=="nmp") {
  celltypes.to.plot <- c("Caudal_Mesoderm", "Somitic_mesoderm", "NMP", "Spinal_cord")
}

#####################
## Load trajectory ##
#####################

trajectory.dt <- fread(io$trajectory_file)

##############################
## Load marker gene and TFs ##
##############################

# TFs <- fread(paste0(io$basedir,"/results/rna_atac/gene_regulatory_networks/TFs.txt"))[[1]]
marker_TFs_all.dt <- fread(io$rna.atlas.marker_TFs.all) %>% 
  .[celltype%in%celltypes.to.plot] %>%
  .[,gene:=toupper(gene)]

marker_TFs_filt.dt <- marker_TFs_all.dt %>% 
  .[score>=opts$min_tf_score]

print(unique(marker_TFs_filt.dt$gene))

##############################
## Load RNA expression data ##
##############################

# SingleCellExperiment at cellular resolution
# rna_cells.sce <- readRDS(args$rna_cells.sce)[,trajectory.dt$cell]

# SingleCellExperiment at metacell resolution
rna_metacells.sce <- readRDS(io$rna_metacells.sce)

# SingleCellExperiment at pseudobulk resolution
rna_pseudobulk.sce <- readRDS(io$rna_pseudobulk.sce)#[,celltypes.to.plot]

##################################
## Load global GRN coefficients ##
##################################

GRN_coef.dt <- fread(io$grn_coef) %>% 
  .[,gene:=toupper(gene)] %>% .[gene%in%unique(marker_TFs_filt.dt$gene) & tf%in%unique(marker_TFs_filt.dt$gene)] %>%
  .[pvalue<0.10 & abs(beta)>=opts$min_coef]

##########################
## Filter TFs and genes ##
##########################

# Filter TFs that have few connections
# TFs <- intersect(GRN_coef.dt[,.N,by="tf"][N>=3,tf], GRN_coef.dt[,.N,by="gene"][N>=2,gene])
TFs <- union(GRN_coef.dt[,.N,by="tf"][N>=3,tf], GRN_coef.dt[,.N,by="gene"][N>=3,gene])
GRN_coef.dt <- GRN_coef.dt[tf%in%TFs & gene%in%TFs]

# Fetch RNA expression matrices
rna_tf.mtx <- logcounts(rna_pseudobulk.sce)[str_to_title(TFs),]; rownames(rna_tf.mtx) <- toupper(rownames(rna_tf.mtx))

# Scale
rna_tf_scaled.mtx <- apply(rna_tf.mtx,1,minmax.normalisation) %>% .[celltypes.to.plot,] %>% t 

####################
## Create network ##
####################

# Create node and edge data.frames
TFs <- unique(c(GRN_coef.dt$tf,GRN_coef.dt$gene))
node_list.dt <- data.table(node_id=1:length(TFs), node_name=TFs)
edge_list.dt <- GRN_coef.dt[,c("tf","gene","beta")] %>% setnames(c("from","to","weight")) %>% .[!from==to]

# Create igraph object
igraph.net <- graph_from_data_frame(d = edge_list.dt)

# igraph::degree(igraph.net)

###########################
## Define network layout ##
###########################

# Convert to network class
network.net <- asNetwork(igraph.net)

node_coords.mtx <- sna::gplot.layout.fruchtermanreingold(network.net, layout.par = NULL)

#########################
## Plot global network ##
#########################

marker_TFs_all.dt <- marker_TFs_all.dt[celltype!="Caudal_Mesoderm"]

# Define node color
tmp <- marker_TFs_all.dt[gene%in%names(V(igraph.net))] %>% .[,.SD[which.max(score)][,"celltype"],by="gene"] %>% setkey(gene)
# tmp[celltype=="Caudal_Mesoderm",celltype:="NMP"]
V(igraph.net)$color <- opts$celltype.colors[tmp[names(V(igraph.net)),celltype]]

# Define edge colors
# E(igraph.net)$color <- "gray70"  # Colour negative correlation edges as blue
E(igraph.net)[which(E(igraph.net)$weight<0)]$color <- "darkblue"  # Colour negative correlation edges as blue
E(igraph.net)[which(E(igraph.net)$weight>0)]$color <- "darkred"   # Colour positive correlation edges as red

# Define edge width
edge_weight <- abs(E(igraph.net)$weight)   # Convert edge weights to absolute values
edge_weight <- 0.60*(minmax.normalisation(edge_weight)+0.001)

p <- ggnet2(
  net = network.net,
  mode = node_coords.mtx, 
  color = V(igraph.net)$color,
  edge.color = E(igraph.net)$color,
  edge.size = edge_weight,
  # node.size = c(2.5,8)[factor(V(igraph.net)$class)],
  # node.alpha = alphas,
  # node.shape = V(igraph.net)$shape,
  label = TRUE,
  label.color = "black",
  # label.size = node.label.size,
  arrow.gap = 0.03,
  arrow.size = 3,
  legend.position = "none"
)

pdf(file.path(io$outdir,"global_network.pdf"), width = 10, height = 8)
print(p)
dev.off()



##################################
## Plot pepressive interactions ##
##################################

# Define edge width
edge_weight <- abs(E(igraph.net)$weight)
edge_weight[edge_weight>=1] <- 1
# edge_weight <- 0.60*(minmax.normalisation(edge_weight)+0.001)
edge_weight[E(igraph.net)$weight>0] <- 0.001

p <- ggnet2(
  net = network.net,
  mode = node_coords.mtx, 
  color = V(igraph.net)$color,
  edge.color = E(igraph.net)$color,
  edge.size = edge_weight,
  # node.size = c(2.5,8)[factor(V(igraph.net)$class)],
  # node.alpha = alphas,
  # node.shape = V(igraph.net)$shape,
  label = TRUE,
  label.color = "black",
  # label.size = node.label.size,
  arrow.gap = 0.03,
  arrow.size = 3,
  legend.position = "none"
)

pdf(file.path(io$outdir,"global_network_repressive_interactions.pdf"), width = 10, height = 8)
print(p)
dev.off()

##################################
## Plot activatory interactions ##
##################################

# Define edge width
edge_weight <- abs(E(igraph.net)$weight)
edge_weight[edge_weight>=1] <- 1
edge_weight <- 0.55*(minmax.normalisation(edge_weight)+0.001)
edge_weight[E(igraph.net)$weight<0] <- 0.001

p <- ggnet2(
  net = network.net,
  mode = node_coords.mtx, 
  color = V(igraph.net)$color,
  edge.color = E(igraph.net)$color,
  edge.size = edge_weight,
  # node.size = c(2.5,8)[factor(V(igraph.net)$class)],
  # node.alpha = alphas,
  # node.shape = V(igraph.net)$shape,
  label = TRUE,
  label.color = "black",
  # label.size = node.label.size,
  arrow.gap = 0.03,
  arrow.size = 3,
  legend.position = "none"
)

pdf(file.path(io$outdir,"global_network_activatory_interactions.pdf"), width = 10, height = 8)
print(p)
dev.off()

####################################################################
## Highlight same TF with both positive and negative interactions ##
####################################################################

tmp <- edge_list.dt %>% copy %>%
  setnames(c("tf","gene","cor")) %>%
  .[,cor:=sign(cor)*edge_weight] %>%
  .[,sign:=factor(c("-","+"))[(cor>0)+1]] %>%
  .[,.N,by=c("tf","sign")] %>% 
  dcast(tf~sign,value.var="N",fill=0)


i <- "CDX2"

igraph.net.test <- igraph.net
V(igraph.net.test)$color <- rep("gray95",length(V(igraph.net)))
V(igraph.net.test)$color[names(V(igraph.net))==i] <- "purple"
V(igraph.net.test)$color[names(V(igraph.net))%in%edge_list.dt[from==i,to]] <- "gray80"

# Define edge colors
E(igraph.net.test)[which(E(igraph.net.test)$weight<0)]$color <- "darkblue"  # Colour negative correlation edges as blue
E(igraph.net.test)[which(E(igraph.net.test)$weight>0)]$color <- "darkred"   # Colour positive correlation edges as red

# Define edge width
edge_weight <- abs(E(igraph.net.test)$weight)
edge_weight[edge_weight>=1] <- 1
for (j in seq_len(length(E(igraph.net.test)))) {
  if (names(tail_of(igraph.net.test,j))==i) {
    print(head_of(igraph.net.test,j))
    edge_weight[[j]] <- 0.8
  } else {
    edge_weight[[j]] <- 0.001
  }
}

p <- ggnet2(
  net = network.net,
  mode = node_coords.mtx, 
  color = V(igraph.net.test)$color,
  edge.color = E(igraph.net.test)$color,
  edge.size = edge_weight,
  label = TRUE,
  label.color = "black",
  arrow.gap = 0.03,
  arrow.size = 3,
  legend.position = "none"
)

pdf(file.path(io$outdir,"global_network_CDX2.pdf"), width = 10, height = 8)
print(p)
dev.off()

################################
## Plot network per cell type ##
################################

for (i in celltypes.to.plot) {
  
  # Define edge attributes
  # for (j in seq_len(length(E(igraph.net)))) {
  #   edge_weight[[j]] <- edge_weight[[j]]*(rna.expr.i[tail_of(igraph.net,j)])
  # }
  
  ############################
  ## Define node attributes ##
  ############################
  
  # Define node color palete
  pal <- grDevices::colorRamp(c("gray90", "purple"))( (1:100)/100 ); pal <- cbind(pal, seq(100, 255, length.out = 100))
  # pal <- colorRampPalette(c("red", "white", "blue"), space = "Lab")( 100 ); pal <- cbind(pal, seq(100, 255, length.out = 100))
  
  # Define node color based on marker strength
  marker_score.vec <- marker_TFs_all.dt$score; names(marker_score.vec) <- toupper(marker_TFs_all.dt$gene)
  marker_score.vec[marker_score.vec<=0.5] <- 0
  V(igraph.net)$color <- colourvalues::colour_values(c(0,marker_score.vec[names(V(igraph.net))]), palette = pal)[-1]
  
  # Load gene expression values per TF
  foo <- rna_tf_scaled.mtx[,i]; names(foo) <- rownames(rna_tf.mtx)
  rna.expr.i <- foo[names(V(igraph.net))]# %>% minmax.normalisation
  
  # Define node color based on expression levels
  V(igraph.net)$color <- colourvalues::colour_values(c(0,rna.expr.i[names(V(igraph.net))]), palette = pal)[-1]
  
  #####################
  ## Plot PAGA graph ##
  #####################
  
  source(here::here("load_paga_graph.R"))
  
  paga.colors <- rep("gray70",length(opts$celltypes)); names(paga.colors) <- opts$celltypes; paga.colors[[i]] <- opts$celltype.colors[[i]]
  paga.alphas <- rep(0.5,length(opts$celltypes)); names(paga.alphas) <- opts$celltypes; paga.alphas[[i]] <- 1.0
  paga.sizes <- rep(6,length(opts$celltypes)); names(paga.sizes) <- opts$celltypes; paga.sizes[[i]] <- 13
  
  p1 <- ggnet2(
    net = net.paga,
    mode = c("x", "y"),
    color = paga.colors,
    node.alpha = paga.alphas,
    node.size = paga.sizes,
    edge.size = 0.15,
    edge.color = "grey",
    label = TRUE,
    label.size = 3.5,
    legend.position = "none"
  )
  
  ##################
  ## Plot network ##
  ##################
  
  p2 <- ggnet2(
    net = network.net,
    mode = node_coords.mtx, 
    color = V(igraph.net)$color,
    edge.color = E(igraph.net)$color,
    edge.size = edge_weight,
    # node.size = c(2.5,8)[factor(V(igraph.net)$class)],
    # node.alpha = alphas,
    # node.shape = V(igraph.net)$shape,
    label = TRUE,
    label.color = "black",
    # label.size = node.label.size,
    arrow.gap = 0.015,
    arrow.size = 3,
    legend.position = "none"
  )
  
  # Save plots  
  # p <- cowplot::plot_grid(plotlist=list(p1,p2,p3,p4), rel_widths = c(2/5,3/5), rel_heights = c(2/3,1/3), nrow=2)
  # p <- cowplot::plot_grid(plotlist=list(p1,p2), rel_widths = c(2/5,3/5), nrow=1)
  pdf(sprintf("%s/%s_network.pdf",io$outdir,i), width = 10, height = 7)
  print(p2)
  dev.off()
  
  # Save network  
  # saveRDS(igraph.net, sprintf("%s/%s_network.rds",io$outdir,i))
}


############################################
## Plot number of correlated genes per TF ##
############################################



#########################
## Centrality measures ##
#########################

# igraph::degree(igraph.net)
# igraph::closeness(igraph.net)
# igraph::eigen_centrality(igraph.net)$vector
# igraph::betweenness(igraph.net)

############################################
## Plot number of correlated TFs per gene ##
############################################

# genes.to.plot <- tmp[,.(N=sum(abs(cor)>=0.15)),by="gene"] %>% setorder(-N) %>% head(n=15) %>% .$gene
# to.plot <- tmp %>%
#   .[,.(N=sum(abs(cor)>=0.15)),by=c("gene","sign")] %>%
#   .[gene%in%genes.to.plot] %>% .[,gene:=factor(gene,levels=genes.to.plot)]
# 
# p4 <- ggbarplot(to.plot, x="gene", y="N", fill="sign", width=0.55) +
#   coord_flip() +
#   scale_fill_manual(values=c("-"="blue", "+"="red"), drop=F) +
#   labs(x="", y="Number of correlated TFs") +
#   theme(
#     axis.text.y = element_text(size=rel(0.75)),
#     axis.text.x = element_text(colour="black",size=rel(0.75)),
#     axis.ticks.x = element_line(size=rel(0.75)),
#     legend.position = "top"
#   )


##########
## TEST ##
##########

####################
## Create network ##
####################

# Create igraph object
igraph_t.net <- graph_from_data_frame(d = edge_list.dt[from=="T" | to=="T"])

# Convert to network class
network_t.net <- asNetwork(igraph_t.net)

node_coords_t.mtx <- sna::gplot.layout.fruchtermanreingold(network_t.net, layout.par = NULL)

# Define node color
# V(igraph_t.net)$color <- ifelse(names(V(igraph_t.net))=="T","gray30","gray90")
tmp <- marker_TFs_all.dt[gene%in%names(V(igraph_t.net))] %>% .[,.SD[which.max(score)][,"celltype"],by="gene"] %>% setkey(gene)
V(igraph_t.net)$color <- opts$celltype.colors[tmp[names(V(igraph_t.net)),celltype]]

# Define edge colors
# E(igraph_t.net)$color <- "gray70"  # Colour negative correlation edges as blue

# Define edge width
edge_weight <- abs(E(igraph_t.net)$weight)   # Convert edge weights to absolute values
edge_weight <- 0.90*(minmax.normalisation(edge_weight)+0.001)
for (j in seq_len(length(E(igraph_t.net)))) {
  edge_weight[[j]] <- ifelse(names(tail_of(igraph_t.net,j))=="T",0.5,0.25)
}

# Define edge colors
# E(igraph.net)$color <- "gray70"  # Colour negative correlation edges as blue
E(igraph_t.net)[which(E(igraph_t.net)$weight<0)]$color <- "darkblue"  # Colour negative correlation edges as blue
E(igraph_t.net)[which(E(igraph_t.net)$weight>0)]$color <- "darkred"   # Colour positive correlation edges as red

# Define edge alpha
# edge_alpha <- rep(0.1,length(edge_weight))
# for (j in seq_len(length(E(igraph_t.net)))) {
#   edge_alpha[[j]] <- ifelse(names(tail_of(igraph_t.net,j))=="T",1,0.05)
# }

ggnet2(
  net = network_t.net,
  mode = node_coords_t.mtx, 
  color = V(igraph_t.net)$color,
  size = 15,
  edge.color = E(igraph_t.net)$color,
  edge.size = edge_weight,
  # edge.alpha = edge_alpha,
  # node.size = c(2.5,8)[factor(V(igraph_t.net)$class)],
  # node.alpha = alphas,
  # node.shape = V(igraph_t.net)$shape,
  label = TRUE,
  label.color = "black",
  # label.size = node.label.size,
  arrow.gap = 0.05,
  arrow.size = 8,
  legend.position = "none"
)

pdf(file.path(io$outdir,"global_network_highlight_T.pdf"), width = 13, height = 6.5)
print(p)
dev.off()
