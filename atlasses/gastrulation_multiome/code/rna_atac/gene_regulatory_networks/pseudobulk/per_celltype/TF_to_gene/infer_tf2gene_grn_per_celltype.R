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
p$add_argument('--celltype',    type="character",                help='celltype to obtain marker genes from')
p$add_argument('--min_coef',  type="double",            default=0.25,      help='Minimum coefficient from the linear regression (absolute values)')
p$add_argument('--min_tf_marker_strength',  type="double",            default=0.25,      help='Minimum TF marker strength')
p$add_argument('--min_gene_marker_strength',  type="double",            default=0.75,      help='Minimum gene marker strength')
p$add_argument('--outdir',       type="character",               help='Output file')
args <- p$parse_args(commandArgs(TRUE))

## START TEST
args$celltype <- c("NMP")
args$min_coef <- 0.25
args$min_tf_marker_strength <- 0.75
args$min_gene_marker_strength <- 0.25
args$outdir <- file.path(io$basedir,"/results/rna_atac/gene_regulatory_networks/pseudobulk/per_celltype/Gut")
## END TEST

##############################
## Load marker gene and TFs ##
##############################

# TFs <- fread(paste0(io$basedir,"/results/rna_atac/gene_regulatory_networks/TFs.txt"))[[1]]
marker_TFs_all.dt <- fread(io$rna.atlas.marker_TFs.all) %>% .[celltype%in%args$celltype]
marker_TFs.dt <- marker_TFs_all.dt %>% .[score>=args$min_tf_marker_strength]

marker_genes_all.dt <- fread(io$rna.atlas.marker_genes.all) %>%.[celltype%in%args$celltype]# %>% .[!gene%in%unique(marker_TFs.dt$gene) & celltype%in%args$celltype]
marker_genes.dt <- marker_genes_all.dt %>% .[score>=args$min_gene_marker_strength]

##############################
## Load pseudobulk RNA data ##
##############################

io$rna.pseudobulk.sce <- file.path(io$basedir,"results_new/rna/pseudobulk/SingleCellExperiment_pseudobulk_celltype.mapped_mnn.rds")
sce <- readRDS(io$rna.pseudobulk.sce)#[,opts$celltypes]

##################################
## Load global GRN coefficients ##
##################################

io$grn_coef <- file.path(io$basedir,'results_new/rna_atac/gene_regulatory_networks/pseudobulk/test/global_chip_GRN_coef.txt.gz')

GRN_coef.dt <- fread(io$grn_coef) %>% 
  .[pvalue<0.10 & abs(beta)>=args$min_coef]

# Remove T as a target gene because it causes issues
GRN_coef.dt <- GRN_coef.dt[gene!="T"]

##########################
## Filter TFs and genes ##
##########################

# Subset to TFs and genes that are markers of the cell type of interest
TFs <- Reduce("intersect",list(toupper(marker_TFs.dt$gene),unique(GRN_coef.dt$tf),toupper(rownames(sce))))
genes <- Reduce("intersect",list(marker_genes.dt$gene,unique(GRN_coef.dt$gene),rownames(sce)))

GRN_coef.dt <- GRN_coef.dt[tf%in%TFs & gene%in%genes]

# Fetch RNA expression matrices
rna_tf.mtx <- logcounts(sce)[str_to_title(TFs),]; rownames(rna_tf.mtx) <- toupper(rownames(rna_tf.mtx))
rna_genes.mtx <- logcounts(sce)[genes,]

# Scale
rna_tf_scaled.mtx <- apply(rna_tf.mtx,1,minmax.normalisation) %>% t
rna_genes_scaled.mtx <- apply(rna_genes.mtx,1,minmax.normalisation) %>% t

###################################
## Create GRN coefficient matrix ##
###################################

GRN_coef.mtx <- GRN_coef.dt %>%
  dcast(tf~gene, value.var="beta", fill=0) %>%
  matrix.please

####################
## Create network ##
####################

# Create node and edge data.frames
node_list.dt <- data.table(node_id=1:nrow(GRN_coef.mtx), node_name=rownames(GRN_coef.mtx))
target_list.dt <- data.table(target_id=1:ncol(GRN_coef.mtx), target_name=colnames(GRN_coef.mtx))
edge_list.dt <- GRN_coef.dt[,c("tf","gene","beta")] %>% setnames(c("from","to","weight"))
node_list_metadata.dt <- data.table(
  label = c(node_list.dt$node_name, target_list.dt$target_name),
  class = c(rep("TF",nrow(node_list.dt)), rep("gene",nrow(target_list.dt)))
)

# Create igraph object
igraph.net <- graph_from_data_frame(d = edge_list.dt, vertices = node_list_metadata.dt)

############################
## Define node attributes ##
############################

# Define groups
V(igraph.net)$group <- factor(V(igraph.net)$class, levels=c("TF","gene"))

# Define node shapes
V(igraph.net)$shape <- stringr::str_replace_all(V(igraph.net)$class,c("gene"="circle","TF"="triangle"))

# Define node label size
node.label.size <- c(1.5,2.5)[factor(V(igraph.net)$class)]
# node.label.size[rna.expr.i<0.05] <- 0

# Add gene expression values as node attributes
i <- args$celltype
igraph.net <- set_vertex_attr(igraph.net, name = "expr", value = c(rna_tf.mtx[,i],rna_genes.mtx[,i])[names(V(igraph.net))])

# Define node colors
# pal <- grDevices::colorRamp(c("gray60", "purple"))( (1:100)/100 ); pal <- cbind(pal, seq(100, 255, length.out = 100))
pal.TFs <- grDevices::colorRamp(c("gray60", "purple"))( (1:100)/100 ); pal.TFs <- cbind(pal.TFs, seq(100, 255, length.out = 100))
pal.genes <- grDevices::colorRamp(c("gray60", "darkgreen"))( (1:100)/100 ); pal.genes <- cbind(pal.genes, seq(100, 255, length.out = 100))

# Load gene expression values per TF
foo <- rna_tf_scaled.mtx[,i]; names(foo) <- rownames(rna_tf.mtx)
bar <- rna_genes_scaled.mtx[,i]; names(bar) <- rownames(rna_genes_scaled.mtx)
rna.expr.i <- c(foo, bar)[names(V(igraph.net))]# %>% minmax.normalisation

# V(igraph.net)$color <- c(
#   colourvalues::colour_values(rna.expr.i[names(V(igraph.net))[V(igraph.net)$class=="TF"]], palette = pal.TFs), 
#   colourvalues::colour_values(rna.expr.i[names(V(igraph.net))[V(igraph.net)$class=="gene"]], palette = pal.genes)
# )

# Define node color based on marker strength
marker_strength.tf <- marker_TFs_all.dt$score; names(marker_strength.tf) <- toupper(marker_TFs_all.dt$gene)
marker_strength.genes <- marker_genes_all.dt$score; names(marker_strength.genes) <- marker_genes_all.dt$gene

V(igraph.net)$color <- c(
  colourvalues::colour_values(c(0,marker_strength.tf[names(V(igraph.net))[V(igraph.net)$class=="TF"]]), palette = pal.TFs)[-1], 
  colourvalues::colour_values(c(0,marker_strength.genes[names(V(igraph.net))[V(igraph.net)$class=="gene"]]), palette = pal.genes)[-1]
)

############################
## Define edge attributes ##
############################

# Define edge colors
# E(igraph.net)$color <- "gray70"  # Colour negative correlation edges as blue
E(igraph.net)[which(E(igraph.net)$weight<0)]$color <- "darkblue"  # Colour negative correlation edges as blue
E(igraph.net)[which(E(igraph.net)$weight>0)]$color <- "darkred"   # Colour positive correlation edges as red

# Define edge width
edge_weight <- abs(E(igraph.net)$weight)   # Convert edge weights to absolute values
for (j in seq_len(length(E(igraph.net)))) {
  # edge_weight[[j]] <- edge_weight[[j]]*(rna.expr.i[tail_of(igraph.net, j)] * rna.expr.i[head_of(igraph.net, j)])
  edge_weight[[j]] <- edge_weight[[j]]*(rna.expr.i[tail_of(igraph.net, j)])
}
# edge_weight[edge_weight==0] <- 0.01
# edge_weight <- minmax.normalisation(edge_weight)

alphas <- rna.expr.i
# alphas[rna.expr.i<0.5] <- 0.25
# alphas[rna.expr.i>0.75] <- 1

###########################
## Define network layout ##
###########################

# Convert to network class
network.net <- asNetwork(igraph.net)

node_coords.mtx <- sna::gplot.layout.fruchtermanreingold(network.net, layout.par = NULL)

########################################################
## Plot PAGA graph to highlight cell type of interest ##
########################################################

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
  mode = node_coords.mtx, 
  net = network.net,
  color = V(igraph.net)$color,
  edge.color = E(igraph.net)$color,
  edge.size = edge_weight,
  node.size = c(2.5,8)[factor(V(igraph.net)$class)],
  # node.alpha = alphas,
  node.shape = V(igraph.net)$shape,
  label = TRUE,
  label.color = "black",
  label.size = node.label.size,
  arrow.gap = 0.01,
  arrow.size = 3,
  legend.position = "none"
)

############################################
## Plot number of correlated genes per TF ##
############################################

tmp <- edge_list.dt %>% copy %>%
  setnames(c("TF","gene","cor")) %>%
  .[,cor:=sign(cor)*edge_weight] %>%
  .[,sign:=factor(c("-","+"))[(cor>0)+1]]

TFs.to.plot <- tmp[,.(N=sum(abs(cor)>=0.15)),by="TF"] %>% setorder(-N) %>% head(n=15) %>% .$TF
to.plot <- tmp %>%
  .[,.(N=sum(abs(cor)>=0.15)),by=c("TF","sign")] %>%
  .[TF%in%TFs.to.plot] %>% .[,TF:=factor(TF,levels=TFs.to.plot)]

p3 <- ggbarplot(to.plot, x="TF", y="N", fill="sign", width=0.55) +
  coord_flip() +
  scale_fill_manual(values=c("-"="blue", "+"="red"), drop=F) +
  labs(x="", y="Number of correlated genes") +
  theme(
    axis.text.y = element_text(size=rel(0.75)),
    axis.text.x = element_text(colour="black",size=rel(0.75)),
    axis.ticks.x = element_line(size=rel(0.75)),
    legend.position = "top"
  )

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

###################
## Combine plots ##
###################

# p <- cowplot::plot_grid(plotlist=list(p1,p2,p3,p4), rel_widths = c(2/5,3/5), rel_heights = c(2/3,1/3), nrow=2)
p <- cowplot::plot_grid(plotlist=list(p1,p2), rel_widths = c(2/5,3/5), nrow=1)

# png(sprintf("%s/%s_network.png",args$outdir,i), width = 1200, height = 850, bg = "white")
pdf(sprintf("%s/%s_network.pdf",args$outdir,i), width = 14, height = 9)
print(p)
dev.off()


##################
## Save network ##
##################

saveRDS(igraph.net, sprintf("%s/%s_network.rds",args$outdir,i))
