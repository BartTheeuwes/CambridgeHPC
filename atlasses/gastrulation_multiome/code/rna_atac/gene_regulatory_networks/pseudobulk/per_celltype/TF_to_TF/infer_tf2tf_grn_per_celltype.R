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
# p$add_argument('--min_gene_marker_strength',  type="double",            default=0.75,      help='Minimum gene marker strength')
p$add_argument('--outdir',       type="character",               help='Output file')
args <- p$parse_args(commandArgs(TRUE))

## START TEST
args$celltype <- c("Neural_crest")
args$min_coef <- 0.25
args$min_tf_marker_strength <- 0.80
# args$min_gene_marker_strength <- 0.25
args$outdir <- file.path(io$basedir,"/results/rna_atac/gene_regulatory_networks/pseudobulk/per_celltype/TF_to_TF/Neural_crest")
## END TEST

# I/O
dir.create(args$outdir, showWarnings = F)

##############################
## Load marker gene and TFs ##
##############################

# TFs <- fread(paste0(io$basedir,"/results/rna_atac/gene_regulatory_networks/TFs.txt"))[[1]]
marker_TFs_all.dt <- fread(io$rna.atlas.marker_TFs.all) %>% .[celltype%in%args$celltype]
marker_TFs.dt <- marker_TFs_all.dt %>% .[score>=args$min_tf_marker_strength]

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
  .[,gene:=toupper(gene)] %>%
  .[pvalue<0.10 & abs(beta)>=args$min_coef]

# Remove T as a target gene because it causes issues
GRN_coef.dt <- GRN_coef.dt[gene!="T"]

##########################
## Filter TFs and genes ##
##########################

# Select TFs that are markers of other ectodermal cell types
# target_TFs <- fread(io$rna.atlas.marker_TFs.up) %>% 
#   .[!celltype%in%args$celltype] %>%
#   .[,gene:=toupper(gene)] %>% .$gene %>% unique

# Subset to TFs and genes that are markers of the cell type of interest
TFs <- Reduce("intersect",list(toupper(marker_TFs.dt$gene),unique(GRN_coef.dt$tf),toupper(rownames(sce))))
# GRN_coef.dt <- GRN_coef.dt[tf%in%TFs & gene%in%c(TFs,target_TFs)]
GRN_coef.dt <- GRN_coef.dt[tf%in%TFs & gene%in%TFs]

# Filter TFs that have no output connections
# TFs <- GRN_coef.dt[,.N,by="tf"] %>% .[N>=1,tf]
# GRN_coef.dt <- GRN_coef.dt[tf%in%TFs]

# Fetch RNA expression matrices
rna_tf.mtx <- logcounts(sce)[str_to_title(TFs),]; rownames(rna_tf.mtx) <- toupper(rownames(rna_tf.mtx))

# Scale
rna_tf_scaled.mtx <- apply(rna_tf.mtx,1,minmax.normalisation) %>% t

##########################
## Create network (OLD) ##
##########################

# # Create node and edge data.frames
# # TFs <- unique(c(GRN_coef.dt$tf,GRN_coef.dt$gene))
# node_list.dt <- data.table(node_id=1:length(unique(GRN_coef.dt$tf)), node_name=unique(GRN_coef.dt$tf))
# target_list.dt <- data.table(target_id=1:length(unique(GRN_coef.dt$gene)), target_name=unique(GRN_coef.dt$gene))
# edge_list.dt <- GRN_coef.dt[,c("tf","gene","beta")] %>% 
#   setnames(c("from","to","weight")) %>%
#   .[!from==to]
# 
# node_list_metadata.dt <- data.table(
#   label = c(unique(GRN_coef.dt$tf), unique(GRN_coef.dt$gene)),
#   class = c(rep("TF",length(unique(GRN_coef.dt$tf))), rep("TF_target",length(unique(GRN_coef.dt$gene))))
# ) %>% .[!(label%in%TFs & class=="TF_target")]
# 
# # stopifnot(node_list_metadata.dt[class=="TF_target",label]%in%target_TFs)
# # stopifnot(node_list_metadata.dt[class=="TF",label]%in%TFs)
# # node_list_metadata.dt[duplicated(node_list_metadata.dt$label)] %>% View
# # stopifnot(unique(edge_list.dt$from) %in% node_list.dt$node_name)
# # stopifnot(unique(edge_list.dt$to) %in% target_list.dt$target_name)
# # stopifnot(unique(edge_list.dt$to) %in% node_list_metadata.dt$label)
# 
# # Create igraph object
# igraph.net <- graph_from_data_frame(d = edge_list.dt, vertices = node_list_metadata.dt)

# V(igraph.net)$color <- colourvalues::colour_values(marker_strength.vec[names(V(igraph.net))], palette = pal)

####################
## Create network ##
####################

# Create node and edge data.frames
TFs <- unique(c(GRN_coef.dt$tf,GRN_coef.dt$gene))
node_list.dt <- data.table(node_id=1:length(TFs), node_name=TFs)
edge_list.dt <- GRN_coef.dt[,c("tf","gene","beta")] %>% setnames(c("from","to","weight")) %>% .[!from==to]

# Create igraph object
igraph.net <- graph_from_data_frame(d = edge_list.dt)

############################
## Define node attributes ##
############################

# Define node shapes
# V(igraph.net)$shape <- stringr::str_replace_all(V(igraph.net)$class,c("gene"="circle","TF"="triangle"))

# Define node label size
# node.label.size <- c(1.5,2.5)[factor(V(igraph.net)$class)]
# node.label.size[rna.expr.i<0.05] <- 0

# Define node colors
pal <- grDevices::colorRamp(c("gray60", "purple"))( (1:100)/100 ); pal <- cbind(pal, seq(100, 255, length.out = 100))

# Define node color based on marker strength
# marker_strength.vec <- marker_TFs_all.dt$score; names(marker_strength.vec) <- toupper(marker_TFs_all.dt$gene)
# marker_strength.vec[marker_strength.vec<=0.5] <- 0
# V(igraph.net)$color <- colourvalues::colour_values(c(0,marker_strength.vec[names(V(igraph.net))]), palette = pal)[-1]

# Load gene expression values per TF
# foo <- rna_tf_scaled.mtx[,args$celltype]; names(foo) <- rownames(rna_tf.mtx)
foo <- rna_tf.mtx[,args$celltype]; names(foo) <- rownames(rna_tf.mtx)
rna.expr.i <- foo[names(V(igraph.net))]# %>% minmax.normalisation

# Define node color based on expression levels
V(igraph.net)$color <- colourvalues::colour_values(c(0,rna.expr.i[names(V(igraph.net))]), palette = pal)[-1]

############################
## Define edge attributes ##
############################

# Define edge colors
E(igraph.net)$color <- "gray70"  # Colour negative correlation edges as blue
# E(igraph.net)[which(E(igraph.net)$weight<0)]$color <- "darkblue"  # Colour negative correlation edges as blue
# E(igraph.net)[which(E(igraph.net)$weight>0)]$color <- "darkred"   # Colour positive correlation edges as red

# Define edge width
edge_weight <- abs(E(igraph.net)$weight)   # Convert edge weights to absolute values
for (j in seq_len(length(E(igraph.net)))) {
  edge_weight[[j]] <- edge_weight[[j]]*(rna.expr.i[tail_of(igraph.net, j)])
}
edge_weight <- minmax.normalisation(edge_weight)+0.01

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

i <- args$celltype
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
  # node.size = c(2.5,8)[factor(V(igraph.net)$class)],
  # node.alpha = alphas,
  # node.shape = V(igraph.net)$shape,
  label = TRUE,
  label.color = "black",
  # label.size = node.label.size,
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
