suppressMessages(library(GGally))
suppressMessages(library(igraph))
suppressMessages(library(network))
suppressMessages(library(sna))
suppressMessages(library(intergraph))

#####################
## Define settings ##
#####################

source(here::here("settings.R"))
source(here::here("utils.R"))

# I/O
io$rna.sce <- file.path(io$basedir, 'results/rna/metacells/all_cells/SingleCellExperiment_metacells.rds')
io$grn_coef <- file.path(io$basedir,'results/rna_atac/gene_regulatory_networks/metacells/all_cells/global_chip_GRN_coef.txt.gz')
io$outdir <-  file.path(io$basedir,"results/rna_atac/gene_regulatory_networks/metacells/all_cells/test"); dir.create(io$outdir, showWarnings = F)

# Options
# opts$celltypes <- setdiff(opts$celltypes, c("Visceral_endoderm","ExE_endoderm","ExE_ectoderm","Parietal_endoderm"))

###########################
## Load GRN coefficients ##
###########################

GRN_coef.dt <- fread(io$grn_coef)

# Remove T as a target gene because it causes issues
GRN_coef.dt <- GRN_coef.dt[gene!="T"]

##############################
## Load RNA expression data ##
##############################

sce <- readRDS(io$rna.sce)

table(sce$celltype)

##########################
## Filter TFs and genes ##
##########################

TFs <- intersect(unique(GRN_coef.dt$tf),toupper(rownames(sce)))
genes <- intersect(unique(GRN_coef.dt$gene),rownames(sce))

# GRN_coef.dt <- GRN_coef.dt[tf%in%TFs & gene%in%genes,]

# Fetch RNA expression matrices
rna_tf.mtx <- logcounts(sce)[str_to_title(unique(GRN_coef.dt$tf)),]; rownames(rna_tf.mtx) <- toupper(rownames(rna_tf.mtx))
rna_targets.mtx <- logcounts(sce)[unique(GRN_coef.dt$gene),]

###################################
## Create GRN coefficient matrix ##
###################################

# remove negative links
# GRN_coef.dt <-  GRN_coef.dt[beta>0]

# Create matrix of positive regression coefficients
GRN_coef.mtx <- GRN_coef.dt %>% 
  dcast(tf~gene, value.var="beta", fill=0) %>%
  matrix.please

print(sum(GRN_coef.mtx>0))

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

# Create networks
# igraph.net <- graph_from_data_frame(d = edge_list.dt)
igraph.net <- graph_from_data_frame(d = edge_list.dt, vertices = node_list_metadata.dt)

# Define groups
V(igraph.net)$group <- factor(V(igraph.net)$class, levels=c("TF","gene"))

# Define node shapes
V(igraph.net)$shape <- stringr::str_replace_all(V(igraph.net)$class,c("gene"="circle","TF"="triangle"))

# Define node color
# V(net)$color <- c(
#   colourvalues::colour_values(rna.expr.i[names(V(net))[V(net)$class=="TF"]], palette = pal.TFs), 
#   colourvalues::colour_values(rna.expr.i[names(V(net))[V(net)$class=="gene"]], palette = pal.genes)
# )
V(igraph.net)$color <- stringr::str_replace_all(V(igraph.net)$class,c("gene"="green","TF"="purple"))

# Define edge colors
E(igraph.net)[which(E(igraph.net)$weight<0)]$color <- "darkblue"  # Colour negative correlation edges as blue
E(igraph.net)[which(E(igraph.net)$weight>0)]$color <- "darkred"   # Colour positive correlation edges as re
E(igraph.net)[which(E(igraph.net)$weight==0)]$color <- "gray"   # Colour positive correlation edges as re

##################
## Plot network ##
##################

# Convert to network class
network.net <- asNetwork(igraph.net)

# Define layout
node_coords.mtx <- sna::gplot.layout.fruchtermanreingold(network.net, layout.par = NULL)

# Define node colours
# node_colors <- opts$celltype.colors[colnames(rna_tf_pseudobulk_filt.sce)[apply(logcounts(rna_tf_pseudobulk_filt.sce[names(V(net)),]),1,which.max)]] %>% unname

# Define node label size
node.label.size <- c(3,5)[factor(V(igraph.net)$class)]

# Define edge weight
edge_weight <- abs(E(igraph.net)$weight)   # Convert edge weights to absolute values

ggnet2(
  net = network.net,
  mode = node_coords.mtx,
  color = V(igraph.net)$color,
  edge.color = E(igraph.net)$color,
  # edge.size = edge_weight,
  node.size = c(0.5,1)[factor(V(igraph.net)$class)],
  label = FALSE,
  label.color = "black",
  label.size = node.label.size,
  arrow.gap = 0.025,
  arrow.size = 2.5,
  legend.position = "none"
) 

############################
## simulate perturbations ##
############################

# gene_KO <- 'Tal1'
# delta_X <- KO_simulate(
#   gene_KO = gene_KO,
#   GRN_coef.mtx = GRN_coef.mtx,
#   rna.mtx = rna.mtx, 
#   n_propagation = 2
# )

n_propagation <- 1
gene_KO <- "TAL1"
stopifnot(gene_KO%in%rownames(rna_tf.mtx))

delta_X_TFs.mtx <- matrix(0, nrow=ncol(rna_tf.mtx), ncol=nrow(rna_tf.mtx))
dimnames(delta_X_TFs.mtx) <- dimnames(t(rna_tf.mtx))
delta_X_TFs.mtx[,gene_KO] <- - rna_tf.mtx[gene_KO,]

delta_X_genes.mtx <- delta_X_TFs.mtx %*% GRN_coef.mtx

# hist(delta_X_genes.mtx[delta_X_genes.mtx!=0])

###############
## PAGA plot ##
###############

# TO-DO: ALLOW FOR POSITIVE AND NEGATIVE VALUES

# gene2plot <- GRN_coef.dt[tf=="TAL1"] %>% sort.abs("beta") %>% head(n=1) %>% .$gene
gene2plot <- GRN_coef.dt[tf=="TAL1"] %>% sort.abs("beta") %>% tail(n=1) %>% .$gene

# Define colormap for the perturbation
opts$min.delta_X <- round(min(delta_X_genes.mtx))
delta_X.col.seq <- round(seq(opts$min.delta_X,0,1), 0)
perturbation.colors <- colorRampPalette(c("red", "gray90"))(length(delta_X.col.seq))

alphas <- rep(0.75,length(opts$celltypes)); names(alphas) <- opts$celltypes
sizes <- rep(6,length(opts$celltypes)); names(sizes) <- opts$celltypes

source(here::here("load_paga_graph.R"))

delta_X.values <- delta_X_genes.mtx[,gene2plot]
delta_X.values[delta_X.values<opts$min.delta_X] <- opts$min.delta_X
delta_X.colors <- round(delta_X.values,0) %>% map(~ perturbation.colors[which(delta_X.col.seq == .)]) %>% unlist

p.paga <- ggnet2(
  net = net.paga,
  mode = c("x", "y"),
  node.size = 0,
  edge.size = 0.15,
  edge.color = "grey",
  label = FALSE,
  label.size = 2.3
)

p1 <- p.paga + geom_text(label = "\u25D0", aes(x=x, y=y), color=delta_X.colors[p.paga$data$label], size=16, family = "Arial Unicode MS") +
  scale_colour_manual(values=delta_X.colors) + 
  labs(title=paste(gene2plot,'change after',gene_KO,'KO') )+
  theme(
    plot.title = element_text(hjust = 0.5)
  )

p2 <- ggnet2(
  net = net.paga,
  mode = c("x", "y"),
  color = opts$celltype.colors[opts$celltypes],
  node.alpha = alphas,
  node.size = sizes,
  edge.size = 0.15,
  edge.color = "grey",
  label = TRUE,
  label.size = 3.5,
  legend.position = "none"
)

cowplot::plot_grid(plotlist=list(p1,p2), nrow=1)


#################
## Exploration ##
#################

tf.to.plot <- "NKX2-5"

# genes.to.plot <- GRN_coef.dt[tf=="TAL1",gene]
genes.to.plot <- GRN_coef.dt[tf==tf.to.plot & pvalue<=0.01 & beta>0,gene]
# View(GRN_coef.dt[tf=="NKX2-5" & pvalue<=0.05 & beta>0])

## Stacked barplots per gene ##
to.plot <- apply(rna_targets.mtx[genes.to.plot,], 1, function(i) i/sum(i)) %>%
  as.data.table(keep.rownames = T) %>% setnames("rn","celltype") %>%
  melt(id.vars="celltype", variable.name="gene")

# to.plot[gene=="Usp2",c("celltype","value")] %>% setorder(-value) %>% print
# sort(rna_targets.mtx["Gamt",])

p <- ggplot(to.plot, aes(x=gene, y=value)) +
  geom_bar(aes(fill=celltype), stat="identity", color="black", position="fill") +
  # facet_wrap(~cor_sign, nrow=2, scales="free_x") +
  scale_fill_manual(values=opts$celltype.colors[names(opts$celltype.colors)%in%unique(to.plot$celltype)]) +
  theme_classic() +
  labs(x="", y="") +
  # guides(x = guide_axis(angle = 90)) +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    axis.text.x = element_blank(),
    # axis.text.x = element_text(color="black", size=rel(0.85)),
    # axis.text.y = element_text(color="black", size=rel(1.0)),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank()
  )

# pdf(sprintf("%s/TF_chromVAR_celltype_stacked_barplots_neural_crest.pdf",io$outdir), width=6, height=6)
print(p)
# dev.off()

to.plot2 <- to.plot %>%
  .[,.(value=mean(value)),by="celltype"] %>% 
  .[,value:=round(minmax.normalisation(value),2)]

# ggpie(to.plot2, x="value", label = "value", fill="celltype", color="black") +
#   scale_fill_manual(values=opts$celltype.colors[names(opts$celltype.colors)%in%unique(to.plot2$celltype)]) +
#   labs(x="", y="") +
#   theme(
#     legend.position = "none"
#   )

ggplot(to.plot2, aes(x="", y=value, fill=celltype)) +
  geom_bar(stat="identity", width=1, color="black") +
  scale_fill_manual(values=opts$celltype.colors) +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(
    legend.position="none"
  )

## CREATE A NULL DISTRIBUTION AND CALCUALTE P-VALUES FOR EACH CELL TYPE ##


tf.to.plot <- "TAL1"
genes.to.plot <- GRN_coef.dt[tf==tf.to.plot & pvalue<=0.01 & beta>0,gene]

foo <- 1:100 %>% map(function(i) {
  set.seed(i)
  rna_targets_random.mtx <- matrix(sample(rna_targets.mtx), nrow=nrow(rna_targets.mtx), ncol=ncol(rna_targets.mtx))
  dimnames(rna_targets_random.mtx) <- dimnames(rna_targets.mtx)
  
  apply(rna_targets_random.mtx[genes.to.plot,], 1, function(i) i/sum(i)) %>%
    as.data.table(keep.rownames = T) %>% setnames("rn","celltype") %>%
    melt(id.vars="celltype", variable.name="gene") %>%
    .[,.(value=mean(value)),by="celltype"] %>%
    .[,iteration:=1]
}) %>% rbindlist
  

to.plot <- opts$celltypes %>% map(function(i) {
  bar <- apply(rna_targets.mtx[genes.to.plot,], 1, function(i) i/sum(i)) %>%
    as.data.table(keep.rownames = T) %>% setnames("rn","celltype") %>%
    melt(id.vars="celltype", variable.name="gene") %>%
    .[,.(value=mean(value)),by="celltype"]
  
  data.table(
    celltype = i,
    log_pvalue = mean(bar[celltype==i,value]<=foo[celltype==i,value]) %>% -log10(.+1e-64)
  )
  
}) %>% rbindlist


ggplot(to.plot, aes(x="", y=log_pvalue, fill=celltype)) +
  geom_bar(stat="identity", width=1, color="black") +
  scale_fill_manual(values=opts$celltype.colors) +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(
    legend.position="none"
  )
