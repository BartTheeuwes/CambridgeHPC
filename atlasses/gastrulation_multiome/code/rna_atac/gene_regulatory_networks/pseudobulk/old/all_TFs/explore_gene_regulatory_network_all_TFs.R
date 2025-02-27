suppressMessages(library(GGally))
suppressMessages(library(igraph))
suppressMessages(library(network))
suppressMessages(library(sna))
suppressMessages(library(intergraph))
require(visNetwork)


#####################
## Define settings ##
#####################

# load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

## START TEST
io$correlation_matrix <- paste0(io$basedir,"/results_new/rna_atac/gene_regulatory_networks/pseudobulk/all_TFs/CISBP_connectivity_matrix_all_TFs.rds")
io$network <- paste0(io$basedir,"/results_new/rna_atac/gene_regulatory_networks/pseudobulk/all_TFs/CISBP_network_all_TFs.rds")
io$outdir <- paste0(io$basedir,"/results_new/rna_atac/gene_regulatory_networks/pseudobulk/all_TFs")

## END TEST

# I/O
dir.create(args$outdir, showWarnings = F)

# Options
# opts$TFs <- list.files(io$virtual_chip.dir, pattern = "*.bed.gz") %>% stringr::str_replace_all(".bed.gz",""))

##############################
## Load pseudobulk RNA data ##
##############################

sce.pseudobulk <- readRDS(io$rna.pseudobulk.sce)

##############
## Load GRN ##
##############

connectivity.mtx <- readRDS(io$correlation_matrix)
net <- readRDS(io$network)

###########################
## Create network layout ##
###########################

set.seed(42)
# layout <- layout.fruchterman.reingold(net)
layout <- layout.sphere(net)

#################################
## Fetch RNA expression values ##
#################################

TFs <- unique(names(V(net)))

# Calculate expression values for each cell type
rna.mtx <- logcounts(sce.pseudobulk)[str_to_title(TFs),,drop=F]; rownames(rna.mtx) <- toupper(rownames(rna.mtx))

rna.mtx.scaled <- apply(rna.mtx,2,minmax.normalisation)# %>% t
# rna.mtx.target_genes.scaled <- apply(rna.mtx.target_genes,1,minmax.normalisation) %>% t
# rna.mtx.target_genes.scaled <- scale(rna.mtx.target_genes, center=TRUE, scale = TRUE)

# Add gene expression values as network attributes
for (i in colnames(sce.pseudobulk)) { 
  net <- set_vertex_attr(net, name = sprintf("expr_%s",i), value = rna.mtx.scaled[,i][names(V(net))])
}
vertex_attr(net,"expr_Epiblast")



########################
## Define edge format ##
########################

# Define edge colors
# E(net)$color <- "gray70"
# E(net)[which(E(net)$weight<0)]$color <- "blue"
# E(net)[which(E(net)$weight>0)]$color <- "red"

# Define edge width
# edge_weight <- abs(E(net)$weight)   # Convert edge weights to absolute values
# for (j in seq_len(length(E(net)))) {
#   edge_weight[[j]] <- edge_weight[[j]]*(rna.expr.i[tail_of(net, j)] * rna.expr.i[head_of(net, j)])
# }
# edge_weight[edge_weight==0] <- 0.01

# alphas <- rna.expr.i

# Discretise a bit
# alphas[rna.expr.i<0.5] <- 0.25
# alphas[rna.expr.i>0.75] <- 1

##################
## Plot network ##
##################

p2 <- ggnet2(
  mode = layout, # "fruchtermanreingold",
  net = net,
  color = V(net)$color,
  edge.color = E(net)$color,
  edge.size = edge_weight,
  node.size = c(15,20)[factor(V(net)$class)],
  # node.alpha = alphas,
  node.shape = V(net)$shape,
  label = TRUE,
  label.color = "black",
  label.size = node.label.size,
  arrow.size = 0.1,
  legend.position = "none"
) 

# test
ggnet2(
  mode = layout, # "fruchtermanreingold",
  net = net,
  color = V(net)$color,
  # edge.color = E(net)$color,
  edge.size = 0.1,
  node.size = 2.5,
  label = TRUE,
  label.color = "black",
  label.size = 1.5,
  arrow.size = 0.1,
  legend.position = "none"
) 

#########################
## Centrality measures ##
#########################

igraph::degree(net) %>% sort %>% tail
igraph::closeness(net) %>% sort %>% tail
igraph::eigen_centrality(net)$vector %>% sort %>% tail
igraph::betweenness(net) %>% sort %>% tail

# p3 <- ggbarplot(to.plot, x="TF", y="N", fill="sign", width=0.55) +
#   coord_flip() +
#   scale_fill_manual(values=c("-"="blue", "+"="red"), drop=F) +
#   labs(x="", y="Number of correlated genes") +
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
# 
# outfile <- sprintf("%s/%s_network.png",args$outdir,i)
# png(outfile, width = 1200, height = 850, bg = "white")
# print(p)
# dev.off()



#####################
## Plot PAGA graph ##
#####################

source(here::here("load_paga_graph.R"))

colors <- rep("gray70",length(opts$celltypes)); names(colors) <- opts$celltypes
alphas <- rep(0.5,length(opts$celltypes)); names(alphas) <- opts$celltypes
sizes <- rep(6,length(opts$celltypes)); names(sizes) <- opts$celltypes

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



################
## visNetwork ##
################

nodes <- names(V(net))

pal <- grDevices::colorRamp(c("gray62", "purple"))( (1:100)/100 )
pal <- cbind(pal, seq(100, 255, length.out = 100))



# Edge width and color
edge_weight <- abs(E(net)$weight)   # Convert edge weights to absolute values
for (j in seq_len(length(E(net)))) {
  edge_weight[[j]] <- edge_weight[[j]]*(rna.mtx.scaled[tail_of(net, j)] * rna.mtx.scaled[head_of(net, j)])
}
edge_weight[edge_weight==0] <- 0.01
E(net)$weight <- edge_weight*10
# names(edge_attr(net))

color_by <- "tmp"
celltype.to.plot <- "NMP"

# Node color
if (color_by=="eigenvalue_centrality") {
  eigenvalue_centrality_scores <- igraph::eigen_centrality(net)$vector[nodes]
  color_nodes <- colourvalues::colour_values(eigenvalue_centrality_scores, palette = pal)
  names(color_nodes) <- nodes
} else if (color_by=="degree_centrality") {
  degree_centrality_scores <- igraph::degree(net)[nodes]
  color_nodes <- colourvalues::colour_values(degree_centrality_scores, palette = pal)
  names(color_nodes) <- nodes
} else if (color_by=="expression") {
  stopifnot(sprintf("expr_%s",celltype.to.plot)%in%names(vertex_attr(net)))
  expression_values <- vertex_attr(net)[[sprintf("expr_%s",celltype.to.plot)]]; names(expression_values) <- nodes
  color_nodes <- colourvalues::colour_values(minmax.normalisation(expression_values), palette = pal)
  names(color_nodes) <- names(expression_values)
} else if (color_by=="tmp") {
  marker_TFs <- fread(io$rna.atlas.marker_TFs.up) %>% .[,gene] %>% unique %>% toupper
  marker_TFs <- marker_TFs[marker_TFs%in%rownames(rna.mtx)]
  color_nodes <- rep("green",length(nodes)); names(color_nodes) <- nodes
  color_nodes[marker_TFs] <- opts$celltype.colors[colnames(rna.mtx)[apply(rna.mtx[marker_TFs,],1,which.max)]]
}

V(net)$color <- color_nodes[nodes]


visIgraph(net) %>%
  visIgraphLayout(randomSeed=42, physics = FALSE) %>%
  # visNodes(shadow = TRUE, size=15) %>%
  visEdges(color = list(highlight = "black", hover = "black"), selectionWidth=10, arrows=list("to"=list("enabled"=F))) %>%
  # visEdges(width=E(net)$width, selectionWidth=2, color = "gray60", smooth = FALSE, arrows=list("to"=list("enabled"=F))) %>%
  visOptions(highlightNearest = list(enabled=TRUE), nodesIdSelection = list(enabled=TRUE)) %>%
  # visOptions(highlightNearest = list(enabled=TRUE), nodesIdSelection = list(enabled=TRUE, style = 'color: black, width: 500px; height: 50px;')) %>%
  # visGroups(groupname = "TF", shape = "circle", size=35, font = list(color="black", size=30), color = list(background="#63B8FF", hover="#4876FF", border="#63B8FF")) %>%
  # visGroups(groupname = "gene", shape = "triangle", size=20, font = list(color="black", size=35), color = list(background="#EE6363", hover="red", border="#EE6363")) %>%
  # visPhysics(
  #   solver = "forceAtlas2Based",
  #   minVelocity = 0.50,
  #   forceAtlas2Based = list(gravitationalConstant = -250),
  #   stabilization = TRUE # By default, vis.js computes coordinates dynamically and waits for stabilization before rendering
  # ) %>%
  visInteraction(
    hover = TRUE,
    dragNodes = TRUE,
    dragView = TRUE,
    zoomView = TRUE
  )
