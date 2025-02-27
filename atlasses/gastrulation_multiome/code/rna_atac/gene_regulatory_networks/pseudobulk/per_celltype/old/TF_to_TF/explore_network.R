library(igraph) 
library(visNetwork) 
library(colourvalues) 

#####################
## Define settings ##
#####################

# load default setings
source(here::here("settings.R"))
source(here::here("utils.R"))

file <- "/Users/argelagr/data/gastrulation_multiome_10x//results/rna_atac/gene_regulatory_networks/pseudobulk/per_celltype/Gut/Gut_network.rds"

##################
## Load network ##
##################

net <- readRDS(file)

# Modify attributes
V(net)$group <- V(net)$class
# V(ig)$label <- paste0("Node",1:length(V(ig)))

V(net)$shape <- stringr::str_replace_all(V(net)$class,c("gene"="triangle","TF"="circle"))


################
## vizNetwork ##
################

# size=c(15,20)[factor(V(net)$class)]
# Plot
visIgraph(net, randomSeed=42) %>%
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE, selectedBy = "group") %>%
  # visNodes(shadow = FALSE, font='20px arial black', size=50)
  visGroups(groupname = "TF", shape = "circle", size=35, font = list(color="black", size=30), color = list(background="#63B8FF", hover="#4876FF", border="#63B8FF")) %>%
  visGroups(groupname = "gene", shape = "triangle", size=20, font = list(color="black", size=35), color = list(background="#EE6363", hover="red", border="#EE6363"))

##################
## Perturbation ##
##################

# x0 <- rep(0,length(V(net)))
# x0 <- rna.expr.i
# x1 <- x0 %*% A

col.palette <- grDevices::colorRamp(c("gray80", "red"))( (1:100)/100 )
col.palette <- cbind(col.palette, seq(100, 255, length.out = 100))

A <- as_adjacency_matrix(net, sparse=F, attr="weight")

tf.ko <- "GATA3"  

diff <- A[tf.ko,]

V(net)$color <- colour_values(diff, palette = col.palette)
V(net)$color[names(V(net)) == tf.ko] <- "red"
E(net)$color <- colour_values(E(net)$weight, palette = col.palette)
E(net)$color[names(tail_of(net, es=E(net)))!=tf.ko] <- "gray80"

visIgraph(net, randomSeed=42)

#########################
## Centrality measures ##
#########################

##########
## TEST ##
##########


# ggnet2(
#   mode = layout, # "fruchtermanreingold",
#   net = net,
#   color = V(net)$color,
#   edge.color = E(net)$color,
#   edge.size = edge_weight,
#   node.size = c(15,20)[factor(V(net)$class)],
#   # node.alpha = alphas,
#   node.shape = V(net)$shape,
#   label = TRUE,
#   label.color = "black",
#   label.size = node.label.size,
#   arrow.size = 0.1,
#   legend.position = "none"
# ) 


# visIgraph(ig, idToLabel = FALSE) %>% 
#   visIgraphLayout(randomSeed=42, physics = TRUE) %>%
#   visOptions(
#     highlightNearest = TRUE,
#     nodesIdSelection = TRUE,
#     selectedBy = "group"
#   ) %>%
#   visEdges(
#     width = 3,
#     color = "black",
#     smooth = FALSE # better for performance
#   ) %>%
#   visNodes(
#     # color = "black",
#     shadow = TRUE
#   ) %>%
#   visGroups(groupname = "A", shape = "circle", size=35, font = list(color="black", size=30), color = list(background="#63B8FF", hover="#4876FF", border="#63B8FF")) %>%
#   visGroups(groupname = "B", shape = "triangle", size=20, font = list(color="black", size=35), color = list(background="#EE6363", hover="red", border="#EE6363")) %>%
#   visPhysics(
#     solver = "forceAtlas2Based", 
#     forceAtlas2Based = list(gravitationalConstant = -500), 
#     stabilization = TRUE # By default, vis.js computes coordinates dynamically and waits for stabilization before rendering
#     
#   ) %>%
#   visInteraction(
#     hover = TRUE,
#     dragNodes = TRUE, 
#     dragView = TRUE, 
#     zoomView = TRUE
#   )