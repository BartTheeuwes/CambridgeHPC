library(igraph)
# library(GGally)
# library(network)
# library(sna)

#####################
## Define settings ##
#####################

# Load default settings
if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
  source("/Users/ricard/gastrulation_multiome_10x/utils.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
  source("/homes/ricard/gastrulation_multiome_10x/utils.R")
} else {
  stop("Computer not recognised")
}

# I/O
io$outdir <- paste0(io$basedir,"/results/test/networks")

##############################
## Load correlation  matrix ##
##############################

corr.mtx <- readRDS("/Users/ricard/data/gastrulation_multiome_10x/results/rna/coexpression/correlation_matrix_tf2gene.rds")

marker_genes.dt <- fread(io$rna.atlas.marker_genes) %>% .[celltype%in%opts$celltypes]
corr.mtx <- corr.mtx[,colnames(corr.mtx) %in% unique(marker_genes.dt$gene)]

# corr.mtx <- corr.mtx[1:100,1:1000]

diag(corr.mtx) <- 0

corr.mtx[is.na(corr.mtx)] <- 0

corr.mtx[abs(corr.mtx)<0.50] <- 0

# Remove unconnected nodes
corr.mtx <- corr.mtx[rowSums(corr.mtx>0)>50,]
corr.mtx <- corr.mtx[,colSums(corr.mtx>0)>3]

# sum(colSums(corr.mtx>0)<3)
# sum(rowSums(corr.mtx>0)<50)

dim(corr.mtx)
mean(corr.mtx==0)

##################
## Prepare data ##
##################

node_list.dt <- data.table(node_id=1:nrow(corr.mtx), node_name=rownames(corr.mtx))
target_list.dt <- data.table(target_id=1:ncol(corr.mtx), target_name=colnames(corr.mtx))

edge_list.dt <- as.data.table(corr.mtx,keep.rownames = T) %>% 
  setnames("rn","from") %>%
  melt(id.vars=c("from"), variable.name="to", value.name="weight") %>%
  .[weight>0]


node_list_metadata.dt <- data.table(
  label = c(node_list.dt$node_name, target_list.dt$target_name),
  class = c(rep("TF",nrow(node_list.dt)), rep("gene",nrow(target_list.dt)))
)


####################
## Create network ##
####################

# edge_list.df <- edge_list.dt[,c("node_id","target_id","value")] %>%
#   setnames(c("from","to","weight")) %>%
#   tibble::as_tibble()

# node_list.df <- node_list.dt %>% setnames(c("id","label")) %>% tibble::as_tibble()
# network <- network(edge_list.df, vertex.attr = node_list.dt, matrix.type = "edgelist", ignore.eval = FALSE)

###################
## Create igraph ##
###################

# If vertices is NULL: the first two columns of d are used as a symbolic edge list and additional columns as edge attributes. The names of the attributes are taken from the names of the columns.
# net <- graph_from_data_frame(d = edge_list.dt)

# If vertices is not NULL, then it must be a data frame giving vertex metadata. The first column of vertices is assumed to contain symbolic vertex names, this will be added to the graphs as the ‘name’ vertex attribute. Other columns will be added as additional vertex attributes.
net <- graph_from_data_frame(d = edge_list.dt, vertices = node_list_metadata.dt)

###################################
## Basic plotting of the network ##
###################################

# plot(net, layout = layout_with_graphopt, edge.arrow.size = 0.2, vertex.cex = 3)

png(sprintf("%s/network_test.png",io$outdir), width=900, height=800)
plot(net, layout = layout_with_graphopt)
dev.off()

###################
## Change layout ##
###################

layout_in_circles <- function(net, group=1) {
  layout <- lapply(split(V(net), group), function(x) {
    layout_in_circle(induced_subgraph(net,x))
  })
  layout <- Map(`*`, layout, seq_along(layout))
  x <- matrix(0, nrow=vcount(net), ncol=2)
  split(x, group) <- layout
  x
}


# plot(net, layout = layout_in_circles(net, group=V(net)$group))

##############################################
## Plot network, change display attributes ##
##############################################

# net <- set_vertex_attr(net, name="color", index=rownames(corr.mtx), value="blue")
# net <- set_vertex_attr(net, name="color", index=colnames(corr.mtx), value="red")

# Modify nodes
V(net)$group <- factor(V(net)$class, levels=c("TF","gene"))
V(net)$label <- names(V(net))# NA
V(net)$size <- c(3.5,5.5)[factor(V(net)$class)]
V(net)$label.color <- "black"
V(net)$color <- c("tomato", "gold")[factor(V(net)$class)]
V(net)$label.cex <- 0.5

# Modify edges
E(net)$arrow.size <- 0.1
E(net)$edge.color <- "gray80"
# E(net)$edge.width <- 0.25
E(net)$width <- E(net)$edge.width <- E(net)$weight#/6

# plot(net, layout = layout.fruchterman.reingold, edge.arrow.size = 0.1, vertex.cex = 3)

# pdf(sprintf("%s/network_test.pdf",io$outdir), width=10, height=9)
png(sprintf("%s/network_test.png",io$outdir), width=900, height=800)
# plot(net, layout = layout_in_circles(net, group=V(net)$group))
plot(net, layout = layout.fruchterman.reingold)
dev.off()

#################################
## Overlay with RNA expression ##
#################################

# Load RNA expr
sce.pseudobulk <- readRDS(io$rna.pseudobulk.sce)
rownames(sce.pseudobulk)[toupper(rownames(sce.pseudobulk))%in%rownames(corr.mtx)] <- toupper(rownames(sce.pseudobulk)[toupper(rownames(sce.pseudobulk))%in%rownames(corr.mtx)])
rna.mtx <- logcounts(sce.pseudobulk)[unlist(dimnames(corr.mtx)),]

# Define palette
pal <- grDevices::colorRamp(c("gray60", "purple"))( (1:100)/100 )
pal <- cbind(pal, seq(100, 255, length.out = 100))

celltypes.to.plot <- colnames(rna.mtx) %>% head(n=3)

# i <- "Allantois"
for (i in celltypes.to.plot) {
   rna.expr.i <- rna.mtx[,i][names(V(net))]
   V(net)$color <- colourvalues::colour_values(rna.expr.i, palette = pal)
   
   # pdf(sprintf("%s/%s_network_test.pdf",io$outdir,i), width=6, height=5)
   png(sprintf("%s/%s_network_test.png",io$outdir,i), width=900, height=800)
   plot(net, layout = layout_in_circles(net, group=V(net)$group), vertex.label=NA, main=i)
   dev.off()
}


##############################################
## Plot network, colour by the celltype with highest expression ##
##############################################

colors <- opts$celltype.colors[apply(rna.mtx, 1,function(x) colnames(rna.mtx)[which.max(x)])]

V(net)$size <- c(1,4)[factor(V(net)$class)]

# png(sprintf("%s/network_test.png",io$outdir), width=1400, height=1200)
pdf(sprintf("%s/network_test.pdf",io$outdir), width=18, height=14)
plot(net, layout = layout.fruchterman.reingold, vertex.color=colors)
# plot(net, layout = layout.mds, vertex.color=colors)
# plot(net, layout = layout.drl, vertex.color=colors)
dev.off()

################
## Create GIF ##
################

library(magick)

## list file names and read in
imgs.files <- c(
  sprintf("%s/network_test.png",io$outdir,i),
  sprintf("%s/%s_network_test.png",io$outdir,celltypes.to.plot)
  )
img_list <- lapply(imgs.files, image_read)

# join the images together
img_joined <- image_join(img_list)

# animate
img_animated <- image_animate(img_joined, fps = 1)

# save
image_write(image = img_animated, path = sprintf("%s/test.gif",io$outdir))


#########################
## Export to cytoscape ##
#########################

library(RCy3)

cytoscapePing()

createNetworkFromIgraph(net, paste0(io$outdir,"cytoscape_test"))

