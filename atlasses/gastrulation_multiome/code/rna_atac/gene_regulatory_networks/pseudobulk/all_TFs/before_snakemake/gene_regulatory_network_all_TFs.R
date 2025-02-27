suppressMessages(library(argparse))
suppressMessages(library(GGally))
suppressMessages(library(igraph))
suppressMessages(library(network))
suppressMessages(library(sna))
suppressMessages(library(intergraph))

#####################
## Define settings ##
#####################

# load default setings
source(here::here("settings.R"))
source(here::here("utils.R"))

# I/O
io$outdir <- file.path(io$basedir,"results/rna_atac/gene_regulatory_networks/pseudobulk/all_TF")

# Options
opts$TFs <- list.files(io$virtual_chip.dir, pattern = "*.bed.gz") %>% stringr::str_replace_all(".bed.gz","")

##############################
## Load pseudobulk RNA data ##
##############################

# Load SingleCellExperiment
sce.pseudobulk <- readRDS(io$rna.pseudobulk.sce)

# Define TFs
opts$TFs <- opts$TFs[str_to_title(opts$TFs)%in%rownames(sce.pseudobulk)]

#############################
## Load tf2tf correlations ##
#############################

tf2tf_cor.se <- readRDS(io$tf_cor.se)[opts$TFs,opts$TFs]

#############################
## Load peak2gene linkages ##
#############################

peak2gene.dt <- fread(io$archR.peak2gene.all) %>% 
  # .[gene%in%opts$genes] %>%
  .[,peak:=sprintf("chr%s:%s-%s",chr,peak.start,peak.end)]

################################
## Load virtual ChIP-seq data ##
################################

virtual_chip.mtx <- readRDS(io$virtual_chip.mtx)[,opts$TFs]  

######################################################
## Use ChIP-seq data to filter TF-TF associations ##
######################################################

opts$cor.threshold <- 0.30

corr.mtx <- as.matrix(assay(tf2tf_cor.se,"cor"))
for (i in opts$TFs) {

  # Select target peaks
  # target_peaks_i <- virtual_chip.mtx[,i]
  
  # Select target genes
  target_genes_i <- intersect(
    x = names(which(abs(corr.mtx[i,])>=opts$cor.threshold)),  # all correlations
    # x = names(which(corr.mtx$r[i,]>=opts$cor.threshold)),   # only positive correlations
    y = unique(peak2gene.dt[peak%in%target_peaks_i,gene])
  )
  
  print(sprintf("%s: %s target genes",i,length(target_genes_i)))
  
  # Update correlation matrix
  corr.mtx[i,!colnames(corr.mtx)%in%target_genes_i] <- NA
}

# Filter TFs and genes with too little connections
corr.mtx <- corr.mtx[rowSums(!is.na(corr.mtx))>=1,,drop=F]
corr.mtx <- corr.mtx[,colSums(!is.na(corr.mtx))>=1,drop=F]

opts$TFs <- rownames(corr.mtx)

####################
## Create network ##
####################

# Prepare data
node_list.dt <- data.table(node_id=1:nrow(corr.mtx), node_name=rownames(corr.mtx))
target_list.dt <- data.table(target_id=1:ncol(corr.mtx), target_name=colnames(corr.mtx))

edge_list.dt <- as.data.table(corr.mtx,keep.rownames = T) %>% 
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

###########################
## Create network layout ##
###########################

i <- args$celltype

set.seed(42)
# layout <- layout.fruchterman.reingold(net)
layout <- layout.sphere(net)

#################################
## Fetch RNA expression values ##
#################################

# Calculate expression values for each cell type
rna.mtx.TFs <- logcounts(sce.pseudobulk)[str_to_title(rownames(corr.mtx)),,drop=F]; rownames(rna.mtx.TFs) <- toupper(rownames(rna.mtx.TFs))
rna.mtx.target_genes <- logcounts(sce.pseudobulk.filt)[rownames(sce.pseudobulk.filt)%in%opts$genes,,drop=F]

rna.mtx.TFs.scaled <- apply(rna.mtx.TFs,1,minmax.normalisation) %>% t
rna.mtx.target_genes.scaled <- apply(rna.mtx.target_genes,1,minmax.normalisation) %>% t
# rna.mtx.target_genes.scaled <- scale(rna.mtx.target_genes, center=TRUE, scale = TRUE)

# Add gene expression values as network attributes
net <- set_vertex_attr(net, name = "expr", value = c(rna.mtx.TFs[,i],rna.mtx.target_genes[,i])[names(V(net))])

#####################
## Plot PAGA graph ##
#####################

source(here::here("load_paga_graph.R"))

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

########################
## Define node format ##
########################

pal <- grDevices::colorRamp(c("gray60", "purple"))( (1:100)/100 )
pal <- cbind(pal, seq(100, 255, length.out = 100))
pal.TFs <- grDevices::colorRamp(c("gray60", "purple"))( (1:100)/100 )
pal.TFs <- cbind(pal.TFs, seq(100, 255, length.out = 100))
pal.genes <- grDevices::colorRamp(c("gray60", "darkgreen"))( (1:100)/100 )
pal.genes <- cbind(pal.genes, seq(100, 255, length.out = 100))

foo <- rna.mtx.TFs.scaled[,i]; names(foo) <- rownames(rna.mtx.TFs)
bar <- rna.mtx.target_genes.scaled[,i]; names(bar) <- rownames(rna.mtx.target_genes.scaled)
rna.expr.i <- c(foo, bar)[names(V(net))]# %>% minmax.normalisation

# Discretise a bit
# rna.expr.i[rna.expr.i<0.20] <- 0
# rna.expr.i[rna.expr.i>0.80] <- 1

# Define node colours
V(net)$color <- c(
  colourvalues::colour_values(rna.expr.i[names(V(net))[V(net)$class=="TF"]], palette = pal.TFs), 
  colourvalues::colour_values(rna.expr.i[names(V(net))[V(net)$class=="gene"]], palette = pal.genes)
)
# V(net)$color <- colourvalues::colour_values(rna.expr.i, palette = pal)

# Define node shapes
V(net)$shape <- stringr::str_replace_all(V(net)$class,c("gene"="circle","TF"="triangle"))

# Define node label size
node.label.size <- c(3,5)[factor(V(net)$class)]
node.label.size[rna.expr.i<0.05] <- 0

########################
## Define edge format ##
########################

# Define edge colors
E(net)$color <- "gray70"  # Colour negative correlation edges as blue
# E(net)[which(E(net)$weight<0)]$color <- "blue"  # Colour negative correlation edges as blue
# E(net)[which(E(net)$weight>0)]$color <- "red"   # Colour positive correlation edges as red

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

# igraph::degree(net)
# igraph::closeness(net)
# igraph::eigen_centrality(net)$vector
# igraph::betweenness(net)

############################################
## Plot number of correlated TFs per gene ##
############################################

genes.to.plot <- tmp[,.(N=sum(abs(cor)>=0.15)),by="gene"] %>% setorder(-N) %>% head(n=15) %>% .$gene
to.plot <- tmp %>%
  .[,.(N=sum(abs(cor)>=0.15)),by=c("gene","sign")] %>%
  .[gene%in%genes.to.plot] %>% .[,gene:=factor(gene,levels=genes.to.plot)]

p4 <- ggbarplot(to.plot, x="gene", y="N", fill="sign", width=0.55) +
  coord_flip() +
  scale_fill_manual(values=c("-"="blue", "+"="red"), drop=F) +
  labs(x="", y="Number of correlated TFs") +
  theme(
    axis.text.y = element_text(size=rel(0.75)),
    axis.text.x = element_text(colour="black",size=rel(0.75)),
    axis.ticks.x = element_line(size=rel(0.75)),
    legend.position = "top"
  )

###################
## Combine plots ##
###################

p <- cowplot::plot_grid(plotlist=list(p1,p2,p3,p4), rel_widths = c(2/5,3/5), rel_heights = c(2/3,1/3), nrow=2)

outfile <- sprintf("%s/%s_network.png",args$outdir,i)
png(outfile, width = 1200, height = 850, bg = "white")
print(p)
dev.off()


##################
## Save network ##
##################

outfile <- sprintf("%s/%s_network.rds",args$outdir,i)
saveRDS(net, outfile)
