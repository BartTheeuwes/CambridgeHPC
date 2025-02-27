suppressMessages(library(argparse))
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
p$add_argument('--outdir',       type="character",                help='Output file')
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

# load default setings
source(here::here("settings.R"))
source(here::here("utils.R"))

## START TEST
args$distance <- 1e5
args$trajectory <- paste0(io$basedir,"/results/rna/trajectories/blood_trajectory/blood_trajectory.txt.gz")
args$trajectory_name <- "blood"
args$outdir <- paste0(io$basedir,"/results/rna_atac/gene_regulatory_networks/trajectories/blood")
## END TEST

#####################
## Load trajectory ##
#####################

trajectory.dt <- fread(args$trajectory)

##################
## Load network ##
##################

# saveRDS(tmp, sprintf("%s/corr_matrix.rds",args$outdir))

###################
## Load metadata ##
###################

sample_metadata <- fread(io$metadata) %>%
  .[pass_atacQC==TRUE & pass_rnaQC==TRUE & doublet_call==FALSE]

cells <- intersect(trajectory.dt$cell,sample_metadata$cell)
trajectory.dt <- trajectory.dt[cell%in%cells] 
sample_metadata <- sample_metadata[cell%in%cells] %>% setkey(cell) %>% .[cells]

opts$celltypes.subset <- opts$celltypes[opts$celltypes%in%unique(sample_metadata$celltype.predicted)]

#######################################
## Load virtual ChIP-seq annotation ###
#######################################

virtual_chip.mtx <- readRDS(io$virtual_chip.mtx)

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

rm(rna.sce)


####################
## Create network ##
####################

#################################
## Plot number of edges per TF ##
#################################

to.plot <- edge_list.dt[,.N,by="from"] %>% setnames("from","TF")

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
# print(p)
# dev.off()

###########################################################################
## Plot network, colouring each gene by celltype-specific RNA expression ##
###########################################################################

# Load PAGA graph
source(here::here("load_paga_graph.R"))

# Define colors
pal <- grDevices::colorRamp(c("gray60", "purple"))( (1:100)/100 )
pal <- cbind(pal, seq(100, 255, length.out = 100))

rna_scaled.mtx <- apply(logcounts(rna_tf_pseudobulk.sce),1,minmax.normalisation) %>% t

celltypes.to.plot <- opts$celltypes.subset

for (i in celltypes.to.plot) {
  
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
    edge_weight[[j]] <- edge_weight[[j]]*(rna.expr.i[tail_of(net, j)] * rna.expr.i[head_of(net, j)])
  }
  edge_weight[edge_weight==0] <- 0.03
  # edge_weight <- minmax.normalisation(edge_weight)
  edge_weight[edge_weight>0.75] <- 0.75
  
  alphas <- rna.expr.i
  # alphas[rna.expr.i<0.5] <- 0.25
  # alphas[rna.expr.i>0.75] <- 1
  
  p2 <- ggnet2(
    mode = layout,
    net = net,
    color = V(net)$color,
    edge.color = E(net)$color,
    edge.size = edge_weight,
    node.size = 7,
    node.alpha = alphas,
    label = TRUE,
    label.color = "black",
    label.size = 3.5,
    arrow.size = 0.1,
    legend.position = "none"
  ) 
  
  ###################
  ## Combine plots ##
  ###################
  
  p <- cowplot::plot_grid(plotlist=list(p1,p2), rel_widths = c(1.5/5,3.5/5), nrow=1)
    
  outfile <- sprintf("%s/%d_%s.png",args$outdir,match(i,opts$celltypes.subset),i)
  png(outfile, width = 1200, height = 550, bg = "white")
  print(p)
  dev.off()
}


