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

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--celltype',    type="character",                help='celltype to obtain marker genes from')
p$add_argument('--rna_sce',          type="character",   help='RNA SingleCellExperiment')
p$add_argument('--tf2tf_cor_mtx',          type="character",   help='')
p$add_argument('--virtual_chip',          type="character",   help='')
p$add_argument('--outdir',       type="character",               help='Output file')
args <- p$parse_args(commandArgs(TRUE))

## START TEST
args$celltype <- c("Gut")
args$outdir <- file.path(io$basedir,"/results_new/rna_atac/gene_regulatory_networks/pseudobulk/per_celltype/Gut")
## END TEST

dir.create(args$outdir, showWarnings = F)

##############################
## Load marker gene and TFs ##
##############################

# TO-DO: INCLUDE DOWNREGULATION EVENTS

io$marker_tfs <- "/Users/argelagr/data/gastrulation_multiome_10x/results_new/rna/differential/TFs/marker_TFs/marker_TFs_upregated_filtered.txt.gz"

marker_TFs.dt <- fread(io$marker_tfs) %>% 
  .[celltype%in%args$celltype & score>=0.70]

##############################
## Load pseudobulk RNA data ##
##############################

# Load SingleCellExperiment
rna_pseudobulk.sce <- readRDS(io$rna.pseudobulk.sce)#[,opts$celltypes]

# Define TFs
opts$TFs <- marker_TFs.dt$gene %>% toupper()
# opts$TFs <- opts$TFs[str_to_title(opts$TFs)%in%rownames(rna_pseudobulk.sce)]
opts$TFs <- opts$TFs[opts$TFs%in%(list.files(io$virtual_chip.dir, pattern = "*.bed.gz") %>% stringr::str_replace_all(".bed.gz",""))]

###############################
## Load TF2TF correlations ##
###############################

# Load
tf2tf_cor.mtx <- readRDS(io$tf2tf_cor.mtx)

# Filter TFs
opts$TFs <- opts$TFs[opts$TFs%in%unique(c(rownames(tf2tf_cor.mtx),colnames(tf2tf_cor.mtx)))]
tf2tf_cor.mtx <- tf2tf_cor.mtx[opts$TFs,opts$TFs]

#############################
## Load peak2gene linkages ##
#############################

peak2gene.dt <- fread(io$archR.peak2gene.all) %>% 
  .[,gene:=toupper(gene)] %>%
  .[gene%in%opts$TFs] %>%
  .[,peak:=sprintf("chr%s:%s-%s",chr,peak.start,peak.end)]

################################
## Load virtual ChIP-seq data ##
################################

# Load
virtual_chip.mtx <- readRDS(io$virtual_chip.mtx)

# Filter peaks and TFs
virtual_chip.mtx <- virtual_chip.mtx[rownames(virtual_chip.mtx)%in%unique(peak2gene.dt$peak),]
virtual_chip.mtx <- virtual_chip.mtx[,opts$TFs]

######################################################
## Use virtual ChIP-seq data to filter TF2TF links ##
######################################################

opts$min.cor.threshold <- 0.30
opts$min.chip.score <- 0.15

corr.mtx <- tf2tf_cor.mtx
# i <- "PAX9"
for (i in opts$TFs) {

  # Select target peaks
  target_peaks_i <- names(which(abs(virtual_chip.mtx[,i])>=opts$min.chip.score))
  
  # Select target genes
  target_genes_i <- intersect(
    x = names(which(abs(corr.mtx[i,])>=opts$min.cor.threshold)),  # all correlations
    # x = names(which(corr.mtx$r[i,]>=opts$min.cor.threshold)),   # only positive correlations
    y = unique(peak2gene.dt[peak%in%target_peaks_i,gene])
  )
  
  print(sprintf("%s: %s target peaks and %s target genes",i,length(target_peaks_i),length(target_genes_i)))
  
  # Update correlation matrix
  corr.mtx[i,!colnames(corr.mtx)%in%target_genes_i] <- NA
}

# Filter TFs and genes with too little connections
# corr.mtx <- corr.mtx[rowSums(!is.na(corr.mtx))>=1,,drop=F]
# corr.mtx <- corr.mtx[,colSums(!is.na(corr.mtx))>=1,drop=F]
# opts$TFs <- rownames(corr.mtx)

diag(corr.mtx) <- NA

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

net <- graph_from_data_frame(d = edge_list.dt, directed = T)

###########################
## Create network layout ##
###########################

i <- args$celltype

set.seed(42)
# layout <- layout.fruchterman.reingold(net)
layout <- layout.sphere(net)

###################################################
## Add RNA expression values as a node attribute ##
###################################################

# Fetch RNA expression values
tf_pseudobulk.sce <- rna_pseudobulk.sce[str_to_title(names(V(net))),]
rownames(tf_pseudobulk.sce) <- toupper(rownames(tf_pseudobulk.sce)) 

# Calculate expression values for each cell type
rna_tf.mtx <- logcounts(tf_pseudobulk.sce)[names(V(net)),,drop=F]

rna_tf_scaled.mtx <- apply(rna_tf.mtx,1,minmax.normalisation) %>% t

# Add gene expression values as network attributes
net <- set_vertex_attr(net, name = "expr", value = rna_tf.mtx[,i])

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

rna.expr.i <- rna_tf_scaled.mtx[,i]; names(rna.expr.i) <- rownames(rna_tf.mtx)

# Define node colours
V(net)$color <- colourvalues::colour_values(rna.expr.i, palette = pal.TFs)

# Define node label size
# node.label.size <- c(3,5)[factor(V(net)$class)]
# node.label.size[rna.expr.i<0.05] <- 0

########################
## Define edge format ##
########################

# Define edge colors
E(net)$color <- "gray70"  # Colour negative correlation edges as blue
E(net)[which(E(net)$weight<0)]$color <- "blue"  # Colour negative correlation edges as blue
E(net)[which(E(net)$weight>0)]$color <- "red"   # Colour positive correlation edges as red

# Define edge width
edge_weight <- abs(E(net)$weight)   # Convert edge weights to absolute values
for (j in seq_len(length(E(net)))) {
  edge_weight[[j]] <- edge_weight[[j]]*(rna.expr.i[tail_of(net,j)] * rna.expr.i[head_of(net,j)])
}
edge_weight[edge_weight==0] <- 0.01
# edge_weight <- minmax.normalisation(edge_weight)

# alphas <- rna.expr.i
# alphas[rna.expr.i<0.5] <- 0.25
# alphas[rna.expr.i>0.75] <- 1

##################
## Plot network ##
##################

ggnet2(
  mode = layout, # "fruchtermanreingold",
  net = net,
  color = V(net)$color,
  edge.color = E(net)$color,
  edge.size = edge_weight,
  # node.size = c(15,20)[factor(V(net)$class)],
  # node.alpha = alphas,
  label = TRUE,
  label.color = "black",
  # label.size = node.label.size,
  arrow.size = 10,
  arrow.gap = 0.02,
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
