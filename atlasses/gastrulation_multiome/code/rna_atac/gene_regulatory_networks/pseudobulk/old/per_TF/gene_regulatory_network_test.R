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
# p$add_argument('--celltypes',    type="character",    nargs="+",  help='Cell type')
# p$add_argument('--remove_ExE_celltypes', action="store_true",   help='Remove ExE cell types?')
# p$add_argument('--motif_annotation',    type="character",   default="Motif_cisbp",      help='Motif annotation')
# p$add_argument('--denoised',     action="store_true",             help='Use denoised RNA data to calculate the coexpression matrix?')
# p$add_argument('--knn',  type="integer",            default=15,      help='Number of kNN')
p$add_argument('--genes',    type="character", nargs="+",    help='TFs to plot (Eomes, Foxa2, etc.)')
p$add_argument('--distance',  type="integer",            default=1e5,      help='Maximum distance for a linkage between a peak and a gene')
p$add_argument('--stringent',     action="store_true",             help='Filter peaks by correlation with TF?')
p$add_argument('--outdir',       type="character",                help='Output file')
args <- p$parse_args(commandArgs(TRUE))

## START TEST
args$distance <- 1e5
args$genes <- c("Mixl1")
args$stringent <- TRUE
args$outdir <- "/Users/ricard/data/gastrulation_multiome_10x/results/rna_atac/gene_regulatory_networks/pseudobulk/per_TF/MIXL1"
## END TEST

#####################
## Define settings ##
#####################

# load default setings
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
# io$archR.pseudobulk.deviations.se <- sprintf("%s/pseudobulk/pseudobulk_DeviationMatrix_%s_summarized_experiment.rds",io$archR.directory,opts$motif_annotation)
# io$archR.pseudobulk.deviations.se <- sprintf("%s/results/atac/archR/chromvar/pseudobulk/chromVAR_deviations_summarized_experiment_%s_pseudobulk_correlated_peaks.rds",io$basedir,opts$motif_annotation)
# io$pca.rna <- paste0(io$basedir,"/results/rna/dimensionality_reduction/all_cells/E7.5_rep1-E7.5_rep2-E8.0_rep1-E8.0_rep2-E8.5_rep1-E8.5_rep2_pca_features2500_pcs30_batchcorrectionbysample.txt.gz")

# Options
opts$celltypes = c(
  "Epiblast",
  "Primitive_Streak",
  "Caudal_epiblast",
  # "PGC",
  # "Anterior_Primitive_Streak",
  "Notochord",
  "Def._endoderm",
  "Gut",
  "Nascent_mesoderm",
  "Mixed_mesoderm",
  "Intermediate_mesoderm",
  "Caudal_Mesoderm",
  "Paraxial_mesoderm",
  "Somitic_mesoderm",
  "Pharyngeal_mesoderm",
  "Cardiomyocytes",
  "Allantois",
  "ExE_mesoderm",
  "Mesenchyme",
  "Haematoendothelial_progenitors",
  "Endothelium",
  "Blood_progenitors_1",
  "Blood_progenitors_2",
  "Erythroid1",
  "Erythroid2",
  "Erythroid3",
  "NMP",
  "Rostral_neurectoderm",
  # "Caudal_neurectoderm",
  "Neural_crest",
  "Forebrain_Midbrain_Hindbrain",
  "Spinal_cord",
  "Surface_ectoderm",
  "Visceral_endoderm",
  "ExE_endoderm",
  "ExE_ectoderm",
  "Parietal_endoderm"
)
# Sanity checks
args$genes <- str_to_title(args$genes)


##############################
## Load pseudobulk RNA data ##
##############################

sce.pseudobulk <- readRDS(io$rna.pseudobulk.sce)[,opts$celltypes]

# Filter genes manually
# sort(grep("Rik$|^Gm|^Rps|^Rpl|^mt-",rownames(sce.pseudobulk), value = T))
sce.pseudobulk.filt <- sce.pseudobulk[!grepl("Rik$|^Gm|^Rps|^Rpl|^mt-",rownames(sce.pseudobulk)),]

# Filter to marker genes
marker_genes.dt <- fread(io$rna.atlas.marker_genes)
sce.pseudobulk.filt <- sce.pseudobulk.filt[rownames(sce.pseudobulk.filt)%in%unique(marker_genes.dt$gene),]

# Select variable genes
sce.pseudobulk.filt <- sce.pseudobulk.filt[apply(logcounts(sce.pseudobulk.filt),1,var)>0.01,]

###############################
## Load motifmatcher results ##
###############################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/load_motifmatchR.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/load_motifmatchR.R")
} else {
  stop("Computer not recognised")
}

###################################################
## Load TF2peak correlation results (pseudobulk) ##
###################################################

# io$tf2peak_cor.dt <- paste0(io$basedir,"/results/rna_atac/rna_vs_acc/pseudobulk/TFexpr_vs_peakAcc/cor_TFexpr_vs_peakAcc.txt.gz")
tf2peak_cor.dt <- fread(io$tf2peak_cor.dt) %>%
  .[!is.na(cor)] %>%
  .[,cor_sign:=c("-","+")[(cor>0)+1]]

# io$tf2peak_cor.se <- paste0(io$basedir,"/results/rna_atac/rna_vs_acc/pseudobulk/TFexpr_vs_peakAcc/cor_TFexpr_vs_peakAcc_SummarizedExperiment.rds")
tf2peak_cor.se <- readRDS(io$tf2peak_cor.se)

#############################
## Load peak2gene linkages ##
#############################

# peak2gene.dt <- fread(io$archR.peak2gene.all) %>%
peak2gene.dt <- fread(io$archR.peak2gene.nearest) %>%
  .[gene%in%rownames(sce.pseudobulk.filt)] %>%
  .[,peak:=sprintf("chr%s:%s-%s",chr,peak.start,peak.end)]

####################################################################################
## Calculate RNA expression correlation between TFs and target genes (pseudobulk) ##
####################################################################################

# TO-DO: INCLUDE CORRELATIONS BETWEEN TFS

stopifnot(args$genes %in% rownames(sce.pseudobulk.filt))
stopifnot(toupper(args$genes) %in% colnames(motifmatcher.se))
stopifnot(toupper(args$genes) %in% colnames(tf2peak_cor.se))

rna.mtx.TFs <- logcounts(sce.pseudobulk)[args$genes,,drop=F]; rownames(rna.mtx.TFs) <- toupper(rownames(rna.mtx.TFs))
rna.mtx.target_genes <- logcounts(sce.pseudobulk.filt)[!rownames(sce.pseudobulk.filt)%in%args$genes,]

TFs <- rownames(rna.mtx.TFs)
target_genes <- rownames(rna.mtx.target_genes)

# Compute correlations
corr.mtx <- psych::corr.test(t(rna.mtx.TFs),t(rna.mtx.target_genes), ci = F)

###########################################################
## Select peaks based on the target genes computed above ##
###########################################################

peak2gene.dt.filt <- peak2gene.dt[gene%in%rownames(rna.mtx.target_genes) & dist<args$distance]
peaks <- unique(peak2gene.dt.filt$peak)

motifmatcher.se.filt <- motifmatcher.se[peaks,TFs]
tf2peak_cor.se.filt <- tf2peak_cor.se[peaks,]
tf2peak_cor.dt.filt <- tf2peak_cor.dt[peak%in%peaks]

####################################################################################################
## Filter peaks by correlation between TF's RNA expression and the peak's chromatin accessibility ##
####################################################################################################

opts$cor.threshold <- 0.40
# opts$pvalue.threshold <- 0.25

corr.mtx.filt <- corr.mtx
corr.mtx.filt <- corr.mtx
for (i in TFs) {
  
  #################
  ## Prune peaks ##
  #################
  
  # Filter peaks by RNA expr & ATAC cor
  # target_peaks_i <- which(dropNA2matrix(assay(tf2peak_cor.se.filt[,i],"pvalue"))[,1]<=opts$pvalue.threshold) %>% names
  # target_peaks_i <- which(dropNA2matrix(assay(tf2peak_cor.se.filt[,i],"cor"))[,1]>=opts$cor.threshold) %>% names
  target_peaks_i <- which(abs(dropNA2matrix(assay(tf2peak_cor.se.filt[,i],"cor")))[,1]>=opts$cor.threshold) %>% names
  
  # Filter peaks by motif presence (no RNAexpr & ATACcor)
  target_peaks_noacc_i <- names(which(assay(motifmatcher.se.filt[,i])[,1]))
  
  print(sprintf("%s: %s/%s target peaks after filtering based on TF expression vs peak accessibility",i,length(target_peaks_i),length(target_peaks_noacc_i)))
  
  ########################
  ## Prune target genes ##
  ########################
  
  putative_target_genes_i <- intersect(
    # x = names(which(corr.mtx$p[i,]<=opts$pvalue.threshold)),  
    # x = names(which(corr.mtx$r[i,]>=opts$cor.threshold)),  
    x = names(which(abs(corr.mtx$r[i,])>=opts$cor.threshold)),  
    y = unique(peak2gene.dt.filt[peak%in%target_peaks_i,gene])
  )
  
  putative_target_genes_noacc_i <- intersect(
    x = names(which(corr.mtx$r[i,]>=opts$cor.threshold)),  
    # x = names(which(corr.mtx$p[i,]<=opts$pvalue.threshold)),  
    y = unique(peak2gene.dt.filt[peak%in%target_peaks_noacc_i,gene])
    )
  
  print(sprintf("%s: %s/%s target genes after filtering based on TF expression vs peak accessibility",i,length(putative_target_genes_i),length(putative_target_genes_noacc_i)))
  
  #################################
  ## Update correlation matrices ##
  #################################
  
  if (args$stringent) {
    corr.mtx.filt$p[i,!colnames(corr.mtx.filt$p)%in%putative_target_genes_i] <- NA
    corr.mtx.filt$r[i,!colnames(corr.mtx.filt$p)%in%putative_target_genes_i] <- NA
  } else {
    corr.mtx.filt$p[i,!colnames(corr.mtx.filt$p)%in%putative_target_genes_noacc_i] <- NA
    corr.mtx.filt$r[i,!colnames(corr.mtx.filt$p)%in%putative_target_genes_noacc_i] <- NA
  }
}

####################
## Create network ##
####################

# TO-DO: FILTER CASES WITH TOO LITTLE GENES

tmp <- corr.mtx.filt$r
# tmp[abs(tmp)<0.25] <- NA
tmp <- tmp[,colSums(!is.na(tmp))>=1,drop=F]

# Save
if (args$stringent) {
  saveRDS(tmp, sprintf("%s/corr_matrix_stringent.rds",args$outdir))
} else {
  saveRDS(tmp, sprintf("%s/corr_matrix_lenient.rds",args$outdir))
}

# Prepare data
node_list.dt <- data.table(node_id=1:nrow(tmp), node_name=rownames(tmp))
target_list.dt <- data.table(target_id=1:ncol(tmp), target_name=colnames(tmp))

edge_list.dt <- as.data.table(tmp,keep.rownames = T) %>% 
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

set.seed(42)
# layout <- layout.fruchterman.reingold(net)
layout <- layout.sphere(net)

###########################################################################
## Plot network, colouring each gene by celltype-specific RNA expression ##
###########################################################################

# Load PAGA graph
if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/load_paga_graph.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/load_paga_graph.R")
} else {
  stop("Computer not recognised")
}

# Define colors
pal <- grDevices::colorRamp(c("gray60", "purple"))( (1:100)/100 )
pal <- cbind(pal, seq(100, 255, length.out = 100))

rna.mtx.TFs.scaled <- apply(rna.mtx.TFs,1,minmax.normalisation) %>% t
rna.mtx.target_genes.scaled <- apply(rna.mtx.target_genes,1,minmax.normalisation) %>% t
# rna.mtx.target_genes.scaled <- scale(rna.mtx.target_genes, center=TRUE, scale = TRUE)

# celltypes.to.plot <- c("Epiblast", "Notochord", "Gut", "Erythroid3")
celltypes.to.plot <- opts$celltypes

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
  foo <- rna.mtx.TFs.scaled[,i]; names(foo) <- rownames(rna.mtx.TFs)
  rna.expr.i <- c(foo, rna.mtx.target_genes.scaled[,i])[names(V(net))]# %>% minmax.normalisation
  
  # Discretise a bit
  rna.expr.i[rna.expr.i<0.33] <- 0
  rna.expr.i[rna.expr.i>0.66] <- 1
  
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
  edge_weight[edge_weight==0] <- 0.01
  # edge_weight <- minmax.normalisation(edge_weight)
  
  alphas <- rna.expr.i
  # alphas[rna.expr.i<0.5] <- 0.25
  # alphas[rna.expr.i>0.75] <- 1
  
  p2 <- ggnet2(
    mode = layout, # "fruchtermanreingold",
    net = net,
    color = V(net)$color,
    edge.color = E(net)$color,
    edge.size = edge_weight,
    node.size = c(5,25)[factor(V(net)$class)],
    node.alpha = alphas,
    label = TRUE,
    label.color = "black",
    label.size = c(3,5)[factor(V(net)$class)],
    arrow.size = 0.1,
    legend.position = "none"
  ) 
  
  ###################
  ## Combine plots ##
  ###################
  
  p <- cowplot::plot_grid(plotlist=list(p1,p2), rel_widths = c(2/5,3/5), nrow=1)
    
  # pdf(sprintf("%s/%s_%s_test.pdf",args$outdir,paste(args$genes,collapse="-"),i), width=13, height=6, bg="white")
  if (args$stringent) {
    outfile <- sprintf("%s/%s_%s_stringent.png",args$outdir,paste(args$genes,collapse="-"),i)
  } else {
    outfile <- sprintf("%s/%s_%s_lenient.png",args$outdir,paste(args$genes,collapse="-"),i)
  }
  png(outfile, width = 1200, height = 550, bg = "white")
  print(p)
  dev.off()
}


################
## Create GIF ##
################

# library(magick)
# 
# celltypes.to.plot = c(
#   "Epiblast",
#   "Primitive_Streak",
#   "Def._endoderm",
#   "Notochord",
#   "Gut",
#   "Nascent_mesoderm",
#   "Intermediate_mesoderm",
#   "Somitic_mesoderm",
#   "Pharyngeal_mesoderm",
#   "Cardiomyocytes",
#   "Allantois",
#   "ExE_mesoderm",
#   "Haematoendothelial_progenitors",
#   "Endothelium",
#   "Blood_progenitors_1",
#   "Blood_progenitors_2",
#   "Erythroid1",
#   "Erythroid2",
#   "Erythroid3"
# )
# 
# sprintf("%s/%s_%s_test.png",args$outdir,paste(args$genes,collapse="-"),celltypes.to.plot) %>%
#   map(image_read) %>%
#   image_join %>%
#   image_animate(fps=1) %>%
#   image_write(quality=100, path = sprintf("%s/test.gif",args$outdir))



###########################
## Test interactive plot with ggiraph ##
###########################

# library(ggiraph)
# 
# i <- "Gut"
# 
# colors <- opts$celltype.colors[opts$celltypes]
# 
# alphas <- rep(1.0,length(opts$celltypes)); names(alphas) <- opts$celltypes
# 
# sizes <- rep(9,length(opts$celltypes)); names(sizes) <- opts$celltypes
# 
# p1 <- ggnet2_interactive(
#   net = net.paga,
#   mode = c("x", "y"),
#   color = colors,
#   node.alpha = alphas,
#   node.size = sizes,
#   edge.size = 0.15,
#   edge.color = "grey",
#   label = TRUE,
#   label.size = 2.5,
#   legend.position = "none"
# )
# 
# 
# 
# # htmlwidget call
# x <- girafe(
#   ggobj = p1, 
#   width_svg = 6, height_svg = 6,
#   options = list(
#     opts_sizing(rescale = FALSE),
#     # opts_tooltip(
#     #   opacity = .8,
#     #   css = "background-color:white;color:black;padding:2px;border-radius:2px;"
#     # ),
#     opts_hover_inv(css = "opacity:0.65;"),
#     # opts_hover(css = "fill:#1279BF;stroke:#1279BF;cursor:pointer;r:15px")
#     opts_hover(css = "cursor:pointer;r:17px")
#   )
# )
# x
# 
