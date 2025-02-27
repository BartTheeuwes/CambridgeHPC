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
p$add_argument('--distance',  type="integer",            default=1e5,      help='Maximum distance for a linkage between a peak and a gene')
p$add_argument('--min_chip_threshold',  type="double",            default=0.25,      help='')
p$add_argument('--cor_threshold',  type="double",            default=0.30,      help='')
p$add_argument('--outdir',       type="character",                help='Output file')
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

## START TEST
args$distance <- 5e4
args$trajectory <- paste0(io$basedir,"/results/rna/trajectories/blood_trajectory/blood_trajectory.txt.gz")
args$trajectory_name <- "blood"
args$min_chip_threshold <- 0.25
args$cor_threshold <- 0.30
args$outdir <- paste0(io$basedir,"/results/rna_atac/gene_regulatory_networks/trajectories/blood")
## END TEST

#####################
## Load trajectory ##
#####################

trajectory.dt <- fread(args$trajectory)

###################
## Load metadata ##
###################

sample_metadata <- fread(io$metadata) %>%
  .[pass_atacQC==TRUE & pass_rnaQC==TRUE & doublet_call==FALSE]

##################
## Filter cells ##
##################

trajectory.dt <- trajectory.dt[cell%in%cells] 
cells <- intersect(trajectory.dt$cell,sample_metadata$cell)
sample_metadata <- sample_metadata[cell%in%cells] %>% setkey(cell) %>% .[cells]

# opts$celltypes.subset <- opts$celltypes[opts$celltypes%in%unique(sample_metadata$celltype.predicted)]
opts$celltypes_trajectory <- c("Haematoendothelial_progenitors", "Blood_progenitors_1", "Blood_progenitors_2", "Erythroid1", "Erythroid2", "Erythroid3")

#######################################
## Load virtual ChIP-seq annotation ###
#######################################

# io$virtual_chip.mtx <- file.path(io$basedir,"results_new/rna_atac/virtual_chipseq/CISBP/motifmatchr_virtual_chip.rds")
# virtual_chip.mtx <- readRDS(io$virtual_chip.mtx) %>% assay()
# virtual_chip.mtx[virtual_chip.mtx<=args$min_chip_threshold] <- 0

####################################
## Load pseudobulk RNA expression ##
####################################

io$rna.pseudobulk.sce <- file.path(io$basedir,"results_new/rna/pseudobulk/SingleCellExperiment_CISBP_pseudobulk_celltype.mapped_mnn.rds")
rna_tf_pseudobulk.sce <- readRDS(io$rna.pseudobulk.sce)[,opts$celltypes_trajectory]
# rna_tf_pseudobulk.sce <- rna_pseudobulk.sce[str_to_title(rownames(rna_tf.sce))]
# rownames(rna_tf_pseudobulk.sce) <- toupper(rownames(rna_tf_pseudobulk.sce))

####################################################
## Load TF2gene links based on in silico ChIP-seq ##
####################################################

tf2gene_chip.dt <- fread(io$tf2gene_virtual_chip) %>%
  .[chip_score>=opts$min_chip_score & dist<=opts$max_distance] %>% 
  .[,c("tf","gene")] %>% unique # Only keep TF-gene links

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

# Smooth single-cell data
pca.rna <- fread(io$pca.rna) %>% matrix.please %>% .[sample_metadata$cell,]
logcounts(rna_tf.sce) <- smoother_aggregate_nearest_nb(mat=as.matrix(logcounts(rna_tf.sce)), D=pdist(pca.rna), k=50)
logcounts(rna_genes.sce) <- smoother_aggregate_nearest_nb(mat=as.matrix(logcounts(rna_genes.sce)), D=pdist(pca.rna), k=50)

#######################
## Feature selection ##
#######################

rna.var <- apply(logcounts(rna_tf.sce),1,var) %>% sort 
TFs.rna <- rna.var[rna.var>0.05] %>% names

############################################################
## Select TFs with TF motif variability in the trajectory ##
############################################################

# chromvar_chip.mtx <- readRDS("/Users/argelagr/data/gastrulation_multiome_10x/results_new/atac/archR/chromvar_chip/pseudobulk/chromVAR_deviations_CISBP_archr_chip.rds")[,opts$celltypes_trajectory] %>% assay("z")
# chromvar_chip.var <- apply(chromvar_chip.mtx,1,var,na.rm=T) %>% sort
# TFs.atac <- chromvar_chip.var[chromvar_chip.var>=200] %>% names

chromvar_chip.mtx <- readRDS(file.path(io$basedir,"results_new/atac/archR/chromvar_chip/chromVAR_chip_deviations_CISBP_archr.rds"))[,sample_metadata$cell] %>% assay("z")
chromvar_chip.var <- apply(chromvar_chip.mtx,1,var,na.rm=T) %>% sort
TFs.atac <- chromvar_chip.var[chromvar_chip.var>=5] %>% names

#################
## Filter data ##
#################

TFs <- Reduce("intersect",list(TFs.rna,colnames(virtual_chip.mtx),TFs.atac))
# TFs <- c("TAL1","RUNX1","GATA1","KLF1","JUN")

rna_tf_filt.sce <- rna_tf.sce[TFs,]
rna_tf_pseudobulk_filt.sce <- rna_tf_pseudobulk.sce[TFs,]
virtual_chip_filt.mtx <- virtual_chip.mtx[,TFs]

peaks <- intersect(rownames(virtual_chip_filt.mtx), unique(peak2gene.dt[gene%in%TFs,peak]))
virtual_chip_filt.mtx <- virtual_chip_filt.mtx[peaks,]

# Print statistics
# print(sprintf("Number of TFs: %s",ncol(virtual_chip.mtx)))
# print(sprintf("Number of peaks: %s",nrow(virtual_chip.mtx)))
# print(sprintf("Number of genes: %s",nrow(rna_genes.sce)))

################################################################
## Calculate RNA expression correlation between TFs and genes ##
################################################################

# Compute correlations
corr.mtx <- psych::corr.test(t(logcounts(rna_tf.sce)),t(logcounts(rna_genes.sce)), ci = F)
diag(corr.mtx$r) <- NA

##############################
## Link TFs to target genes ##
##############################


# opts$pvalue.threshold <- 0.25

corr.mtx.filt <- corr.mtx
for (i in TFs) {
  
  # Select target peaks (note that we only take positive correlations into account)
  target_peaks_i <- names(which(virtual_chip.mtx[,i]>=args$min_chip_threshold))
  
  if (length(target_peaks_i)>=1) {
    
    # Select target genes
    target_genes_i <- intersect(
      x = names(which(abs(corr.mtx$r[i,])>=args$cor_threshold)),  # all correlations
      # x = names(which(corr.mtx$r[i,]>=opts$cor.threshold)),   # only positive correlations
      y = unique(peak2gene.dt[peak%in%target_peaks_i,gene])
    )
    
    print(sprintf("%s: %s target peaks & %s target genes",i,length(target_peaks_i),length(target_genes_i)))
    
    # Update correlation matrices
    corr.mtx.filt$p[i,!colnames(corr.mtx.filt$p)%in%target_genes_i] <- NA
    corr.mtx.filt$r[i,!colnames(corr.mtx.filt$p)%in%target_genes_i] <- NA
    
  } else {
    corr.mtx.filt$p[i,] <- NA
    corr.mtx.filt$r[i,] <- NA
  }
}

tmp <- corr.mtx.filt$r
diag(tmp) <- NA

# Filter edges
# tmp[abs(tmp)<=0.50] <- NA

# Filter nodes
# tmp <- tmp[rowSums(!is.na(tmp))>=1,]
# tmp <- tmp[,colSums(!is.na(tmp))>=1]

# Save
saveRDS(tmp, file.path(args$outdir,"tf2gene_matrix.rds"))

####################
## Create network ##
####################

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

# Create network
net <- graph_from_data_frame(d = edge_list.dt)

# ####################
# ## Create network ##
# ####################
# 
# # Create node and edge data.frames
# node_list.dt <- data.table(node_id=1:nrow(GRN_coef.mtx), node_name=rownames(GRN_coef.mtx))
# target_list.dt <- data.table(target_id=1:ncol(GRN_coef.mtx), target_name=colnames(GRN_coef.mtx))
# edge_list.dt <- GRN_coef.dt[,c("tf","gene","beta")] %>% setnames(c("from","to","weight"))
# node_list_metadata.dt <- data.table(
#   label = c(node_list.dt$node_name, target_list.dt$target_name),
#   class = c(rep("TF",nrow(node_list.dt)), rep("gene",nrow(target_list.dt)))
# )
# 
# # Create igraph object
# igraph.net <- graph_from_data_frame(d = edge_list.dt, vertices = node_list_metadata.dt)



