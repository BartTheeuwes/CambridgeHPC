source(here::here("settings.R"))
source(here::here("utils.R"))

suppressPackageStartupMessages(library(furrr))

options(future.globals.maxSize = 5000 * 1024^2)

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--distance',  type="integer",            default=1e5,      help='Maximum distance for a linkage between a peak and a gene')
p$add_argument('--ncores',  type="integer",            default=1,      help='Number of cores')
p$add_argument('--min_chip_threshold',  type="double",            default=0.25,      help='')
# p$add_argument('--cor_threshold',  type="double",            default=0.30,      help='')
p$add_argument('--outdir',       type="character",                help='Output file')
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

## START TEST
args$distance <- 5e4
args$trajectory <- paste0(io$basedir,"/results_new/rna/trajectories/blood/blood_trajectory.txt.gz")
args$trajectory_name <- "blood"
args$min_chip_threshold <- 0.25
# args$cor_threshold <- 0.30
args$ncores <- 1
args$outdir <- paste0(io$basedir,"/results_new/rna_atac/gene_regulatory_networks/trajectories/blood")
## END TEST

# I/O
io$tf2gene_virtual_chip <- file.path(io$basedir,"results_new/rna_atac/gene_regulatory_networks/pseudobulk/TF2gene_after_virtual_chip.txt.gz")

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

cells <- intersect(trajectory.dt$cell,sample_metadata$cell)
trajectory.dt <- trajectory.dt[cell%in%cells] 
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
  .[chip_score>=args$min_chip_threshold & dist<=args$distance] %>% 
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

# Extract TFs
TFs <- intersect(unique(tf2gene_chip.dt$tf),toupper(rownames(rna.sce)))
rna_tf.sce <- rna.sce[rownames(rna.sce)%in%str_to_title(TFs),]
rownames(rna_tf.sce) <- toupper(rownames(rna_tf.sce))

# Extract genes
genes <- intersect(unique(tf2gene_chip.dt$gene),rownames(rna.sce))
rna_genes.sce <- rna.sce[rownames(rna.sce)%in%genes,]

# Smooth single-cell data
# pca.rna <- fread(io$pca.rna) %>% matrix.please %>% .[cells,]
# logcounts(rna_tf.sce) <- smoother_aggregate_nearest_nb(mat=as.matrix(logcounts(rna_tf.sce)), D=pdist(pca.rna), k=50)
# logcounts(rna_genes.sce) <- smoother_aggregate_nearest_nb(mat=as.matrix(logcounts(rna_genes.sce)), D=pdist(pca.rna), k=50)

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

TFs <- Reduce("intersect",list(TFs.rna,unique(tf2gene_chip.dt$tf),TFs.atac))
# TFs <- c("TAL1","RUNX1","GATA1","KLF1","JUN")

rna_tf_filt.sce <- rna_tf.sce[TFs,]
rna_tf_pseudobulk_filt.sce <- rna_tf_pseudobulk.sce[TFs,]
tf2gene_chip.dt <- tf2gene_chip.dt[tf%in%TFs]

# Print statistics
# print(sprintf("Number of TFs: %s",ncol(virtual_chip.mtx)))
# print(sprintf("Number of peaks: %s",nrow(virtual_chip.mtx)))
# print(sprintf("Number of genes: %s",nrow(rna_genes.sce)))

#################################################
## Building GRN using linear regression models ##
#################################################

if (args$ncores>1){
  plan(multicore, workers=args$ncores)
  genes_split <- split(genes, cut(seq_along(genes), args$ncores, labels = FALSE)) 
} else {
  plan(sequential)
  genes_split <- list(genes)
}

GRN_lm.dt <- genes_split %>% future_map(function(genes) {
  tmp <- tf2gene_chip.dt[gene%in%genes]
  genes %>% map(function(i) {
    # print(sprintf("%s (%d/%d)",i,match(i,genes),length(genes)))
    tfs <- tmp[gene==i,tf]
    tfs %>% map(function(j) {
      x <- logcounts(rna_tf_filt.sce)[j,]
      y <- logcounts(rna_genes.sce)[i,]
      lm.fit <- lm(y~x)
      df <- data.frame(tf=j, gene=i, beta=round(coef(lm.fit)[[2]],3), pvalue=summary(lm.fit)$coefficients[2,4])
    }) %>% rbindlist
  }) %>% rbindlist %>% return
}) %>% rbindlist

#################################################
## Building GRN using correlation coefficients ##
#################################################

# Compute correlations
corr.mtx <- psych::corr.test(t(as.matrix(logcounts(rna_tf.sce))),t(as.matrix(logcounts(rna_genes.sce))), ci = F)
diag(corr.mtx$r) <- NA

GRN_cor.dt <- genes_split %>% future_map(function(genes) {
  tmp <- tf2gene_chip.dt[gene%in%genes]
  genes %>% map(function(i) {
    
    corr.mtx$r[,genes]
    
    # print(sprintf("%s (%d/%d)",i,match(i,genes),length(genes)))
    tfs <- tmp[gene==i,tf]
    tfs %>% map(function(j) {
      x <- logcounts(rna_tf_filt.sce)[j,]
      y <- logcounts(rna_genes.sce)[i,]
      lm.fit <- lm(y~x)
      df <- data.frame(tf=j, gene=i, cor=round(corr.mtx$r[j,i]), pvalue=corr.mtx$p[j,i])
    }) %>% rbindlist
  }) %>% rbindlist %>% return
}) %>% rbindlist

#####################################################
## Create network from linear regression estimates ##
#####################################################

args$beta_threshold <- 0.30
tmp <- GRN_lm.dt[abs(beta)>=args$beta_threshold]

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
lm.net <- graph_from_data_frame(d = edge_list.dt)


#####

args$cor_threshold <- 0.30

##########
## Save ##
##########

fwrite(GRN_cor.dt, file.path(args$outdir,"TF2gene_cor.tsv.gz"), sep="\t", quote=F, na="NA")
saveRDS(cor.net, file.path(args$outdir,"TF2gene_network_cor.rds"))

fwrite(GRN_lm.dt, file.path(args$outdir,"TF2gene_lm.tsv.gz"), sep="\t", quote=F, na="NA")
saveRDS(lm.net, file.path(args$outdir,"TF2gene_network_lm.rds"))
