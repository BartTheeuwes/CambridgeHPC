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
# p$add_argument('--sce',       type="character",                help='SingleCellExperiment')
p$add_argument('--outdir',       type="character",                help='Output file')
p$add_argument('--peak2gene_links',       type="character",                help='')
p$add_argument('--tf2tf_cor',       type="character",                help='')
p$add_argument('--virtual_chip_mtx',       type="character",                help='')
p$add_argument('--motif_annotation',       type="character",                help='Motif annotation')
args <- p$parse_args(commandArgs(TRUE))

# cor.threshold

#####################
## Define settings ##
#####################

# load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

## START TEST
args$distance <- 2e4
args$motif_annotation <- "CISBP"
args$tf2tf_cor <- file.path(io$basedir,"results_new/rna/coexpression/correlation_matrix_tf2tf_pseudobulk.rds")
args$virtual_chip_mtx <- file.path(io$basedir,"results_new/rna_atac/virtual_chipseq/CISBP/virtual_chip.mtx")
args$peak2gene_links <- file.path(io$basedir,"results_new/atac/archR/peak_calling/peaks2genes/peaks2genes_all.txt.gz")
args$outdir <- paste0(io$basedir,"/results_new/rna_atac/gene_regulatory_networks/pseudobulk/all_TFs")
## END TEST

# I/O
dir.create(args$outdir, showWarnings = F)

# Options
# opts$TFs <- list.files(io$virtual_chip.dir, pattern = "*.bed.gz") %>% stringr::str_replace_all(".bed.gz",""))

##############################
## Load pseudobulk RNA data ##
##############################

# Load SingleCellExperiment
sce.pseudobulk <- readRDS(io$rna.pseudobulk.sce)

# Define TFs
# opts$TFs <- opts$TFs[str_to_title(opts$TFs)%in%rownames(sce.pseudobulk)]

#############################
## Load tf2tf correlations ##
#############################

tf2tf_cor.mtx <- readRDS(args$tf2tf_cor)#[opts$TFs,opts$TFs]

#############################
## Load peak2gene linkages ##
#############################

peak2gene.dt <- fread(io$archR.peak2gene.all) %>% 
  .[dist<=args$distance]
  # .[gene%in%opts$genes] %>%
  # .[,peak:=sprintf("chr%s:%s-%s",chr,peak.start,peak.end)]

################################
## Load virtual ChIP-seq data ##
################################

virtual_chip.mtx <- readRDS(io$virtual_chip.mtx)#[,opts$TFs]  

################
## Parse data ##
################

TFs <- intersect(colnames(virtual_chip.mtx),colnames(tf2tf_cor.mtx))
virtual_chip.mtx <- virtual_chip.mtx[,TFs]
tf2tf_cor.mtx <- tf2tf_cor.mtx[TFs,TFs]

peak2gene.dt <- peak2gene.dt %>% .[,gene:=toupper(gene)] %>% .[gene%in%TFs]

##############################################################
## Use the virtual ChIP-seq data to find TF-TF associations ##
##############################################################

opts$cor.threshold <- 0.40
opts$chip.threshold <- 0.25

tf2tf_cor_chip.mtx <- tf2tf_cor.mtx

for (i in rownames(tf2tf_cor_chip.mtx)) {

  # Select target peaks
  target_peaks_i <- names(which(virtual_chip.mtx[,i]>=opts$chip.threshold))
  
  # Select target genes
  target_genes_i <- intersect(
    x = names(which(abs(tf2tf_cor_chip.mtx[i,])>=opts$cor.threshold)),  # all correlations
    # x = names(which(tf2tf_cor_chip.mtx$r[i,]>=opts$cor.threshold)),   # only positive correlations
    y = unique(peak2gene.dt[peak%in%target_peaks_i,gene])
  )
  
  print(sprintf("%s: %s target genes",i,length(target_genes_i)))
  
  # Update correlation matrix
  tf2tf_cor_chip.mtx[i,!colnames(tf2tf_cor_chip.mtx)%in%target_genes_i] <- NA
}

# Filter TFs and genes with too little connections
tf2tf_cor_chip.mtx <- tf2tf_cor_chip.mtx[rowSums(!is.na(tf2tf_cor_chip.mtx))>=1,,drop=F]
tf2tf_cor_chip.mtx <- tf2tf_cor_chip.mtx[,colSums(!is.na(tf2tf_cor_chip.mtx))>=1,drop=F]

# opts$TFs <- rownames(tf2tf_cor_chip.mtx)

sum(!is.na(tf2tf_cor_chip.mtx))

####################
## Create network ##
####################

# Prepare data
node_list.dt <- data.table(node_id=1:nrow(tf2tf_cor_chip.mtx), node_name=rownames(tf2tf_cor_chip.mtx))
target_list.dt <- data.table(target_id=1:ncol(tf2tf_cor_chip.mtx), target_name=colnames(tf2tf_cor_chip.mtx))

edge_list.dt <- as.data.table(tf2tf_cor_chip.mtx,keep.rownames = T) %>% 
  setnames("rn","from") %>%
  melt(id.vars=c("from"), variable.name="to", value.name="weight") %>%
  .[!is.na(weight)]

# Filter out self-regulatory events
edge_list.dt <- edge_list.dt[from!=to]

# net <- graph_from_data_frame(d = edge_list.dt, vertices = node_list_metadata.dt)
net <- graph_from_data_frame(d = edge_list.dt)

##################
## Save network ##
##################

saveRDS(tf2tf_cor_chip.mtx, file.path(args$outdir,sprintf("%s_connectivity_matrix_all_TFs.rds",args$motif_annotation)))
saveRDS(net, file.path(args$outdir,sprintf("%s_network_all_TFs.rds",args$motif_annotation)))

# save connectivity matrix as sparse matrix
# mtx_to_save[is.na(mtx_to_save)] <- 0
# mtx_to_save <- Matrix::Matrix(mtx_to_save)
# Matrix::writeMM(mtx_to_save, file.path(args$outdir,sprintf("%s_all_TFs.mtx",args$motif_annotation)))