suppressMessages(library(argparse))

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--tf2gene_matrix',       type="character",                help='Output file')
p$add_argument('--trajectory',       type="character",                help='Output file')
p$add_argument('--trajectory_name',       type="character",                help='Output file')
p$add_argument('--outdir',       type="character",                help='Output file')
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

# load default setings
source(here::here("settings.R"))
source(here::here("utils.R"))

## START TEST
args$tf2gene_matrix <- file.path(io$basedir,"results/rna_atac/gene_regulatory_networks/trajectories/blood/tf2gene_matrix.rds")
args$trajectory <- file.path(io$basedir,"results/rna/trajectories/blood_trajectory/blood_trajectory.txt.gz")
args$trajectory_name <- "blood"
args$outdir <- file.path(io$basedir,"results/rna_atac/gene_regulatory_networks/trajectories/blood/TF_vs_gene_expr")
## END TEST

dir.create(args$outdir, showWarnings = F)

#####################
## Load trajectory ##
#####################

trajectory.dt <- fread(args$trajectory)

###################
## Load metadata ##
###################

sample_metadata <- fread(io$metadata) %>%
  .[pass_atacQC==TRUE & pass_rnaQC==TRUE & doublet_call==FALSE]

cells <- intersect(trajectory.dt$cell,sample_metadata$cell)
trajectory.dt <- trajectory.dt[cell%in%cells] 
sample_metadata <- sample_metadata[cell%in%cells] %>% setkey(cell) %>% .[cells]

opts$celltypes.subset <- opts$celltypes[opts$celltypes%in%unique(sample_metadata$celltype.predicted)]

#########################
## Load TF2gene matrix ##
#########################

tf2gene.mtx <- readRDS(args$tf2gene_matrix)

# write.table(tf2gene.mtx, file.path(io$basedir,"results/rna_atac/gene_regulatory_networks/trajectories/blood/tf2gene_matrix.csv.gz"), sep=",", quote=F, col.names=T, row.names=T)

# tmp <- tf2gene.mtx
# tmp[!is.na(tmp)] <- 1
# tmp[is.na(tmp)] <- 0
# table(tmp)
# rownames(tmp) <- str_to_title(rownames(tmp))
# write.table(tmp, "/Users/argelagr/data/gastrulation_multiome_10x/results/rna_atac/GRN/chip_GRN_ricard.csv", sep=",", quote=T, col.names=T, row.names=T)
# 
# tmp["KLF1","H19"]

#####################################
## Load single-cell RNA expression ##
#####################################

rna.sce <- load_SingleCellExperiment(
  file = io$rna.sce, 
  cells = sample_metadata$cell, 
  normalise = TRUE, 
  remove_non_expressed_genes = TRUE
)

rna_tf.sce <- rna.sce[rownames(rna.sce)%in%str_to_title(rownames(tf2gene.mtx)),]
rownames(rna_tf.sce) <- toupper(rownames(rna_tf.sce))

rna_genes.sce <- rna.sce[rownames(rna.sce)%in%colnames(tf2gene.mtx),]

##############################################
## Load single-cell chromatin accessibility ##
##############################################

io$archR.GeneScoreMatrix.se <- file.path(io$basedir,"results/atac/archR/gene_scores/GeneScoreMatrix_tss.rds")

atac_gene_scores.se <- readRDS(io$archR.GeneScoreMatrix.se)
dim(atac_gene_scores.se)
[sample_metadata$cell,]


rna_tf.sce <- rna.sce[rownames(rna.sce)%in%str_to_title(rownames(tf2gene.mtx)),]
rownames(rna_tf.sce) <- toupper(rownames(rna_tf.sce))

rna_genes.sce <- rna.sce[rownames(rna.sce)%in%colnames(tf2gene.mtx),]


#############################
## Smooth single-cell data ##
#############################

pca.rna <- fread(io$pca.rna) %>% matrix.please %>% .[sample_metadata$cell,]
logcounts(rna_tf.sce) <- smoother_aggregate_nearest_nb(mat=as.matrix(logcounts(rna_tf.sce)), D=pdist(pca.rna), k=50)
logcounts(rna_genes.sce) <- smoother_aggregate_nearest_nb(mat=as.matrix(logcounts(rna_genes.sce)), D=pdist(pca.rna), k=50)

#######################################
## Plot TF vs target gene expression ##
#######################################

TFs.to.plot <- rownames(tf2gene.mtx)

# i <- "KLF1"
for (i in TFs.to.plot) {
  genes.to.plot <- names(which(!is.na(tf2gene.mtx[i,])))
  stopifnot(genes.to.plot%in%rownames(rna_genes.sce))
  
  # j <- "Nars"
  for (j in genes.to.plot) {
    
    to.plot <- data.table(
      cell = colnames(rna_tf.sce),
      TF = logcounts(rna_tf.sce[i,])[1,],
      gene = logcounts(rna_genes.sce[j,])[1,]
    ) %>% merge(sample_metadata[,c("cell","celltype.predicted")]) %>% 
      setnames("celltype.predicted","celltype") %>%
      merge(trajectory.dt[,c("cell","PC1")],by="cell") %>%
      melt(id.vars=c("cell","PC1","celltype"), measure.vars=c("TF","gene"), variable.name="class") %>%
      .[,value:=minmax.normalisation(value),by="class"]
    
    p <- ggplot(to.plot, aes(x=PC1, y=value, group=class)) +
      geom_point(aes(fill=class), size=1.25, shape=21, stroke=0.1) +
      stat_smooth(aes(color=class), method="loess",alpha=1, span=0.5) +
      # geom_rug(aes(color=celltype), sides="b") +
      # scale_color_manual(values=opts$celltype.colors) +
      guides(color="none") +
      scale_fill_discrete(labels = c(i,j)) +
      labs(x="Pseudotime", y="RNA expression") +
      theme_classic() +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        legend.position="top"
      )
    
    pdf(file.path(args$outdir,sprintf("%s_vs_%s_%s_trajectory.pdf",i,j,args$trajectory_name)))
    print(p)
    dev.off()
    
  }
    
  
}