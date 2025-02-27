source(here::here("settings.R"))
source(here::here("utils.R"))

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--sce',        type="character",        help='SingleCellExperiment file')
p$add_argument('--metadata',    type="character",           help='Cell metadata')
p$add_argument('--chromvar',  type="character",   help='Motif annotation')
p$add_argument('--trajectory',       type="character",    help='File with the trajectory')
p$add_argument('--trajectory_name',   type="character",    help='Name of the trajectory')
p$add_argument('--outdir',          type="character",            help='Output directory')
args <- p$parse_args(commandArgs(TRUE))

## START TEST
args <- list()
io$basedir <- file.path(io$basedir,"test")
args$sce <- file.path(io$basedir,"results/rna/metacells/trajectories/nmp/SingleCellExperiment_metacells.rds")
args$metadata <- file.path(io$basedir,"results/rna/metacells/trajectories/nmp/metacells_metadata.txt.gz")
args$chromvar <- file.path(io$basedir,"results/atac/archR/chromvar_chip/metacells/chromVAR_chip_CISBP_archr.rds")
args$trajectory <- file.path(io$basedir,"results/rna/metacells/trajectories/nmp/metacell_trajectory.txt.gz")
args$trajectory_name <- "nmp"
args$outdir <- file.path(io$basedir,"results/rna_atac/rna_vs_chromvar_chip/metacells/trajectories")
## END TEST

#####################
## Define settings ##
#####################

#####################
## Load trajectory ##
#####################

trajectory.dt <- fread(args$trajectory) %>% setnames(c("metacell","V1","V2"))

########################
## Load cell metadata ##
########################

metadata.dt <- fread(args$metadata)

# Select cells
# cells <- intersect(trajectory.dt$cell,metadata.dt$cell)
# trajectory.dt <- trajectory.dt[cell%in%cells] 
# metadata.dt <- metadata.dt[cell%in%cells] %>% setkey(cell) %>% .[cells]

# Define cell types
opts$celltypes <- opts$celltypes[opts$celltypes%in%unique(metadata.dt$celltype)]

###################
## Load chromVAR ##
###################

chromvar.mtx <- readRDS(args$chromvar) %>% assay(.,"z")
chromvar.mtx <- chromvar.mtx[,colnames(chromvar.mtx) %in% metadata.dt$metacell]

# Cap chromVAR scores
# chromvar.mtx[chromvar.mtx>=20] <- 20
# chromvar.mtx[chromvar.mtx<=(-10)] <- (-10)
# hist(chromvar.mtx)

##############
## Load RNA ##
##############

rna.sce <- load_SingleCellExperiment(
  file = args$sce, 
  cells = metadata.dt$cell, 
  normalise = TRUE, 
  remove_non_expressed_genes = FALSE
)


# Extract RNA expression matrix for TFs
TFs <- intersect(toupper(rownames(rna.sce)),rownames(chromvar.mtx))
rna_tf.mtx <- logcounts(rna.sce)[str_to_title(TFs),] %>% as.matrix
rownames(rna_tf.mtx) <- toupper(rownames(rna_tf.mtx))

# Cap RNA expression
rna_tf.mtx[rna_tf.mtx>=4] <- 4

#################
## Smooth data ##
#################
  
stopifnot(colnames(rna_tf.mtx) == trajectory.dt$cell)
stopifnot(colnames(chromvar.mtx) == trajectory.dt$cell)

if (args$denoise & args$knn>1) {
  
  # RNA
  pca.rna <- fread(args$pca_rna) %>% matrix.please %>% .[metadata.dt$cell,]
  rna_tf.mtx <- smoother_aggregate_nearest_nb(mat=rna_tf.mtx, D=pdist(pca.rna), k=args$knn)
  # rna_tf.mtx <- smoother_aggregate_nearest_nb(mat=as.matrix(logcounts(rna_tf.mtx)), D=pdist(matrix.please(trajectory.dt)), k=args$knn)

  # ATAC chromVAR
  pca.atac <- fread(args$pca_atac) %>% matrix.please %>% .[metadata.dt$cell,]
  chromvar.mtx <- smoother_aggregate_nearest_nb(mat=chromvar.mtx, D=pdist(pca.atac), k=args$knn)
  # chromvar.mtx <- smoother_aggregate_nearest_nb(mat=as.matrix(assay(chromvar.mtx,"z")), D=pdist(matrix.please(trajectory.dt)), k=args$knn)
}
  
##################################
## Create data.tables and merge ##
##################################

chromvar.dt <- chromvar.mtx %>% t %>%
  as.data.table(keep.rownames = T) %>%
  setnames("rn","cell") %>%
  melt(id.vars=c("cell"), variable.name="gene", value.name="chromvar_zscore") %>%
  .[,chromvar_zscore:=round(chromvar_zscore,2)]

rna_tf.dt <- rna_tf.mtx %>%
  as.data.table(keep.rownames = T) %>%
  setnames("rn","gene") %>%
  melt(id.vars="gene", variable.name="cell", value.name="expr") %>%
  .[,expr:=round(expr,2)]

# Merge
chromvar_rna_dt <- merge(rna_tf.dt, chromvar.dt, by = c("cell","gene"))

# rm(rna_dt); rm(chromvar_dt)
length(unique(chromvar_rna_dt$gene))
length(unique(chromvar_rna_dt$cell))

# fwrite(chromvar_rna_dt, io$rna_chromvar.file, quote=F, sep="\t", na="NA")

# Load precomputed
# chromvar_rna_dt <- fread(io$rna_chromvar.file)

######################################################
## Plot dynamics of individual TFs along pseudotime ##
######################################################

to.plot <- chromvar_rna_dt %>%
  merge(metadata.dt[,c("cell","celltype.mapped_mnn")], by="cell") %>%
  setnames("celltype.mapped_mnn","celltype") %>%
  merge(trajectory.dt,by="cell") %>%
  .[!is.na(chromvar_zscore)] %>% .[chromvar_zscore<=0,chromvar_zscore:=0] %>%
  # .[,value_scaled:=(value-min(value))/(max(value)-min(value)), by="modality"]
  melt(id.vars=c("cell","PC1","celltype","gene"), measure.vars=c("chromvar_zscore","expr"), variable.name="modality") %>%
  .[,celltype:=factor(celltype,levels=opts$celltypes)] %>%
  .[,modality:=factor(modality,levels=c("expr","chromvar_zscore"))] %>%
  .[,value_scaled:=value/max(value),by="modality"]

facet.labels <- c(expr = "RNA expression", chromvar_zscore = "Motif accessibility")


genes.to.plot <- to.plot %>% .[modality=="expr",.(var(value_scaled)),by="gene"] %>% .[V1>0.005,gene]
# genes.to.plot <- unique(to.plot$gene)
# genes.to.plot <- c("TAL1", "GATA1", "JUN", "RUNX1","GATA2", "KLF1", "PBX1")

for (i in genes.to.plot) {
    
  p1 <- ggplot(to.plot[gene==i], aes(x=PC1, y=value_scaled)) +
    # ggrastr::geom_point_rast(aes(fill=celltype), size=1.25, shape=21, stroke=0.1) +
    geom_point(aes(fill=celltype), size=1.75, shape=21, stroke=0.1) +
    stat_smooth(method="loess", color="black", alpha=0.75, span=0.5) +
    geom_rug(aes(color=celltype), sides="b") +
    facet_wrap(~modality, nrow=2, scales="free_y",  labeller = as_labeller(facet.labels)) +
    scale_color_manual(values=opts$celltype.colors) +
    scale_fill_manual(values=opts$celltype.colors) +
    guides(fill="none", color="none") +
    coord_cartesian(ylim=c(-0.02,1.02)) +
    # scale_fill_manual(values=opts$celltype.colors) +
    # scale_fill_brewer(palette="Dark2") +
    labs(x="Pseudotime", y=i) +
    theme_classic() +
    theme(
      strip.text = element_text(size=rel(1.25), color="black"),
      axis.text.x = element_blank(),
      axis.text.y = element_text(color="black"),
      axis.ticks.x = element_blank(),
      legend.title = element_blank(),
      legend.position="none"
    )
  
  p2 <- ggplot(to.plot[gene==i], aes_string(x="celltype", y="value_scaled", fill="celltype")) +
    geom_jitter(aes(color=celltype), size=0.5, alpha=0.5, width = 0.2) +
    geom_boxplot(outlier.shape=NA, alpha=0.75) +
    facet_wrap(~modality, nrow=2, scales="free_y",  labeller = as_labeller(facet.labels)) +
    coord_cartesian(ylim=c(-0.02,1.02)) +
    scale_fill_manual(values=opts$celltype.colors) +
    scale_color_manual(values=opts$celltype.colors) +
    labs(x="", y="") +
    theme_classic() +
    guides(x = guide_axis(angle = 90)) +
    theme(
      strip.text = element_text(size=rel(1.1), color="black"),
      legend.position = "none",
      axis.text.x = element_blank(),
      axis.text.y = element_text(color="black"),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  
  p <- cowplot::plot_grid(plotlist=list(p1,p2), nrow = 1, rel_widths = c(1/2,1/2))
  
  # png(file.path(args$outdir,sprintf("individual_genes/%s_rna_chromvar_chip_vs_pseudotime_%s_knn%s.png",i,args$trajectory_name,args$knn)), width = 800, height = 500)
  pdf(file.path(args$outdir,sprintf("individual_genes/%s_rna_chromvar_chip_vs_pseudotime_%s_knn%s.pdf",i,args$trajectory_name,args$knn)), width = 10, height = 6)
  print(p)
  dev.off()
}

##########################
## Correlation analysis ##
##########################

stop("FILTER GENES BY VARIABILITY BEFORE CORRELATION ANAYSUS")

# genes.to.test <- chromvar_rna_dt %>%
#   .[,.(expr=mean(expr),chromvar_zscore=mean(chromvar_zscore)),by=c("gene")] %>% 
#   .[expr>opts$min.expr,gene]

cor.dt <- copy(chromvar_rna_dt) %>%
  .[!is.na(chromvar_zscore)] %>%
  .[,c("chromvar_zscore","expr"):=list(chromvar_zscore + rnorm(n=.N,mean=0,sd=1e-5), expr + rnorm(n=.N,mean=0,sd=1e-5))] %>% # add some noise 
  .[, .(PC1 = unlist(cor.test(chromvar_zscore, expr)[c("estimate", "p.value")])), by = c("gene")] %>%
  .[, para := rep(c("r","p"), .N/2)] %>% 
  data.table::dcast(gene ~ para, value.var = "PC1") %>%
  .[,"padj_fdr" := list(p.adjust(p, method="fdr"))] %>%
  .[, sig := padj_fdr<=0.01] %>%
  setorder(padj_fdr, na.last = T)

# Volcano plot
negative_hits <- cor.dt[sig==TRUE & r<0,gene]
positive_hits <- cor.dt[sig==TRUE & r>0,gene]
all <- nrow(cor.dt)

xlim <- max(abs(cor.dt$r), na.rm=T)
ylim <- max(-log10(cor.dt$padj_fdr+1e-100), na.rm=T)

p <- ggplot(cor.dt, aes(x=r, y=-log10(padj_fdr+1e-100))) +
  # geom_hline(yintercept = -log10(opts$threshold_fdr), color="blue") +
  geom_segment(aes(x=0, xend=0, y=0, yend=ylim-1), color="orange", size=0.5) +
  # ggrastr::geom_point_rast(aes(color=sig, size=sig)) +
  geom_jitter(aes(color=abs(r), alpha=abs(r)), width = 0.02, height=2) +
  ggrepel::geom_text_repel(data=head(cor.dt[r>0],n=15), aes(x=r, y=-log10(padj_fdr+1e-100), label=gene), size=3,  max.overlaps=Inf) +
  ggrepel::geom_text_repel(data=head(cor.dt[r<0],n=5), aes(x=r, y=-log10(padj_fdr+1e-100), label=gene), size=3,  max.overlaps=Inf) +
  scale_color_gradient(low = "gray80", high = "red") +
  scale_alpha_continuous(range=c(0.25,1)) +
  # scale_size_manual(values=c(0.75,1.25)) +
  scale_x_continuous(limits=c(-xlim-0.15,xlim+0.15)) +
  scale_y_continuous(limits=c(0,ylim+6)) +
  annotate("text", x=0, y=ylim+6, size=4, label=sprintf("(%d)", all)) +
  annotate("text", x=-xlim, y=5, size=4, label=sprintf("%d (-)",length(negative_hits))) +
  annotate("text", x=xlim, y=5, size=4, label=sprintf("%d (+)",length(positive_hits))) +
  labs(x="Pearson correlation", y=expression(paste("-log"[10],"(p.value)"))) +
  theme_classic() +
  theme(
    axis.text = element_text(size=rel(0.75), color='black'),
    axis.title = element_text(size=rel(1.0), color='black'),
    legend.position="none"
  )

pdf(file.path(args$outdir,sprintf("volcano_correlation_%s_knn%s.pdf",args$trajectory_name,args$knn)), width = 9, height = 6)
print(p)
dev.off()

fwrite(cor.dt, file.path(args$outdir,sprintf("correlation_results.txt.gz")), quote=F, sep="\t", na="NA")

##########
## Test ##
##########

# to.plot <- chromvar_rna_dt[gene=="Sox18"] 
# 
# ggscatter(to.plot, x="expr", y="chromvar_zscore", color="celltype.predicted", size=1) +
#   stat_smooth(method="lm", color="black", alpha=0.75, span=0.5) +
#   scale_color_manual(values=opts$celltype.colors) +
#   theme(
#     axis.text = element_text(size=rel(0.75)),
#     legend.position = "none"
#   )
