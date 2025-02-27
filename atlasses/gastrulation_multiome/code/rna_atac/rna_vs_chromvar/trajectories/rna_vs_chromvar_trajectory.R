# suppressMessages(library(scater))
# suppressMessages(library(edgeR))
suppressMessages(library(argparse))

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--knn',    type="integer", default=25,    help='Number of kNN')
p$add_argument('--celltypes',  type="character",    nargs="+", help='Cell types')
p$add_argument('--motif_annotation',  type="character", default="Motif_cisbp", help='Motif annotation')
p$add_argument('--denoise', action="store_true", help='apply kNN denoising?')
p$add_argument('--trajectory',   type="character",    help='File with the trajectory')
p$add_argument('--trajectory_name',   type="character",    help='Name of the trajectory')
p$add_argument('--outdir',   type="character",    help='Output directory')
args <- p$parse_args(commandArgs(TRUE))

## START TEST
# args <- list()
# args$celltypes = c("Forebrain_Midbrain_Hindbrain", "Neural_crest")
# args$motif_annotation <- "Motif_cisbp"
# args$denoise <- TRUE
# args$knn <- 25
# args$trajectory <- "/Users/ricard/data/gastrulation_multiome_10x/results/rna/trajectories/neural_crest_trajectory/neural_crest_trajectory.txt.gz"
# args$trajectory_name <- "neural_crest"
# args$outdir <- "/Users/ricard/data/gastrulation_multiome_10x/results/rna_atac/rna_vs_chromvar/trajectories/neural_crest_test"
## END TEST

#####################
## Define settings ##
#####################

source(here::here("settings.R"))
source(here::here("utils.R"))

# I/O
# io$pca.rna <- paste0(io$basedir,"/results/rna/dimensionality_reduction/all_cells/E7.5_rep1-E7.5_rep2-E8.0_rep1-E8.0_rep2-E8.5_rep1-E8.5_rep2_pca_features2500_pcs30_batchcorrectionbysample.txt.gz")
# io$pca.atac <- paste0(io$basedir,"/results/atac/archR/dimensionality_reduction/PeakMatrix/all_cells/E7.5_rep1-E7.5_rep2-E8.0_rep1-E8.0_rep2-E8.5_rep1-E8.5_rep2_lsi_features50000_ndims50.txt.gz")
# io$outdir <- paste0(io$basedir,"/results/rna_atac/rna_vs_chromvar/trajectories/blood_trajectory")

dir.create(paste0(args$outdir,"/individual_genes"), showWarnings = F)

# Options
if (isFALSE(args$denoise)) args$knn <- 0

########################
## Load cell metadata ##
########################

sample_metadata <- fread(io$metadata) %>%
  .[pass_atacQC==TRUE & pass_rnaQC==TRUE & doublet_call==FALSE] %>%
  .[celltype.predicted%in%opts$celltypes & sample%in%opts$samples]

###############
## Load data ##
###############

# Load RNA-based trajectory
trajectory.dt <- fread(args$trajectory) %>% .[,PC1:=round(PC1,2)]

# Filter cells
cells <- intersect(sample_metadata$cell,trajectory.dt$cell)
sample_metadata <- sample_metadata[cell%in%cells]
trajectory.dt <- trajectory.dt[cell%in%cells] %>% setkey(cell) %>% .[cells]

################################
## Load RNA and chromVAR data ##
################################

io$rna_chromvar.file <- sprintf("%s/chromvar_rna.txt.gz",args$outdir)

if (file.exists(io$rna_chromvar.file)) {

  sprintf("Found precomputed values, loading %s...",io$rna_chromvar.file)
  chromvar_rna_dt <- fread(io$rna_chromvar.file)

} else {
  
  source(here::here("rna_atac/load_rna_atac_single_cells.R"))
  
  rm(list=c("rna.sce","atac.peakMatrix.se"))
  
  rna_tf.sce <- rna_tf.sce[,cells]
  atac.chromvar.se <- atac.chromvar.se[,cells]
  
  #################
  ## Filter data ##
  #################
  
  ## RNA ##
  
  # Remove non-expressed genes
  # rna.sce <- rna.sce[Matrix::rowMeans(logcounts(rna.sce))>0.01]
  # rna_tf.sce <- rna_tf.sce[Matrix::rowMeans(logcounts(rna_tf.sce))>0.01]
  
  # Subset genes that correlate with the trajectory
  # tmp <- cor(trajectory.dt$PC1,t(as.matrix(logcounts(rna.sce))))
  # rna.sce <- rna.sce[abs(tmp)>=0.15,]
  # tmp <- cor(trajectory.dt$PC1,t(as.matrix(logcounts(rna_tf.sce))))[1,]
  # rna_tf.sce <- rna_tf.sce[abs(tmp)>=0.15,]
  

  ## ATAC ##
  
  # atac.peakMatrix.se <- atac.peakMatrix.se[Matrix::rowMeans(assay(atac.peakMatrix.se))>0.01]
  # Subset peaks that correlate with the trajectory
  # tmp <- cor(trajectory.dt$PC1,t(as.matrix(assay(atac.peakMatrix.se))))
  # atac.peakMatrix.se <- atac.peakMatrix.se[abs(tmp)>=0.10,]
  
  #################
  ## Smooth data ##
  #################
  
  stopifnot(colnames(rna_tf.sce) == trajectory.dt$cell)
  stopifnot(colnames(atac.chromvar.se) == trajectory.dt$cell)
  
  if (args$denoise & args$knn>1) {
    
    # RNA
    # rna.mtx <- smoother_aggregate_nearest_nb(mat=as.matrix(logcounts(rna.sce)), D=pdist(matrix.please(trajectory.dt)), k=args$knn)
    pca.rna <- fread(io$pca.rna) %>% matrix.please %>% .[sample_metadata$cell,]
    rna_tf.mtx <- smoother_aggregate_nearest_nb(mat=as.matrix(logcounts(rna_tf.sce)), D=pdist(pca.rna), k=args$knn)
    # rna_tf.mtx <- smoother_aggregate_nearest_nb(mat=as.matrix(logcounts(rna_tf.sce)), D=pdist(matrix.please(trajectory.dt)), k=args$knn)

    # ATAC (chromVAR)
    pca.atac <- fread(io$pca.atac) %>% matrix.please %>% .[sample_metadata$cell,]
    atac.chromvar.mtx <- smoother_aggregate_nearest_nb(mat=as.matrix(assay(atac.chromvar.se,"z")), D=pdist(pca.atac), k=args$knn)
    # atac.chromvar.mtx <- smoother_aggregate_nearest_nb(mat=as.matrix(assay(atac.chromvar.se,"z")), D=pdist(matrix.please(trajectory.dt)), k=args$knn)

    # ATAC peaks
    # atac.peak.mtx <- smoother_aggregate_nearest_nb(mat=as.matrix(assay(atac.peakMatrix.se)), D=pdist(matrix.please(trajectory.dt)), k=args$knn)

  } else {
    # rna.mtx <- as.matrix(logcounts(rna.sce))
    rna_tf.mtx <- as.matrix(logcounts(rna_tf.sce))
    # atac.peak.mtx <- as.matrix(assay(atac.peakMatrix.se))
    atac.chromvar.mtx <- as.matrix(assay(atac.chromvar.se,"z"))
  }

  # colnames(rna.mtx) <- colnames(rna.sce)
  colnames(rna_tf.mtx) <- colnames(rna_tf.sce)
  # colnames(atac.peak.mtx) <- colnames(atac.peakMatrix.se)
  colnames(atac.chromvar.mtx) <- colnames(atac.chromvar.se)
  
  ##################################
  ## Create data.tables and merge ##
  ##################################
  
  chromvar.dt <- atac.chromvar.mtx %>% t %>%
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

  fwrite(chromvar_rna_dt, io$rna_chromvar.file, quote=F, sep="\t", na="NA")
}

# Load precomputed
chromvar_rna_dt <- fread(io$rna_chromvar.file)

# Add metadata
chromvar_rna_dt <- chromvar_rna_dt %>% 
  merge(sample_metadata[,c("cell","celltype.predicted")], by="cell") %>%
  merge(trajectory.dt,by="cell")

######################################################
## Plot dynamics of individual TFs along pseudotime ##
######################################################

# genes.to.plot <- chromvar_rna_dt %>%
#   .[,.(expr=mean(expr),chromvar_zscore=mean(chromvar_zscore)),by=c("gene")] %>% 
#   .[expr>opts$min.expr,gene]

genes.to.plot <- unique(chromvar_rna_dt$gene)# %>% head(n=3)

for (i in genes.to.plot) {
  
  to.plot <- chromvar_rna_dt[gene==i] %>% 
    setnames("celltype.predicted","celltype") %>%
    # .[chromvar_zscore>7,chromvar_zscore:=7] %>%
    melt(id.vars=c("cell","PC1","celltype"), measure.vars=c("chromvar_zscore","expr"), variable.name="modality")# %>%
    # .[,value_scaled:=value/max(value),by="modality"]
  # .[,value_scaled:=(value-min(value))/(max(value)-min(value)), by="modality"]
    
  p1 <- ggplot(to.plot, aes(x=PC1, y=value)) +
    # ggrastr::geom_point_rast(aes(fill=celltype), size=1.25, shape=21, stroke=0.1) +
    geom_point(aes(fill=celltype), size=1.75, shape=21, stroke=0.1) +
    stat_smooth(method="loess", color="black", alpha=0.75, span=0.5) +
    geom_rug(aes(color=celltype), sides="b") +
    facet_wrap(~modality, nrow=2, scales="free_y") +
    scale_color_manual(values=opts$celltype.colors) +
    scale_fill_manual(values=opts$celltype.colors) +
    guides(fill=F, color=F) +
    # scale_fill_manual(values=opts$celltype.colors) +
    # scale_fill_brewer(palette="Dark2") +
    labs(x="Pseudotime", y=i) +
    theme_classic() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.title = element_blank(),
      legend.position="top"
    )
  
  p2 <- ggboxplot(to.plot, x="celltype", y="value", fill="celltype", outlier.shape=NA) +
    scale_fill_manual(values=opts$celltype.colors) +
    labs(x="", y="") +
    # stat_compare_means(comparisons = list( c("Somitic_mesoderm", "NMP"), c("NMP", "Spinal_cord") ), label="p.signif", hide.ns=T) +
    facet_wrap(~modality, nrow=2, scales="free_y") +
    # coord_cartesian(ylim=c(0,5)) +
    theme_classic() +
    guides(x = guide_axis(angle = 90)) +
    theme(
      # axis.text.x = element_text(color="black", size=rel(0.75)),
      # legend.title = element_blank(),
      legend.position = "none",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  
  p <- cowplot::plot_grid(plotlist=list(p1,p2), nrow = 1, rel_widths = c(1/2,1/2))
  
  # pdf(sprintf("%s/individual_genes/%s_rna_chromvar_vs_pseudotime_blood.pdf",args$outdir,i), width=8, height=6)
  png(sprintf("%s/individual_genes/%s_rna_chromvar_vs_pseudotime_%s_knn%s.png",args$outdir,i,args$trajectory_name,args$knn), width = 800, height = 600)
  print(p)
  dev.off()
}

#################################################
## Calculate directionality of the TF activity ##
#################################################

tmp <- trajectory.dt %>% 
  merge(sample_metadata[,c("cell","celltype.predicted")],by="cell") %>%
  .[,.(PC1=mean(PC1)),by="celltype.predicted"] %>% setorder(-PC1)
celltype.start <- tmp[nrow(tmp),celltype.predicted]
celltype.end <- tmp[1,celltype.predicted]

directionality_chromvar.dt <- chromvar_rna_dt %>%
  .[,.(chromvar_zscore=mean(chromvar_zscore)),by=c("gene","celltype.predicted")] %>% 
  dcast(gene~celltype.predicted, value.var=c("chromvar_zscore")) %>%
  .[,c("gene",celltype.start,celltype.end), with=F] %>%
  setnames(c("gene","start","end")) %>%
  .[,chromvar_diff:=end-start] %>% 
  .[,chromvar_sign:=c("Down","Up")[as.numeric(chromvar_diff>0)+1]] %>%
  .[,c("gene","chromvar_diff","chromvar_sign")]

directionality_rna.dt <- chromvar_rna_dt %>%
  .[,.(expr=mean(expr)),by=c("gene","celltype.predicted")] %>% 
  dcast(gene~celltype.predicted, value.var=c("expr")) %>%
  .[,c("gene",celltype.start,celltype.end), with=F] %>%
  setnames(c("gene","start","end")) %>%
  .[,rna_diff:=end-start] %>% 
  .[,rna_sign:=c("Down","Up")[as.numeric(rna_diff>0)+1]] %>%
  .[,c("gene","rna_diff","rna_sign")]

##########################
## Correlation analysis ##
##########################

opts$threshold_fdr <- 0.10

# genes.to.test <- chromvar_rna_dt %>%
#   .[,.(expr=mean(expr),chromvar_zscore=mean(chromvar_zscore)),by=c("gene")] %>% 
#   .[expr>opts$min.expr,gene]

cor.dt <- chromvar_rna_dt %>%
  # .[gene%in%genes.to.test,] %>%
  .[,c("chromvar_zscore","expr"):=list(chromvar_zscore + rnorm(n=.N,mean=0,sd=1e-5), expr + rnorm(n=.N,mean=0,sd=1e-5))] %>% # add some noise 
  .[, .(PC1 = unlist(cor.test(chromvar_zscore, expr)[c("estimate", "p.value")])), by = c("gene")] %>%
  .[, para := rep(c("r","p"), .N/2)] %>% 
  data.table::dcast(gene ~ para, value.var = "PC1") %>%
  .[,"padj_fdr" := list(p.adjust(p, method="fdr"))] %>%
  .[, sig := padj_fdr <= opts$threshold_fdr]

# Add directionality info
cor.dt <- cor.dt %>% 
  merge(directionality_rna.dt,by="gene") %>%
  merge(directionality_chromvar.dt,by="gene") %>%
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
  geom_point(aes(color=sig, size=sig)) +
  ggrepel::geom_text_repel(data=head(cor.dt[sig==T],n=25), aes(x=r, y=-log10(padj_fdr+1e-100), label=gene), size=3,  max.overlaps=Inf) +
  scale_color_manual(values=c("black","red")) +
  scale_size_manual(values=c(0.75,1.25)) +
  scale_x_continuous(limits=c(-xlim-0.5,xlim+0.5)) +
  scale_y_continuous(limits=c(0,ylim+6)) +
  annotate("text", x=0, y=ylim+6, size=4, label=sprintf("(%d)", all)) +
  annotate("text", x=-xlim-0.5, y=ylim+6, size=4, label=sprintf("%d (-)",length(negative_hits))) +
  annotate("text", x=xlim+0.5, y=ylim+6, size=4, label=sprintf("%d (+)",length(positive_hits))) +
  labs(x="Pearson correlation", y=expression(paste("-log"[10],"(p.value)"))) +
  theme_classic() +
  theme(
    axis.text = element_text(size=rel(0.75), color='black'),
    axis.title = element_text(size=rel(1.0), color='black'),
    legend.position="none"
  )

pdf(sprintf("%s/volcano_correlation_%s_knn%s.pdf",args$outdir,args$trajectory_name,args$knn), width = 9, height = 6)
print(p)
dev.off()

fwrite(cor.dt, sprintf("%s/correlation_results.txt.gz",args$outdir), quote=F, sep="\t", na="NA")

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
