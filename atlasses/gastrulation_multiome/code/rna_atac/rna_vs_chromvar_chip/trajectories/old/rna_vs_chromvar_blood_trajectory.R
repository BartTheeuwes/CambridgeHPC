suppressMessages(library(scater))
suppressMessages(library(edgeR))
suppressMessages(library(argparse))

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--knn',    type="integer", default=25,    help='Number of kNN')
p$add_argument('--celltypes',  type="character",    nargs="+", help='Cell types')
p$add_argument('--motif_annotation',  type="character", default="Motif_cisbp", help='Motif annotation')
p$add_argument('--denoise', action="store_true", help='apply kNN denoising?')
p$add_argument('--outdir',   type="character",    help='Output directory')
args <- p$parse_args(commandArgs(TRUE))

## START TEST
args <- list()
args$celltypes = c("Haematoendothelial_progenitors", "Blood_progenitors_1", "Blood_progenitors_2", "Erythroid1", "Erythroid2", "Erythroid3")
args$motif_annotation <- "Motif_cisbp"
args$smooth <- TRUE
args$knn <- 25
args$outdir <- paste0(io$basedir,"/results/rna_atac/rna_vs_chromvar/trajectories/blood_trajectory")
## END TEST

########################
## Load ArchR Project ##
########################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/atac/archR/load_archR_project.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/atac/archR/load_archR_project.R")
} else {
  stop("Computer not recognised")
}


#####################
## Define settings ##
#####################

# I/O
io$pca.rna <- paste0(io$basedir,"/results/rna/dimensionality_reduction/all_cells/E7.5_rep1-E7.5_rep2-E8.0_rep1-E8.0_rep2-E8.5_rep1-E8.5_rep2_pca_features2500_pcs30_batchcorrectionbysample.txt.gz")
io$pca.atac <- paste0(io$basedir,"/results/atac/archR/dimensionality_reduction/PeakMatrix/all_cells/E7.5_rep1-E7.5_rep2-E8.0_rep1-E8.0_rep2-E8.5_rep1-E8.5_rep2_lsi_features50000_ndims30.txt.gz")
# io$outdir <- paste0(io$basedir,"/results/rna_atac/rna_vs_chromvar/trajectories/blood_trajectory")

dir.create(paste0(args$outdir,"/individual_genes"), showWarnings = F)

#####################
## Update metadata ##
#####################

sample_metadata <- fread(io$metadata) %>%
  .[pass_atacQC==TRUE & pass_rnaQC==TRUE & doublet_call==FALSE] %>%
  .[celltype.predicted%in%opts$celltypes & sample%in%opts$samples]

##################
## Subset ArchR ##
##################

ArchRProject.filt <- ArchRProject[sample_metadata$archR_cell,]

###############
## Load data ##
###############

# Load RNA-based trajectory
io$pseudotime <- paste0(io$basedir,"/results/rna/trajectories/blood_trajectory/blood_trajectory.txt.gz")
trajectory.dt <- fread(io$pseudotime)# %>%
  # merge(sample_metadata[,c("cell","archR_cell")]) %>% 
  # .[,cell:=NULL] %>% setnames("archR_cell","cell")

# Load highly variable along the trajectory
io$hvgs <- paste0(io$basedir,"/results/rna/trajectories/blood_trajectory/hvgs.rds")
hvgs <- readRDS(io$hvgs)

################################
## Load RNA and chromVAR data ##
################################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/rna_atac/rna_vs_chromvar/load_rna_chromvar_single_cells.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/rna_atac/rna_vs_chromvar/load_rna_chromvar_single_cells.R")
} else {
  stop("Computer not recognised")
}

# Merge
chromvar_rna_dt <- merge(
  rna_dt,
  chromvar_dt %>% 
    merge(sample_metadata[,c("cell","celltype.predicted")], by="cell"),
  by = c("cell","gene")
) %>% merge(trajectory.dt,by="cell")

length(unique(chromvar_rna_dt$gene))
length(unique(chromvar_rna_dt$cell))

fwrite(chromvar_rna_dt, sprintf("%s/chromvar_rna_blood_nodenoising.txt.gz",args$outdir), quote=F, sep="\t", na="NA")
# chromvar_rna_dt <- fread(sprintf("%s/chromvar_rna_blood.txt.gz",args$outdir))

######################################################
## Plot dynamics of individual TFs along pseudotime ##
######################################################

# genes.to.plot <- chromvar_rna_dt %>%
#   .[,.(expr=mean(expr),chromvar_zscore=mean(chromvar_zscore)),by=c("gene")] %>% 
#   .[expr>opts$min.expr,gene]

genes.to.plot <- unique(chromvar_rna_dt$gene)

for (i in genes.to.plot) {
  
  to.plot <- chromvar_rna_dt[gene==i] %>% 
    setnames("celltype.predicted","celltype") %>%
    # .[chromvar_zscore>7,chromvar_zscore:=7] %>%
    melt(id.vars=c("cell","V1","celltype"), measure.vars=c("chromvar_zscore","expr"), variable.name="modality")# %>%
    # .[,value_scaled:=value/max(value),by="modality"]
  # .[,value_scaled:=(value-min(value))/(max(value)-min(value)), by="modality"]
    
  p1 <- ggplot(to.plot, aes(x=V1, y=value)) +
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
  png(sprintf("%s/individual_genes/%s_rna_chromvar_vs_pseudotime_blood.png",args$outdir,i), width = 800, height = 600)
  print(p)
  dev.off()
}

#################################################
## Calculate directionality of the TF activity ##
#################################################

directionality.dt <- chromvar_rna_dt %>%
  .[,.(chromvar_zscore=mean(chromvar_zscore)),by=c("gene","celltype.predicted")] %>% 
  dcast(gene~celltype.predicted, value.var=c("chromvar_zscore")) %>%
  .[,diff:=Erythroid3-Haematoendothelial_progenitors] %>% 
  .[,sign:=c("Down","Up")[as.numeric(diff>0)+1]] %>%
  .[,c("gene","diff","sign")]

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
  .[, .(V1 = unlist(cor.test(chromvar_zscore, expr)[c("estimate", "p.value")])), by = c("gene")] %>%
  .[, para := rep(c("r","p"), .N/2)] %>% 
  data.table::dcast(gene ~ para, value.var = "V1") %>%
  .[,"padj_fdr" := list(p.adjust(p, method="fdr"))] %>%
  .[, sig := padj_fdr <= opts$threshold_fdr] %>% 
  setorder(padj_fdr, na.last = T)

# Add directionality info
cor.dt <- cor.dt %>% merge(directionality.dt,by=c("gene"))

# Volcano plot
negative_hits <- cor.dt[sig==TRUE & r<0,gene]
positive_hits <- cor.dt[sig==TRUE & r>0,gene]
all <- nrow(cor.dt)

xlim <- max(abs(cor.dt$r), na.rm=T)
ylim <- max(-log10(cor.dt$padj_fdr+1e-100), na.rm=T)

p <- ggplot(cor.dt, aes(x=r, y=-log10(padj_fdr+1e-100))) +
  # geom_hline(yintercept = -log10(opts$threshold_fdr), color="blue") +
  geom_segment(aes(x=0, xend=0, y=0, yend=ylim-1), color="orange", size=0.5) +
  ggrastr::geom_point_rast(aes(color=sig, size=sig)) +
  ggrepel::geom_text_repel(data=head(cor.dt[sig==T],n=40), aes(x=r, y=-log10(padj_fdr+1e-100), label=gene), size=3,  max.overlaps=50) +
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

pdf(sprintf("%s/volcano_pearson_correlation.pdf",args$outdir), width = 9, height = 6)
print(p)
dev.off()

fwrite(cor.dt, sprintf("%s/correlation_results_blood.txt.gz",args$outdir), quote=F, sep="\t", na="NA")

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
