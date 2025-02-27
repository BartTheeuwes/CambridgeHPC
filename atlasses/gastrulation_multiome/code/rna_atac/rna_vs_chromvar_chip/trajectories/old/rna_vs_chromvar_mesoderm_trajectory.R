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
io$outdir <- paste0(io$basedir,"/results/rna_atac/rna_vs_chromvar/trajectories/mesoderm_trajectory")

# Options

# opts$celltypes = c(
#   "Epiblast",
#   "Primitive_Streak",
#   "Nascent_mesoderm",
#   "Mixed_mesoderm",
#   "Intermediate_mesoderm"
#   # "Caudal_Mesoderm",
#   # "Paraxial_mesoderm",
#   # "Somitic_mesoderm",
#   # "Pharyngeal_mesoderm"
# )

opts$stage_celltype = c(
  "E7.5_Epiblast",
  "E7.5_Primitive_Streak",
  "E7.5_Nascent_mesoderm",
  "E7.5_Mixed_mesoderm",
  # "E8.5_Mixed_mesoderm",
  "E8.5_Intermediate_mesoderm"
)

opts$stage_celltype.colors <- c(
  "E7.5_Epiblast" = "#635547",
  "E7.5_Primitive_Streak" = "#DABE99",
  "E7.5_Nascent_mesoderm" = "#C594BF",
  "E7.5_Mixed_mesoderm" = "#DFCDE4",
  "E8.5_Intermediate_mesoderm" = "#139992"
)

opts$min.expr <- 0.1

#####################
## Update metadata ##
#####################

sample_metadata <- fread(io$metadata) %>%
  .[pass_atacQC==TRUE & pass_rnaQC==TRUE] %>%
  # .[celltype.mapped%in%opts$celltypes & sample%in%opts$samples]
  .[,stage_celltype:=paste(stage,celltype.mapped,sep="_")] %>%
  .[stage_celltype%in%opts$stage_celltype] %>%
  .[,stage_celltype:=factor(stage_celltype,levels=opts$stage_celltype)]

##################
## Subset ArchR ##
##################

ArchRProject.filt <- ArchRProject[sample_metadata$archR_cell,]

###############
## Load data ##
###############

# Load RNA-based trajectory
io$pseudotime <- "/Users/ricard/data/gastrulation_multiome_10x/results/rna/trajectories/mesoderm_trajectory/mesoderm_trajectory.txt.gz"
trajectory.dt <- fread(io$pseudotime)# %>%
  # merge(sample_metadata[,c("cell","archR_cell")]) %>% 
  # .[,cell:=NULL] %>% setnames("archR_cell","cell")

# Load highly variable along the NMP trajectory
io$hvgs <- "/Users/ricard/data/gastrulation_multiome_10x/results/rna/trajectories/mesoderm_trajectory/hvgs.rds"
hvgs <- readRDS(io$hvgs)

###############################
## Load SingleCellExperiment ##
###############################

sce <- load_SingleCellExperiment(
  file = io$sce, 
  cells = sample_metadata$cell, 
  normalise = TRUE, 
  remove_non_expressed_genes = FALSE
)

##########################
## Load chromVAR scores ##
##########################

atac.deviation.se <- readRDS(io$archR.deviations.se)
atac.deviation.mtx <- atac.deviation.se %>%
  .[,colnames(atac.deviation.se)%in%sample_metadata$archR_cell] %>% 
  assay(.,"z")
dim(atac.deviation.mtx)

#################
## Smooth data ##
#################

opts$knn <- 15
opts$npcs <- 5

# RNA
# foo <- trajectory.dt %>% tibble::column_to_rownames("cell") %>% as.matrix
pca.rna <- irlba::prcomp_irlba(t(as.matrix(logcounts(sce[hvgs,]))), n=opts$npcs)$x
rownames(pca.rna) <- colnames(sce)
rna.matrix.smoothed <- smoother_aggregate_nearest_nb(mat=as.matrix(logcounts(sce)), D=pdist(pca.rna), k=opts$knn)
colnames(rna.matrix.smoothed) <- colnames(sce)

# ATAC
pca.atac <- irlba::prcomp_irlba(t(atac.deviation.mtx), n=opts$npcs)$x
rownames(pca.atac) <- colnames(atac.deviation.mtx)
atac.deviation.mtx.smoothed <- smoother_aggregate_nearest_nb(mat=atac.deviation.mtx, D=pdist(pca.atac), k=opts$knn)
colnames(atac.deviation.mtx.smoothed) <- colnames(atac.deviation.mtx)

################################
## Load motif2gene annotation ##
################################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/atac/archR/load_motif_annotation.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/atac/archR/load_motif_annotation.R")
} else {
  stop("Computer not recognised")
}

motif2gene.dt <- motif2gene.dt %>%
  .[motif%in%rownames(atac.deviation.se)] %>% setkey(motif) %>% .[rownames(atac.deviation.se)] %>%
  .[gene%in%rownames(sce)] %>%
  .[,N:=length(unique(motif)),by="gene"] %>% .[N==1] %>% .[,N:=NULL]

sce <- sce[motif2gene.dt$gene,]
atac.deviation.mtx <- atac.deviation.mtx %>% .[motif2gene.dt$motif,]

################
## parse data ##
################

# chromvar_dt <- atac.deviation.mtx %>% 
chromvar_dt <- atac.deviation.mtx.smoothed %>% 
  .[motif2gene.dt$motif,] %>%
  as.matrix %>% t %>% as.data.table(keep.rownames = T) %>%
  setnames("rn","archR_cell") %>% 
  melt(id.vars=c("archR_cell"), variable.name="motif", value.name="chromvar_zscore")

# rna_dt <- sce %>%
#   .[motif2gene.dt$gene,] %>%
#   logcounts %>% as.matrix %>% 
rna_dt <- rna.matrix.smoothed[motif2gene.dt$gene,] %>%
  as.data.table(keep.rownames = T) %>%
  setnames("rn","gene") %>%
  melt(id.vars="gene", variable.name="cell", value.name="expr")

############################
## Merge RNA and chromVAR ##
############################

chromvar_rna_dt <- merge(
  rna_dt,
  chromvar_dt %>% 
    merge(sample_metadata[,c("cell","archR_cell","celltype.mapped","stage_celltype")], by="archR_cell") %>%
    merge(motif2gene.dt[,c("motif","gene")],by="motif") %>% .[,archR_cell:=NULL],
  by = c("cell","gene")
) %>% merge(trajectory.dt,by="cell")

fwrite(chromvar_rna_dt, sprintf("%s/chromvar_rna_mesoderm.txt.gz",io$outdir), quote=F, sep="\t", na="NA")

##########
## Plot ##
##########

genes.to.plot <- chromvar_rna_dt %>%
  .[,.(expr=mean(expr),chromvar_zscore=mean(chromvar_zscore)),by=c("gene","motif")] %>% 
  .[expr>opts$min.expr,gene]

# genes.to.plot <- unique(chromvar_rna_dt$gene) %>% head(n=3)

for (i in genes.to.plot) {
  
  to.plot <- chromvar_rna_dt[gene==i] %>% 
    setnames("celltype.mapped","celltype") %>%
    # .[chromvar_zscore>7,chromvar_zscore:=7] %>%
    melt(id.vars=c("cell","V1","celltype","stage_celltype"), measure.vars=c("chromvar_zscore","expr"), variable.name="modality")# %>%
    # .[,value_scaled:=value/max(value),by="modality"]
  # .[,value_scaled:=(value-min(value))/(max(value)-min(value)), by="modality"]
    
  p1 <- ggplot(to.plot, aes(x=V1, y=value)) +
    geom_point(aes(fill=celltype), size=1.25, shape=21, stroke=0.1) +
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
  
  p2 <- ggboxplot(to.plot, x="stage_celltype", y="value", fill="stage_celltype", outlier.shape=NA) +
    scale_fill_manual(values=opts$stage_celltype.colors) +
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
  
  pdf(sprintf("%s/individual_genes/%s_rna_chromvar_vs_pseudotime_mesoderm.pdf",io$outdir,i), width=8, height=6)
  print(p)
  dev.off()
}

#################################################
## Calculate directionality of the TF activity ##
#################################################

directionality.dt <- chromvar_rna_dt %>%
  .[,.(chromvar_zscore=mean(chromvar_zscore)),by=c("gene","motif","stage_celltype")] %>% 
  dcast(gene+motif~stage_celltype, value.var=c("chromvar_zscore")) %>%
  .[,diff:=E8.5_Intermediate_mesoderm-E7.5_Epiblast] %>% 
  .[,sign:=c("Down","Up")[as.numeric(diff>0)+1]] %>%
  .[,c("motif","gene","diff","sign")]

##########################
## Correlation analysis ##
##########################

opts$threshold_fdr <- 0.10

genes.to.test <- chromvar_rna_dt %>%
  .[,.(expr=mean(expr),chromvar_zscore=mean(chromvar_zscore)),by=c("gene","motif")] %>% 
  .[expr>opts$min.expr,gene]

cor.dt <- chromvar_rna_dt[gene%in%genes.to.test,] %>%
  .[,c("chromvar_zscore","expr"):=list(chromvar_zscore + rnorm(n=.N,mean=0,sd=1e-5), expr + rnorm(n=.N,mean=0,sd=1e-5))] %>% # add some noise 
  .[, .(V1 = unlist(cor.test(chromvar_zscore, expr)[c("estimate", "p.value")])), by = c("gene","motif")] %>%
  .[, para := rep(c("r","p"), .N/2)] %>% 
  data.table::dcast(gene+motif ~ para, value.var = "V1") %>%
  .[,"padj_fdr" := list(p.adjust(p, method="fdr"))] %>%
  .[, sig := padj_fdr <= opts$threshold_fdr] %>% 
  setorder(padj_fdr, na.last = T)

# Add directionality info
cor.dt <- cor.dt %>% merge(directionality.dt,by=c("gene","motif"))

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

pdf(sprintf("%s/volcano_plots/volcano_pearson_correlation.pdf",io$outdir), width = 9, height = 6)
# png(sprintf("%s/volcano_plots/volcano_pearson_correlation.png",io$outdir), width = 800, height = 500)
print(p)
dev.off()

fwrite(cor.dt, sprintf("%s/correlation_results_mesoderm.txt.gz",io$outdir), quote=F, sep="\t", na="NA")

##########
## Test ##
##########

# to.plot <- chromvar_rna_dt[gene=="Sox18"] 
# 
# ggscatter(to.plot, x="expr", y="chromvar_zscore", color="celltype.mapped", size=1) +
#   stat_smooth(method="lm", color="black", alpha=0.75, span=0.5) +
#   scale_color_manual(values=opts$celltype.colors) +
#   theme(
#     axis.text = element_text(size=rel(0.75)),
#     legend.position = "none"
#   )
