
#####################
## Define settings ##
#####################

source(here::here("settings.R"))
source(here::here("utils.R"))

# Options
opts$trajectory_name <- "blood"
opts$knn <- 50

# I/O
io$trajectory <- sprintf("%s/results/rna/trajectories/%s_trajectory/%s_trajectory.txt.gz",io$basedir,opts$trajectory_name,opts$trajectory_name)
io$outdir <- sprintf("%s/results/rna_atac/rna_vs_chromvar/trajectories/%s_trajectory_knn%s/explore",io$basedir,opts$trajectory_name,opts$knn); dir.create(io$outdir, showWarnings = F)

#####################
## Load metadata ##
#####################

sample_metadata <- fread(io$metadata) %>%
  .[pass_atacQC==TRUE & pass_rnaQC==TRUE & doublet_call==FALSE] %>%
  .[celltype.predicted%in%opts$celltypes & sample%in%opts$samples]

###############
## Load data ##
###############

# Load trajectory
trajectory.dt <- fread(io$trajectory) %>% .[,PC1:=round(PC1,2)]

# Filter cells
cells <- intersect(sample_metadata$cell,trajectory.dt$cell)
sample_metadata <- sample_metadata[cell%in%cells]
trajectory.dt <- trajectory.dt[cell%in%cells] %>% setkey(cell) %>% .[cells]

#################################################
## Load precomputed RNA and chromVAR estimates ##
#################################################

io$rna_chromvar.file <- sprintf("%s/results/rna_atac/rna_vs_chromvar/trajectories/%s_trajectory_knn%s/chromvar_rna.txt.gz",io$basedir,opts$trajectory_name,opts$knn)
chromvar_rna_dt <- fread(io$rna_chromvar.file) %>% 
  # merge(sample_metadata[,c("cell","celltype.predicted")], by="cell") %>%
  # merge(trajectory.dt,by="cell") %>%
  setnames("celltype.predicted","celltype")

##########################
## Correlation analysis ##
##########################

# Load precomputed correlation estimates
io$cor_rna_vs_chromvar <- sprintf("%s/results/rna_atac/rna_vs_chromvar/trajectories/%s_trajectory_knn%s/correlation_results.txt.gz",io$basedir,opts$trajectory_name,opts$knn)
cor.dt <- fread(io$cor_rna_vs_chromvar)


cor.dt %>% 
  .[abs(r)<0.25,sig:=F] %>%
  .[,log_p:=-log10(padj_fdr)] %>%
  .[is.infinite(log_p),log_p:=runif(.N,min=275, max=350)] %>%
  .[,dot_size:=minmax.normalisation(log_p)]
  # .[padj_fdr==1e-275,padj_fdr:=sample(size=.N)]

# Volcano plot
negative_hits <- cor.dt[sig==TRUE & r<0,gene]
positive_hits <- cor.dt[sig==TRUE & r>0,gene]
all <- nrow(cor.dt)

xlim <- max(abs(cor.dt$r), na.rm=T)
ylim <- max(cor.dt$log_p, na.rm=T)

p <- ggplot(cor.dt, aes(x=r, y=log_p)) +
  geom_segment(aes(x=0, xend=0, y=0, yend=ylim-1), color="orange", size=0.5) +
  geom_jitter(aes(fill=sig, size=dot_size), shape=21, color="black", width=0.1, height=0.05, alpha=0.8) +
  ggrepel::geom_text_repel(data=cor.dt[r>0.75][sample(.N,6)], aes(x=r, y=log_p, label=gene), size=4,  max.overlaps=Inf) +
  ggrepel::geom_text_repel(data=cor.dt[r<(-0.6)][sample(.N,6)], aes(x=r, y=log_p, label=gene), size=4,  max.overlaps=Inf) +
  scale_fill_manual(values=c("gray30","red")) + 
  scale_size_continuous(range = c(0.25,2)) +
  scale_x_continuous(limits=c(-xlim-0.5,xlim+0.5)) +
  scale_y_continuous(limits=c(0,ylim+25)) +
  # annotate("text", x=0, y=ylim+25, size=4, label=sprintf("(%d)", all)) +
  annotate("text", x=-xlim-0.3, y=25, size=4, label=sprintf("%d (-)",length(negative_hits))) +
  annotate("text", x=xlim+0.3, y=25, size=4, label=sprintf("%d (+)",length(positive_hits))) +
  labs(x="Pearson correlation", y=expression(paste("-log"[10],"(p.value)"))) +
  guides(size="none") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size=rel(0.9), color='black'),
    axis.text.y = element_text(size=rel(0.75), color='black'),
    axis.title = element_text(size=rel(1.0), color='black'),
    legend.position="none",
    legend.title = element_blank()
  )

pdf(sprintf("%s/volcano_correlation_%s_knn%s.pdf",io$outdir,opts$trajectory_name,opts$knn), width = 7, height = 6)
print(p)
dev.off()


######################################################
## Plot dynamics of individual TFs along pseudotime ##
######################################################

facet.labels <- c(expr = "RNA expression", chromvar_zscore = "Motif accessibility (chromVAR+)")

# genes.to.plot <- unique(chromvar_rna_dt$gene)# %>% head(n=3)
# genes.to.plot <- cor.dt[sig==T & abs(r)>=0.5 & rna_sign=="Up",gene]

genes.to.plot <- c("TAL1", "GATA1", "KLF1", "RUNX1", "JUN")  # Blood
genes.to.plot <- c("POU3F1", "MYB", "MYT1L", "OTX2 ", "NANOG", "SP8", "T", "FOXB1", "LEF1", "FOXC1","GATA4") # mesoderm
c("POU3F1", "NFIB", "BCL11A", "MYT1L", "ARID3B", "ELF3", "PRDM1", "FOXA1", "LHX1", "GATA5", "PAX9", "RFX6") # endoderm
for (i in genes.to.plot) {
  
  to.plot <- chromvar_rna_dt[gene==i] %>% 
    melt(id.vars=c("cell","PC1","celltype"), measure.vars=c("chromvar_zscore","expr"), variable.name="modality") %>%
    .[,modality:=factor(modality,levels=c("expr","chromvar_zscore"))]
    # .[,value_scaled:=value/max(value),by="modality"]
  # .[,value_scaled:=(value-min(value))/(max(value)-min(value)), by="modality"]
    
  p1 <- ggplot(to.plot, aes(x=PC1, y=value)) +
    # ggrastr::geom_point_rast(aes(fill=celltype), size=1.25, shape=21, stroke=0.1) +
    geom_point(aes(fill=celltype), size=1.75, shape=21, stroke=0.1) +
    stat_smooth(method="loess", color="black", alpha=0.75, span=0.5) +
    geom_rug(aes(color=celltype), sides="b") +
    facet_wrap(~modality, nrow=2, scales="free_y",  labeller = as_labeller(facet.labels)) +
    scale_color_manual(values=opts$celltype.colors) +
    scale_fill_manual(values=opts$celltype.colors) +
    guides(fill=F, color=F) +
    labs(x="", y="") +
    theme_classic() +
    theme(
      axis.text.y = element_text(size=rel(0.8), color="black"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.title = element_blank(),
      legend.position="top"
    )
  
  # p2 <- ggboxplot(to.plot, x="celltype", y="value", fill="celltype", outlier.shape=NA) +
  #   scale_fill_manual(values=opts$celltype.colors) +
  #   labs(x="", y="") +
  #   # stat_compare_means(comparisons = list( c("Somitic_mesoderm", "NMP"), c("NMP", "Spinal_cord") ), label="p.signif", hide.ns=T) +
  #   facet_wrap(~modality, nrow=2, scales="free_y") +
  #   # coord_cartesian(ylim=c(0,5)) +
  #   theme_classic() +
  #   guides(x = guide_axis(angle = 90)) +
  #   theme(
  #     # axis.text.x = element_text(color="black", size=rel(0.75)),
  #     # legend.title = element_blank(),
  #     legend.position = "none",
  #     axis.text.x = element_blank(),
  #     axis.title.x = element_blank(),
  #     axis.ticks.x = element_blank()
  #   )
  
  # p <- cowplot::plot_grid(plotlist=list(p1,p2), nrow = 1, rel_widths = c(1/2,1/2))
  
  pdf(sprintf("%s/%s_rna_chromvar_vs_pseudotime_%s_knn%s.pdf",io$outdir,i,opts$trajectory_name,opts$knn), width=5, height=6)
  print(p1)
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
  .[,chromvar_sign:=c("down","up")[as.numeric(chromvar_diff>0)+1]] %>%
  .[,c("gene","chromvar_diff","chromvar_sign")]

directionality_rna.dt <- chromvar_rna_dt %>%
  .[,.(expr=mean(expr)),by=c("gene","celltype.predicted")] %>% 
  dcast(gene~celltype.predicted, value.var=c("expr")) %>%
  .[,c("gene",celltype.start,celltype.end), with=F] %>%
  setnames(c("gene","start","end")) %>%
  .[,rna_diff:=end-start] %>% 
  .[,rna_sign:=c("down","up")[as.numeric(rna_diff>0)+1]] %>%
  .[,c("gene","rna_diff","rna_sign")]
