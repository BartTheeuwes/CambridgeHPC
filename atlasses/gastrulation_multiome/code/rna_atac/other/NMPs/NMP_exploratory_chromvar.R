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
io$outdir <- paste0(io$basedir,"/results/rna_atac/NMPs/chromvar")

# Options
opts$samples <- c(
  # "E7.5_rep1",
  # "E7.5_rep2",
  "E8.5_rep1",
  "E8.5_rep2"
)

opts$celltypes = c(
  # "Intermediate_mesoderm",
  # "Caudal_Mesoderm",
  # "Paraxial_mesoderm",
  "Somitic_mesoderm",
  "NMP",
  # "Forebrain_Midbrain_Hindbrain",
  "Spinal_cord"
)

#####################
## Update metadata ##
#####################

sample_metadata <- fread(io$metadata) %>%
  .[pass_atacQC==TRUE & pass_rnaQC==TRUE] %>%
  .[sample%in%opts$samples & celltype.mapped%in%opts$celltypes] %>%
  .[,celltype.mapped:=factor(celltype.mapped,levels=opts$celltypes)] 
stopifnot(sample_metadata$archR_cell %in% rownames(ArchRProject))

##################
## Subset ArchR ##
##################

ArchRProject.filt <- ArchRProject[sample_metadata$archR_cell,]

###############
## Load data ##
###############

# Load RNA-based trajectory
io$nmp.pseudotime <- "/Users/ricard/data/gastrulation_multiome_10x/results/rna/NMP_trajectory/NMP_trajectory.txt.gz"
nmp.trajectory <- fread(io$nmp.pseudotime)# %>%
  # merge(sample_metadata[,c("cell","archR_cell")]) %>% 
  # .[,cell:=NULL] %>% setnames("archR_cell","cell")

# Load highly variable along the NMP trajectory
io$nmp.hvgs <- "/Users/ricard/data/gastrulation_multiome_10x/results/rna/NMP_trajectory/hvgs.rds"
hvgs <- readRDS(io$nmp.hvgs)

###############################
## Load SingleCellExperiment ##
###############################

sce <- load_SingleCellExperiment(io$sce, cells = sample_metadata$cell, normalise = TRUE, remove_non_expressed_genes = FALSE)

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
# foo <- nmp.trajectory %>% tibble::column_to_rownames("cell") %>% as.matrix
pca.rna <- irlba::prcomp_irlba(t(as.matrix(logcounts(sce[hvgs,]))), n=opts$npcs)$x
rownames(pca.rna) <- colnames(sce)
rna.matrix.smoothed <- smoother_aggregate_nearest_nb(mat=as.matrix(logcounts(sce)), D=pdist(pca.rna), k=opts$knn)
colnames(rna.matrix.smoothed) <- colnames(sce)

# ATAC
pca.atac <- irlba::prcomp_irlba(t(atac.deviation.mtx), n=opts$npcs)$x
rownames(pca.atac) <- colnames(atac.deviation.mtx)
atac.deviation.mtx.smoothed <- smoother_aggregate_nearest_nb(mat=atac.deviation.mtx, D=pdist(pca.atac), k=opts$knn)
colnames(atac.deviation.mtx.smoothed) <- colnames(atac.deviation.mtx)

# Visualise smoothing
hist(rna.matrix.smoothed[1:500,1:500])
hist(as.matrix(logcounts(sce))[1:500,1:500])
hist(as.matrix(atac.deviation.mtx[1:500,1:500]))
hist(as.matrix(atac.deviation.mtx.smoothed[1:500,1:500]))

##################################
## Create motif2gene annotation ##
##################################

motif2gene.dt <- getPeakAnnotation(ArchRProject.filt, name="Motif")$motifSummary %>%
  as.data.table(keep.rownames = T) %>% setnames("rn","motif") %>% .[,strand:=NULL] %>% setnames("name","gene") %>%
  .[motif%in%rownames(atac.deviation.mtx)] %>% setkey(motif) %>% .[rownames(atac.deviation.mtx)]

motif2gene.dt[,gene:=gsub("Tcfap","Tfap",gene)] %>%
  .[,gene:=gsub("Tcfe","Tfe",gene)] %>%
  .[,gene:=gsub("Nkx1","Nkx1-",gene)] %>%
  .[,gene:=gsub("Nkx2","Nkx2-",gene)] %>%
  .[,gene:=gsub("Nkx3","Nkx3-",gene)] %>%
  .[,gene:=gsub("Nkx6","Nkx6-",gene)] %>%
  .[,gene:=gsub("Foxf1a","Foxf1",gene)] %>%
  .[,gene:=gsub("Hmga1rs1","Hmga1-rs1",gene)] %>%
  .[,gene:=gsub("Mycl1$","Mycl",gene)] %>%
  .[,gene:=gsub("Dux$","Duxf3",gene)] %>%
  .[,gene:=gsub("Duxbl$","Duxbl1",gene)] %>%
  .[,gene:=gsub("Pit1$","Prop1",gene)] %>%
  .[,gene:=gsub("ENSMUSG00000079994","Sox1",gene)]
# 
motif2gene.dt <- motif2gene.dt[gene%in%rownames(sce)]
      
motif2gene.dt <- motif2gene.dt[,N:=length(unique(motif)),by="gene"] %>% .[N==1] %>% .[,N:=NULL]

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

##########
## Plot ##
##########

chromvar_rna_dt <- merge(
  rna_dt,
  chromvar_dt %>% 
    merge(sample_metadata[,c("cell","archR_cell","celltype.mapped")], by="archR_cell") %>%
    merge(motif2gene.dt[,c("motif","gene")],by="motif") %>% .[,archR_cell:=NULL],
  by = c("cell","gene")
) %>% merge(nmp.trajectory,by="cell")


genes.to.plot <- chromvar_rna_dt %>%
  .[,.(expr=mean(expr),chromvar_zscore=mean(chromvar_zscore)),by=c("gene","motif")] %>% 
  .[expr>0.1,gene]

genes.to.plot <- unique(chromvar_rna_dt$gene)

for (i in genes.to.plot) {
  
  to.plot <- chromvar_rna_dt[gene==i] %>% 
    setnames("celltype.mapped","celltype") %>%
    # .[chromvar_zscore>7,chromvar_zscore:=7] %>%
    melt(id.vars=c("cell","V1","V2","celltype"), measure.vars=c("chromvar_zscore","expr"), variable.name="modality")# %>%
    # .[,value_scaled:=value/max(value),by="modality"]
  # .[,value_scaled:=(value-min(value))/(max(value)-min(value)), by="modality"]
    
  p1 <- ggplot(to.plot, aes(x=V1, y=value)) +
    geom_point(aes(fill=celltype), size=1.5, shape=21, stroke=0.1) +
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
    stat_compare_means(comparisons = list( c("Somitic_mesoderm", "NMP"), c("NMP", "Spinal_cord") ), label="p.signif", hide.ns=T) +
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
  
  pdf(sprintf("%s/%s_rna_chromvar_vs_pseudotime.pdf",io$outdir,i), width=7, height=6)
  print(p)
  dev.off()
}

##########################
## Correlation analysis ##
##########################

opts$threshold_fdr <- 0.10

to.plot <- chromvar_rna_dt %>% copy %>%
  .[,c("chromvar_zscore","expr"):=list(chromvar_zscore + rnorm(n=.N,mean=0,sd=1e-5), expr + rnorm(n=.N,mean=0,sd=1e-5))] %>% # add some noise 
  .[, .(V1 = unlist(cor.test(chromvar_zscore, expr)[c("estimate", "p.value")])), by = c("gene","motif")] %>%
  .[, para := rep(c("r","p"), .N/2)] %>% 
  data.table::dcast(gene+motif ~ para, value.var = "V1") %>%
  .[,"padj_fdr" := list(p.adjust(p, method="fdr"))] %>%
  # .[,"log_padj_fdr" := list(-log10(padj_fdr))] %>%
  .[, sig := padj_fdr <= opts$threshold_fdr] %>% 
  setorder(padj_fdr, na.last = T)

head(to.plot,n=3)

# Volcano plot
negative_hits <- to.plot[sig==TRUE & r<0,gene]
positive_hits <- to.plot[sig==TRUE & r>0,gene]
all <- nrow(to.plot)

xlim <- max(abs(to.plot$r), na.rm=T)
ylim <- max(-log10(to.plot$padj_fdr+1e-100), na.rm=T)

p <- ggplot(to.plot, aes(x=r, y=-log10(padj_fdr+1e-100))) +
  # geom_hline(yintercept = -log10(opts$threshold_fdr), color="blue") +
  geom_segment(aes(x=0, xend=0, y=0, yend=ylim-1), color="orange", size=0.5) +
  ggrastr::geom_point_rast(aes(color=sig, size=sig)) +
  ggrepel::geom_text_repel(data=head(to.plot[sig==T],n=40), aes(x=r, y=-log10(padj_fdr+1e-100), label=gene), size=3,  max.overlaps=50) +
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
