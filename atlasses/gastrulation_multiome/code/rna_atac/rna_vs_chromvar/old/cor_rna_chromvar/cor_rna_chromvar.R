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
io$outdir <- paste0(io$basedir,"/results/rna_atac/rna_vs_chromvar/all_cells")

# Options
opts$samples <- c(
  "E7.5_rep1",
  "E7.5_rep2",
  "E8.0_rep1",
  "E8.0_rep2",
  "E8.5_rep1",
  "E8.5_rep2"
)

opts$celltypes = c(
  "Epiblast",
  "Primitive_Streak",
  "Caudal_epiblast",
  "PGC",
  "Anterior_Primitive_Streak",
  "Notochord",
  "Def._endoderm",
  "Gut",
  "Nascent_mesoderm",
  "Mixed_mesoderm",
  "Intermediate_mesoderm",
  "Caudal_Mesoderm",
  "Paraxial_mesoderm",
  "Somitic_mesoderm",
  "Pharyngeal_mesoderm",
  "Cardiomyocytes",
  "Allantois",
  "ExE_mesoderm",
  "Mesenchyme",
  "Haematoendothelial_progenitors",
  "Endothelium",
  "Blood_progenitors_1",
  "Blood_progenitors_2",
  "Erythroid1",
  "Erythroid2",
  "Erythroid3",
  "NMP",
  "Rostral_neurectoderm",
  "Caudal_neurectoderm",
  "Neural_crest",
  "Forebrain_Midbrain_Hindbrain",
  "Spinal_cord",
  "Surface_ectoderm",
  "Visceral_endoderm",
  "ExE_endoderm",
  "ExE_ectoderm",
  "Parietal_endoderm"
)

# opts$motif_annotation <- "Motif_JASPAR2020_human"
opts$motif_annotation <- "Motif_cisbp"

###################
## Load metadata ##
###################

sample_metadata <- fread(io$metadata) %>%
  .[pass_atacQC==TRUE & pass_rnaQC==TRUE & doublet_call==FALSE] %>%
  .[sample%in%opts$samples & celltype.predicted%in%opts$celltypes] %>%
  .[,celltype.predicted:=factor(celltype.predicted,levels=opts$celltypes)] 

##################
## Subset ArchR ##
##################

stopifnot(sample_metadata$cell %in% rownames(ArchRProject))
ArchRProject.filt <- ArchRProject[sample_metadata$cell,]

###############################
## Load SingleCellExperiment ##
###############################

sce <- load_SingleCellExperiment(
  file = io$rna.sce, 
  cells = sample_metadata$cell, 
  normalise = TRUE, 
  remove_non_expressed_genes = FALSE
)

rownames(sce) <- toupper(rownames(sce))

##########################
## Load chromVAR scores ##
##########################

opts$motif.annotation <- "Motif_cisbp"
chromvar.se <- readRDS(sprintf("%s/deviations_summarized_experiment_%s.rds",io$archr.chromvar.dir,opts$motif.annotation))

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
  .[gene%in%rownames(chromvar.se) & gene%in%rownames(sce)] %>%
  .[,N:=length(unique(motif)),by="gene"] %>% .[N==1] %>% .[,N:=NULL]

sce <- sce[motif2gene.dt$gene,]
chromvar.se <- chromvar.se[motif2gene.dt$gene,]

#################
## Smooth data ##
#################

opts$knn <- 25
opts$npcs <- 30

## RNA ##

# PCA
io$pca.rna <- paste0(io$basedir,"/results/rna/dimensionality_reduction/all_cells/E7.5_rep1-E7.5_rep2-E8.0_rep1-E8.0_rep2-E8.5_rep1-E8.5_rep2_pca_features2500_pcs30_batchcorrectionbysample.txt.gz")
pca.rna <- fread(io$pca.rna) %>% matrix.please %>% .[sample_metadata$cell,]
# pca.rna <- irlba::prcomp_irlba(t(as.matrix(logcounts(sce[hvgs,]))), n=opts$npcs)$x
# rownames(pca.rna) <- colnames(sce)

# kNN denoising
rna.matrix.smoothed <- smoother_aggregate_nearest_nb(mat=as.matrix(logcounts(sce)), D=pdist(pca.rna), k=opts$knn)
colnames(rna.matrix.smoothed) <- colnames(sce)

## chromVAR ##

# PCA
io$pca.atac <- paste0(io$basedir,"/results/atac/archR/dimensionality_reduction/PeakMatrix/all_cells/E7.5_rep1-E7.5_rep2-E8.0_rep1-E8.0_rep2-E8.5_rep1-E8.5_rep2_lsi_features50000_ndims30.txt.gz")
pca.atac <- fread(io$pca.atac) %>% matrix.please %>% .[sample_metadata$cell,]
# pca.atac <- irlba::prcomp_irlba(t(chromvar.deviation.mtx), n=opts$npcs)$x
# rownames(pca.atac) <- colnames(chromvar.deviation.mtx)

# KNN denoising
chromvar.deviation.mtx <- smoother_aggregate_nearest_nb(mat=as.matrix(assay(chromvar.se,"z")), D=pdist(pca.atac), k=opts$knn)
colnames(chromvar.deviation.mtx) <- colnames(chromvar.se)

################
## parse data ##
################

# chromvar_dt <- as.matrix(assay(chromvar.se)) %>%
chromvar_dt <- chromvar.deviation.mtx %>%
  as.data.table(keep.rownames = F) %>%
  .[,gene:=rownames(chromvar.deviation.mtx)] %>%
  melt(id.vars=c("gene"), variable.name="cell", value.name="chromvar_zscore")

# rna_dt <- sce %>% logcounts %>% as.matrix %>%
rna_dt <- rna.matrix.smoothed %>%
  as.data.table(keep.rownames = F) %>%
  .[,gene:=rownames(rna.matrix.smoothed)] %>%
  melt(id.vars="gene", variable.name="cell", value.name="expr")

##########
## Plot ##
##########

chromvar_rna_dt <- merge(
  rna_dt,
  chromvar_dt %>% 
    merge(sample_metadata[,c("cell","celltype.predicted")], by="cell"),
    # merge(motif2gene.dt[,c("motif","gene")],by="motif"),
  by = c("cell","gene")
)

##########################
## Correlation analysis ##
##########################

opts$threshold_fdr <- 0.10

cor.dt <- chromvar_rna_dt %>% copy %>%
  .[,c("chromvar_zscore","expr"):=list(chromvar_zscore + rnorm(n=.N,mean=0,sd=1e-5), expr + rnorm(n=.N,mean=0,sd=1e-5))] %>% # add some noise 
  .[, .(V1 = unlist(cor.test(chromvar_zscore, expr)[c("estimate", "p.value")])), by = c("gene")] %>%
  .[, para := rep(c("r","p"), .N/2)] %>% 
  data.table::dcast(gene ~ para, value.var = "V1") %>%
  .[,"padj_fdr" := list(p.adjust(p, method="fdr"))] %>%
  # .[,"log_padj_fdr" := list(-log10(padj_fdr))] %>%
  .[, sig := padj_fdr <= opts$threshold_fdr] %>% 
  setorder(padj_fdr, na.last = T)

head(cor.dt,n=3)

# Save
fwrite(cor.dt, paste0(io$outdir,"/cor_rna_vs_chromvar_allcells_smoothed.txt.gz"), sep="\t", quote=F)

##################
## Volcano plot ##
##################

to.plot <- cor.dt

negative_hits <- to.plot[sig==TRUE & r<0,gene]
positive_hits <- to.plot[sig==TRUE & r>0,gene]
all <- nrow(to.plot)

xlim <- max(abs(to.plot$r), na.rm=T)
ylim <- max(-log10(to.plot$padj_fdr+1e-100), na.rm=T)

p <- ggplot(to.plot, aes(x=r, y=-log10(padj_fdr+1e-100))) +
  geom_segment(aes(x=0, xend=0, y=0, yend=ylim-1), color="orange", size=0.5) +
  geom_point(aes(color=sig, size=sig)) +
  # ggrastr::geom_point_rast(aes(color=sig, size=sig)) +
  ggrepel::geom_text_repel(data=head(to.plot[sig==T],n=50), aes(x=r, y=-log10(padj_fdr+1e-100), label=gene), size=3,  max.overlaps=100) +
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


#####################################
## Scatter plot of individual hits ##
#####################################

genes.to.plot <- unique(chromvar_rna_dt$gene)
# genes.to.plot <- cor.dt[sig==T & abs(r)>0.25,gene]

for (i in genes.to.plot) {
  
  outfile <- sprintf("%s/individual_genes/%s_rna_vs_chromvar_knn%s_allcells.png",io$outdir,i,opts$knn)
  
  if (file.exists(outfile)) {
    print(sprintf("file or %s already exists...",i))
  } else {
    
    to.plot <- chromvar_rna_dt[gene==i]
    # p1 <- ggplot(to.plot, aes(x=expr, y=chromvar_zscore)) +
    #   geom_point(aes(fill=celltype.predicted), size=1.5, shape=21, stroke=0.1) +
    #   stat_smooth(method="lm", color="black", alpha=0.75, span=0.5) +
    #   scale_fill_manual(values=opts$celltype.colors) +
    #   guides(fill=F) +
    #   labs(x=sprintf("%s RNA expression",i), y=sprintf("%s Motif accessibility (z-score)",i)) +
    #   theme_classic()
    
    p1 <- ggscatter(to.plot, x="chromvar_zscore", y="expr", color="celltype.predicted", size=1, 
                    add="reg.line", add.params = list(color="black", fill="lightgray"), conf.int=TRUE) +
      stat_cor(method = "pearson", label.x.npc = "middle", label.y.npc = "bottom") +
      scale_colour_manual(values=opts$celltype.colors) +
      labs(y=sprintf("%s expression",i), x=sprintf("%s Motif accessibility (z-score)",i)) +
      guides(color=F) +
      theme(
        axis.text = element_text(size=rel(0.7))
      )
    
    
    to.plot2 <- to.plot %>% melt(id.vars=c("cell","celltype.predicted","gene"))
    p2 <- ggboxplot(to.plot2, x="celltype.predicted", y="value", fill="celltype.predicted", outlier.shape=NA) +
      facet_wrap(~variable, nrow=2, scales="free_y") +
      scale_fill_manual(values=opts$celltype.colors) +
      geom_hline(yintercept=0, linetype="dashed") +
      labs(x="", y="") +
      theme_classic() +
      guides(x = guide_axis(angle = 90)) +
      theme(
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()
      )
    
    
    p <- cowplot::plot_grid(plotlist=list(p1,p2), nrow = 1, rel_widths = c(1/2,1/2))
    
    # pdf(sprintf("%s/individual_genes/%s_rna_vs_chromvar.pdf",io$outdir,i), width=10, height=4)
    png(outfile, width = 1000, height = 500)
    print(p)
    dev.off()
  }
}
