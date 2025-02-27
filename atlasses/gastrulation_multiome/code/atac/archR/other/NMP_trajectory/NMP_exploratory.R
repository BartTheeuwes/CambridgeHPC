library(ggpubr)

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
io$outdir <- paste0(io$basedir,"/results/atac/archR/NMP_trajectory")

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
  .[pass_atacQC==TRUE] %>%
  .[sample%in%opts$samples & celltype.predicted%in%opts$celltypes]
stopifnot(sample_metadata$archR_cell %in% rownames(ArchRProject))

##################
## Subset ArchR ##
##################

ArchRProject.filt <- ArchRProject[sample_metadata$archR_cell,]

##############################
## Dimensionality reduction ##
##############################

ArchRProject.filt <- addIterativeLSI(
  ArchRProject.filt,
  dimsToUse = 1:5,
  useMatrix = "PeakMatrix", 
  name = "IterativeLSI", 
  depthCol = "nFrags",
  varFeatures = 3000,
  force = TRUE
)


ArchRProject.filt <- addUMAP(
  ArchRProject.filt, 
  reducedDims = "IterativeLSI",
  dimsToUse = 2:5,
  name = "UMAP",
  metric = "cosine",
  nNeighbors = 15, 
  minDist = 0.15, 
  seed = 42,
  force = TRUE
)


# to.plot <- getReducedDims(ArchRProject.filt, "IterativeLSI", dimsToUse=c(2,3)) %>%
to.plot <- getEmbedding(ArchRProject.filt,"UMAP") %>%
  as.data.table(keep.rownames = T) %>%
  setnames(c("archR_cell","umap1","umap2")) %>% merge(sample_metadata,by="archR_cell")

p <- ggplot(to.plot, aes(x=umap1, y=umap2)) +
  geom_point(aes(fill=celltype.predicted), size=2, shape=21, color="black", stroke=0.05) +
  scale_fill_manual(values=opts$celltype.colors) +
  guides(fill = guide_legend(override.aes = list(size=4))) +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank()
  )

# pdf(sprintf("%s/archr_umap_celltype_%s_LSIiter%s_nfeatures%s_neighb%s_mindist%s.pdf",io$outdir,opts$matrix,opts$lsi.iterations, opts$lsi.varFeatures, opts$umap.neighbours,opts$umap.minDist))
print(p)
# dev.off()


########################################
## Overlay ATAC metrics onto RNA UMAP ##
########################################

io$nmp.pseudotime <- "/Users/ricard/data/gastrulation_multiome_10x/results/rna/NMP_trajectory/NMP_trajectory.txt.gz"
nmp.trajectory <- fread(io$nmp.pseudotime) %>%
  merge(sample_metadata[,c("cell","archR_cell")]) %>% 
  .[,cell:=NULL] %>% setnames("archR_cell","cell")

####################################
## Overlay ATAC chromVAR metrics  ##
####################################

# Load Deviation Matrix as a SummarizedExperiment object
atac.deviation.mtx <- readRDS(io$archR.deviations.se) %>% assay(.,"z")

# Subset cells
atac.deviation.mtx <- atac.deviation.mtx[,colnames(atac.deviation.mtx) %in% sample_metadata$archR_cell]
dim(atac.deviation.mtx)

# Filter motifs with no variability
atac.deviation.mtx <- atac.deviation.mtx[apply(atac.deviation.mtx,1,var)>5,]
dim(atac.deviation.mtx)

chromvar_dt <- atac.deviation.mtx %>% t %>% as.data.table(keep.rownames = T) %>%
  setnames("rn","cell") %>%
  merge(nmp.trajectory,by="cell") %>%
  melt(id.vars=c("cell","V1","V2"), variable.name="motif")

motifs.to.plot <- unique(chromvar_dt$motif) %>% as.character# %>% head

# ADD MOTIF


for (i in motifs.to.plot) {
  
  to.plot <- chromvar_dt[motif==i] %>% 
    setnames("cell","archR_cell") %>%
    merge(sample_metadata,by="archR_cell") %>%
    setnames("celltype.predicted","celltype") %>%
    .[,celltype:=factor(celltype,levels=opts$celltypes)]
  
  p1 <- ggplot(to.plot, aes(x=V1, y=value)) +
    geom_point(aes(fill=celltype), size=1.5, shape=21, stroke=0.1) +
    stat_smooth(method="loess", color="black", alpha=0.75, span=0.5) +
    geom_rug(aes(color=celltype), sides="b") +
    scale_color_manual(values=opts$celltype.colors) +
    scale_fill_manual(values=opts$celltype.colors) +
    # scale_fill_brewer(palette="Dark2") +
    labs(x="Pseudotime", y=sprintf("%s accessibility score",i)) +
    theme_classic() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position="none"
    )
  
  my_comparisons <- list( c("Somitic_mesoderm", "NMP"), c("NMP", "Spinal_cord") )
  
  p2 <- ggboxplot(to.plot, x="celltype", y="value", fill="celltype", outlier.shape=NA) +
    scale_fill_manual(values=opts$celltype.colors) +
    labs(x="", y=sprintf("%s accessibility score",i)) +
    stat_compare_means(comparisons = my_comparisons, label="p.signif", hide.ns=T) +
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
  
  p <- cowplot::plot_grid(plotlist=list(p1,p2), nrow = 1, rel_widths = c(2/3,1/3))
  
  
  pdf(sprintf("%s/chromvar/%s_chromvar.pdf",io$outdir,i), width=10, height=4)
  print(p)
  dev.off()
}

##############################
## Overlay ATAC gene scores ##
##############################

# Load HVGs
hvgs <- readRDS("/Users/ricard/data/gastrulation_multiome_10x/results/rna/NMP_trajectory/hvgs.rds")

# ArchR::addGeneScoreMatrix()
gene.score.se <- getMatrixFromProject(ArchRProject.filt, useMatrix = "GeneScoreMatrix")
gene.score.se.filt <- gene.score.se %>%
  .[,sample_metadata$archR_cell] %>%
  .[Matrix::rowMeans(assay(.))>0.5,] 
gene.score.matrix <- assay(gene.score.se.filt) %>% as.matrix %>% t
colnames(gene.score.matrix) <- rowData(gene.score.se.filt)$name
dim(gene.score.matrix)

# Subset 
gene.score.matrix <- gene.score.matrix[,colnames(gene.score.matrix)  %in% hvgs]
dim(gene.score.matrix)

acc_dt <- gene.score.matrix %>% as.data.table(keep.rownames = T) %>%
  setnames("rn","cell") %>%
  merge(nmp.trajectory,by="cell") %>%
  melt(id.vars=c("cell","V1","V2"), variable.name="gene")

# genes.to.plot <- unique(acc_dt$gene)[grep("Hb[a|b]",unique(acc_dt$gene))] %>% as.character
genes.to.plot <- c("T","Wnt3a","Rspo3","Fgf8")
genes.to.plot <- genes.to.plot[genes.to.plot%in%unique(acc_dt$gene)]
genes.to.plot <- unique(acc_dt$gene)

for (i in genes.to.plot) {
  
  to.plot <- acc_dt[gene==i] %>% 
    setnames("cell","archR_cell") %>%
    merge(sample_metadata,by="archR_cell") %>%
    setnames("celltype.predicted","celltype") %>%
    .[,celltype:=factor(celltype,levels=opts$celltypes)]
  
  p1 <- ggplot(to.plot, aes(x=V1, y=value)) +
    geom_point(aes(fill=celltype), size=1.5, shape=21, stroke=0.1) +
    stat_smooth(method="loess", color="black", alpha=0.75, span=0.5) +
    geom_rug(aes(color=celltype), sides="b") +
    scale_color_manual(values=opts$celltype.colors) +
    scale_fill_manual(values=opts$celltype.colors) +
    # scale_fill_brewer(palette="Dark2") +
    labs(x="Pseudotime", y=sprintf("%s accessibility score",i)) +
    theme_classic() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position="none"
    )
  
  my_comparisons <- list( c("Somitic_mesoderm", "NMP"), c("NMP", "Spinal_cord") )
  
  p2 <- ggboxplot(to.plot, x="celltype", y="value", fill="celltype", outlier.shape=NA) +
    scale_fill_manual(values=opts$celltype.colors) +
    labs(x="", y=sprintf("%s accessibility score",i)) +
    stat_compare_means(comparisons = my_comparisons, label="p.signif", hide.ns=T) +
    coord_cartesian(ylim=c(0,5)) +
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
  
  p <- cowplot::plot_grid(plotlist=list(p1,p2), nrow = 1, rel_widths = c(2/3,1/3))
  
  
  pdf(sprintf("%s/%s_gene_accessibility.pdf",io$outdir,i), width=10, height=4)
  print(p)
  dev.off()
}


##########################
## Load pseudobulk ATAC ##
##########################

pseudobulk.atac.peak.se <- readRDS(io$archR.pseudobulk.peaks)[,opts$celltypes]

# Define peak names
row.ranges.dt <- rowData(pseudobulk.atac.peak.se) %>% as.data.table %>% 
  setnames("seqnames","chr") %>%
  .[,c("chr","start","end")] %>%
  .[,idx:=sprintf("%s_%s_%s",chr,start,end)]
rownames(pseudobulk.atac.peak.se) <- row.ranges.dt$idx

################################
## Differential accessibility ##
################################

# Spinal cord vs NMP
markersPeaks.spinal_vs_nmp <- getMarkerFeatures(
  ArchRProject.filt, 
  useMatrix = "PeakMatrix",
  groupBy = "celltype.predicted",
  useGroups = c("Spinal_cord"),
  bgdGroups = c("NMP")
)
markersPeaks.spinal_vs_nmp.dt <- getMarkers(markersPeaks.spinal_vs_nmp, cutOff = "FDR <= Inf & Log2FC >= -Inf") %>%
  as.data.table %>% .[,celltype:=i] %>% .[,c("group","group_name"):=NULL] %>%
  setnames("seqnames","chr") %>%
  .[,idx:=sprintf("%s:%s-%s",chr,start,end)] %>%
  .[,sig:=FDR<0.10 & abs(MeanDiff)>0.1] %>% 
  .[,sign:="Up in Spinal_cord"] %>% 
  .[Log2FC<0,sign:=c("Up in NMP")]
markersPeaks.spinal_vs_nmp.dt[,sum(sig),by="sign"]

# Somitic mesoderm vs NMP

markersPeaks.somitic_vs_nmp <- getMarkerFeatures(
  ArchRProject.filt, 
  useMatrix = "PeakMatrix",
  groupBy = "celltype.predicted",
  useGroups = c("Somitic_mesoderm"),
  bgdGroups = c("NMP")
)
markersPeaks.somitic_vs_nmp.dt <- getMarkers(markersPeaks.somitic_vs_nmp, cutOff = "FDR <= Inf & Log2FC >= -Inf") %>%
  as.data.table %>% .[,celltype:=i] %>% .[,c("group","group_name"):=NULL] %>%
  setnames("seqnames","chr") %>%
  .[,idx:=sprintf("%s:%s-%s",chr,start,end)] %>%
  .[,sig:=FDR<0.10 & abs(MeanDiff)>0.1] %>% 
  .[,sign:="Up in Somitic_mesoderm"] %>% 
  .[Log2FC<0,sign:=c("Up in NMP")]
markersPeaks.somitic_vs_nmp.dt[,sum(sig),by="sign"]

# Find peaks that are exclusively found in NMPs
nmp.peaks <- intersect(
  markersPeaks.somitic_vs_nmp.dt[sign=="Up in NMP" & sig==TRUE & abs(MeanDiff)>0.15,idx],
  markersPeaks.spinal_vs_nmp.dt[sign=="Up in NMP" & sig==TRUE & abs(MeanDiff)>0.15,idx]
)
length(nmp.peaks)
head(nmp.peaks)



######################
## Motif enrichment ##
######################

markersPeaks.nmp <- getMarkerFeatures(
  ArchRProject.filt, 
  useMatrix = "PeakMatrix",
  groupBy = "celltype.predicted",
  useGroups = c("NMP"),
  bgdGroups = c("Somitic_mesoderm","Spinal_cord")
)
foo <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & MeanDiff <= -0.25")[[1]] %>%
  as.data.table %>% .[,idx:=sprintf("%s:%s-%s",seqnames,start,end)]


# rownames(markersPeaks.somitic_vs_nmp) <- markersPeaks.somitic_vs_nmp.dt$idx

motifsUp <- peakAnnoEnrichment(
  ArchRProject.filt,
  seMarker = markersPeaks.nmp,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.01 & MeanDiff <= -0.25",
)
names(assays(motifsUp))

dt <- motifsUp@assays@data$mlog10Padj %>%
  as.data.table(keep.rownames = T) %>%
  setnames(c("motif","log10_padj")) %>%
  setorder(-log10_padj)


for (i in unique(dt$celltype)) {
  to.plot <- dt[celltype==i] %>% 
    setorder(-value) %>%
    .[,TF:=stringr::str_split(TF,"_") %>% map_chr(1)] %>%
    .[,.(value=mean(value)), by=c("celltype","TF")] %>%
    head(n=25) %>%
    .[,TF:=factor(TF,levels=rev(TF))]
  
  p <- ggplot(to.plot, aes(x=TF, y=value)) +
    geom_point(size=2) +
    geom_segment(aes_string(xend="TF"), size=0.75, yend=0) +
    labs(x="q-value (hypergeometric test)") +
    coord_flip() +
    theme_bw()
  
  pdf(sprintf("%s/archr_motif_enrichment_lineplot_%s.pdf",io$outdir,i), width=4, height=10)
  print(p)
  dev.off()
}

##############
## Epiblast ##
##############

a <- fread("/Users/ricard/data/gastrulation_multiome_10x/results/atac/archR/differential/PeakMatrix_Epiblast_vs_Spinal_cord.txt.gz") %>%
  .[,idx:=sprintf("%s:%s-%s",chr,start,end)] %>%
  .[,sig:=FDR<0.10 & abs(MeanDiff)>0.25] %>%
  .[,sign:="Up in Epiblast"] %>% 
  .[Log2FC<0,sign:=c("Up in Spinal cord")] %>%
  setorder(FDR, na.last = T)

a[,sum(sig),by="sign"]
a[sign=="Up in Epiblast"]

markersPeaks.nmp <- getMarkerFeatures(
  ArchRProject.filt, 
  useMatrix = "PeakMatrix",
  groupBy = "celltype.predicted",
  useGroups = c("NMP"),
  bgdGroups = c("Somitic_mesoderm","Spinal_cord")
)