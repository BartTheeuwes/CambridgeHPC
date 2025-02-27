
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
io$outdir <- paste0(io$basedir,"/results/atac/archR/blood")

# Options
opts$samples <- c(
  "E7.5_rep1",
  "E7.5_rep2",
  "E8.5_rep1",
  "E8.5_rep2"
)

opts$celltypes <- c(
  "Allantois",
  # "ExE_mesoderm",
  # "Mesenchyme",
  "Haematoendothelial_progenitors",
  "Blood_progenitors_1",
  "Blood_progenitors_2",
  "Erythroid1",
  "Erythroid2",
  "Erythroid3",
  "Endothelium"
)

#####################
## Update metadata ##
#####################

sample_metadata <- fread(io$metadata) %>%
  .[pass_atacQC==TRUE] %>%
  .[sample%in%opts$samples & celltype.mapped%in%opts$celltypes]
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
  dimsToUse = 1:10,
  useMatrix = "PeakMatrix", 
  name = "IterativeLSI", 
  depthCol = "nFrags",
  varFeatures = 10000,
  force = TRUE
)

# Harmony
if (opts$batch.correction) {
  
  ArchRProject.filt@cellColData$stage <- substr(ArchRProject.filt@cellColData$sample,1,4)
  
  ArchRProject.filt <- addHarmony(
    ArchRProj = ArchRProject.filt,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "stage",
    force = TRUE
  )
  
}

ArchRProject.filt <- addUMAP(
  ArchRProject.filt, 
  reducedDims = "Harmony",
  name = "UMAP",
  metric = "cosine",
  nNeighbors = 25, 
  minDist = 0.3, 
  seed = 42,
  force = TRUE
)

to.plot <- getEmbedding(ArchRProject.filt,"UMAP") %>%
  as.data.table(keep.rownames = T) %>%
  setnames(c("archR_cell","umap1","umap2")) %>% merge(sample_metadata,by="archR_cell")

p <- ggplot(to.plot, aes(x=umap1, y=umap2)) +
  geom_point(aes(fill=celltype.mapped), size=2, shape=21, color="black", stroke=0.05) +
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


p <- ggplot(to.plot, aes(x=umap1, y=umap2)) +
  # geom_point(aes(fill=log2(nFrags_atac)), size=2, shape=21, color="black", stroke=0.05) +
  geom_point(aes(fill=BlacklistRatio_atac), size=2, shape=21, color="black", stroke=0.05) +
  theme_classic() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank()
  )

# pdf(sprintf("%s/archr_umap_celltype_%s_LSIiter%s_nfeatures%s_neighb%s_mindist%s.pdf",io$outdir,opts$matrix,opts$lsi.iterations, opts$lsi.varFeatures, opts$umap.neighbours,opts$umap.minDist))
print(p)
# dev.off()

####################################################
## Overlay ATAC metrics onto RNA-based trajectory ##
####################################################

# io$blood.trajectory <- "/Users/ricard/data/gastrulation_multiome_10x/results/rna/blood_trajectory/umap_blood.txt.gz"
# rna.umap <- fread(io$blood.trajectory) %>%
#   merge(sample_metadata[,c("cell","archR_cell")]) %>% 
#   .[,cell:=NULL] %>% setnames("archR_cell","cell")

gene.score.se <- getMatrixFromProject(ArchRProject.filt, useMatrix = "GeneScoreMatrix_nodistal")[,sample_metadata$archR_cell]

gene.score.matrix <- assay(gene.score.se) %>% as.matrix %>% t
colnames(gene.score.matrix) <- rowData(gene.score.se)$name
dim(gene.score.matrix)

acc_dt <- gene.score.matrix %>% as.data.table(keep.rownames = T) %>%
  setnames("rn","cell") %>%
  merge(rna.umap,by="cell") %>%
  melt(id.vars=c("cell","UMAP1","UMAP2"), variable.name="gene")

genes.to.plot.acc <- colnames(gene.score.matrix)[grep("Hb[a|b]",colnames(gene.score.matrix))] %>% as.character

for (i in genes.to.plot.acc) {
  p <- ggplot(acc_dt[gene==i], aes(x=UMAP1, y=UMAP2)) +
    geom_point(aes(fill=value), size=2, shape=21, color="black", stroke=0.05) +
    scale_fill_gradient(low = "gray80", high = "red") +
    # guides(fill = guide_legend(override.aes = list(size=4))) +
    theme_classic() +
    theme(
      legend.title = element_blank(),
      legend.position = "right",
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank()
    )
  
  # pdf(sprintf("%s/archr_umap_celltype_%s_LSIiter%s_nfeatures%s_neighb%s_mindist%s.pdf",io$outdir,opts$matrix,opts$lsi.iterations, opts$lsi.varFeatures, opts$umap.neighbours,opts$umap.minDist))
  print(p)
  # dev.off()
}

########################################
## Boxplots of chromatin accessibility ##
########################################

tmp <- sample_metadata %>%
  .[,c("archR_cell","celltype.mapped")] %>%
  setnames(c("cell","celltype")) %>%
  .[,celltype:=factor(celltype,levels=opts$celltypes)]

# Box plots
to.plot <- acc_dt[gene%in%genes.to.plot] %>% merge(tmp,by="cell")
  
p <- ggboxplot(to.plot, x="celltype", y="value", fill="celltype", outlier.shape=NA) +
  facet_wrap(~gene, scales="fixed") +
  scale_fill_manual(values=opts$celltype.colors) +
  labs(x="", y="Gene accessibility score") +
  # guides(fill = guide_legend(override.aes = list(size=4))) +
  coord_cartesian(ylim=c(0,5.5)) +
  theme_classic() +
  guides(x = guide_axis(angle = 90)) +
  theme(
    axis.text.x = element_text(color="black", size=rel(0.75)),
    # legend.title = element_blank(),
    legend.position = "none"
    # axis.text = element_blank(),
    # axis.title = element_blank(),
    # axis.ticks = element_blank()
  )

print(p)

########################################
## Plot RNA expression of hemoglobins ##
########################################

opts$rna.cells <- fread(io$metadata) %>%
  .[pass_rnaQC==TRUE & sample%in%opts$samples & celltype.mapped%in%opts$celltypes,cell]

sce <- load_SingleCellExperiment(io$sce, normalise = TRUE, cells = opts$rna.cells)
colData(sce) <- sample_metadata %>% as.data.frame %>% tibble::column_to_rownames("cell") %>%
  .[colnames(sce),] %>% DataFrame()

genes.to.plot.rna <- rownames(sce)[grep("Hb[a|b]",rownames(sce))] %>% as.character

rna_dt <- as.matrix(logcounts(sce[genes.to.plot.rna,])) %>% t %>% as.data.table(keep.rownames="cell") %>% 
  melt(id.vars="cell", value.name="value", variable.name="gene")

to.plot <- rna_dt %>% 
  merge(sample_metadata[,c("cell","celltype.mapped")],by="cell") %>%
  setnames("celltype.mapped","celltype") %>%
  .[,celltype:=factor(celltype,levels=opts$celltypes)]

p <- ggboxplot(to.plot, x="celltype", y="value", fill="celltype", outlier.shape=NA) +
  facet_wrap(~gene, scales="fixed") +
  scale_fill_manual(values=opts$celltype.colors) +
  labs(x="", y="RNA expression") +
  # guides(fill = guide_legend(override.aes = list(size=4))) +
  # coord_cartesian(ylim=c(0,5.5)) +
  theme_classic() +
  guides(x = guide_axis(angle = 90)) +
  theme(
    axis.text.x = element_text(color="black", size=rel(0.75)),
    # legend.title = element_blank(),
    legend.position = "none"
    # axis.text = element_blank(),
    # axis.title = element_blank(),
    # axis.ticks = element_blank()
  )

print(p)


##############################################################################
## Correlation between RNA expression and chromatin accessibility for Hbb's ##
##############################################################################

# Match RNA genes with accessibiliy genes
genes.to.plot <- intersect(genes.to.plot.acc, genes.to.plot.rna)

to.plot <- merge(
  rna_dt[gene%in%genes.to.plot,c("cell","gene","value")] %>% setnames("value","rna") %>%
    merge(sample_metadata[,c("cell","archR_cell","celltype.mapped")]) %>% .[,cell:=NULL] %>% setnames("archR_cell","cell"),
  acc_dt[gene%in%genes.to.plot,c("cell","gene","value")] %>% setnames("value","acc"),
  by=c("cell","gene")
)

to.plot <- to.plot[celltype.mapped%in%c("Erythroid1","Erythroid2","Erythroid3")]

p <- ggscatter(to.plot, x="rna", y="acc", fill="celltype.mapped", shape=21, size=1) + 
               # add="reg.line", add.params = list(color="blue", fill="lightgray"), conf.int=TRUE) +
  facet_wrap(~gene, scales="fixed") +
  labs(x="RNA expression", y="Chromatin accessibility") +
  coord_cartesian(ylim=c(0,10)) +
  # geom_abline(intercept=0, slope=1) +
  scale_fill_manual(values=opts$celltype.colors) +
  theme(
    axis.text = element_text(colour="black",size=rel(0.65)),  
    legend.position = "none"
  )
p


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