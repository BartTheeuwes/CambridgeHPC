library(ggpubr)

########################
## Load ArchR project ##
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

io$metadata <- paste0(io$basedir,"/sample_metadata.txt.gz")
io$outdir <- paste0(io$basedir,"/results/atac/archR/celltype_hierarchies")
io$atac.pseudobulk.peaks <- paste0(io$archR.directory,"/pseudobulk/pseudobulk_PeakMatrix_summarized_experiment.rds")
io$atac.marker.peaks <- paste0(io$basedir, "/results/atac/archR/marker_peaks/all_peaks/markersPeaks.tsv.gz")

opts$samples <- c(
  "E7.5_rep1",
  "E7.5_rep2"
  # "E8.5_rep1",
  # "E8.5_rep2"
)

########################
## Load cell metadata ##
########################

sample_metadata <- fread(io$metadata) %>%
  .[pass_atacQC==TRUE] %>%
  .[sample%in%opts$samples & stage=="E7.5" & celltype.predicted%in%c("Epiblast","Primitive_Streak","Nascent_mesoderm")]
table(sample_metadata$celltype.predicted)

##################ç
## Subset ArchR ##
##################

ArchRProject.filt <- ArchRProject[sample_metadata$archR_cell]
ArchRProject.filt@sampleColData <- ArchRProject.filt@sampleColData[opts$samples,,drop=F]

###################################
## Load precomputed marker peaks ##
###################################

opts$target.celltypes <- c(
  "Notochord",
  "Gut",
  "Mixed_mesoderm",
  "Intermediate_mesoderm",
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
  "Neural_crest",
  "Forebrain_Midbrain_Hindbrain",
  "Spinal_cord",
  "Surface_ectoderm"
)

markersPeaks.dt <- fread(io$atac.marker.peaks)

markersPeaks.dt.filt <- markersPeaks.dt %>%
  .[celltype%in%opts$target.celltypes] %>%
  .[FDR<0.05 & MeanDiff>0.15] %>%
  .[,idx:=sprintf("%s:%s-%s",chr,start,end)]
markersPeaks.dt.filt[,.N,by="celltype"]

###############################################################
## Load peak accessibility per Epiblast cells (single cells) ##
###############################################################

atac.peak.se <- getMatrixFromProject(ArchRProject.filt, useMatrix="PeakMatrix")

peak_names <- rowRanges(atac.peak.se) %>% as.data.table %>% .[,id:=sprintf("%s:%s-%s",seqnames,start,end)] %>% .$id
rownames(atac.peak.se) <- peak_names

atac.peak.se.filt <- atac.peak.se[markersPeaks.dt.filt$idx,sample_metadata$archR_cell]
dim(atac.peak.se.filt)

#############################################################
## Load peak accessibility per Epiblast cells (pseudobulk) ##
#############################################################

pseudobulk.atac.peak.se <- readRDS(io$atac.pseudobulk.peaks)[,"Epiblast"]

row.ranges.dt <- rowData(pseudobulk.atac.peak.se) %>% as.data.table %>% 
  setnames("seqnames","chr") %>%
  .[,idx:=sprintf("%s:%s-%s",chr,start,end)]
rownames(pseudobulk.atac.peak.se) <- row.ranges.dt$idx

#####################
## Plot single cell ##
#####################

foo <- assay(atac.peak.se.filt) %>% as.matrix %>% as.data.table(keep.rownames = T) %>%
  setnames("rn","idx") %>% 
  melt(id.vars="idx", variable.name="cell") %>%
  merge(
    markersPeaks.dt.filt[,c("idx","celltype")] %>% setnames("celltype","celltype_marker"), 
    by = "idx", 
    allow.cartesian = TRUE
  ) %>%
  .[,.(value=mean(value)),by=c("cell","celltype_marker")]


mtx <- foo %>% dcast(cell~celltype_marker) %>%
  matrix.please

pca <- prcomp(mtx, rank.=5)
pca.var.explained <- 100*(pca$sdev**2 / sum(pca$sdev**2)) %>% round(4)

to.plot <- data.table(
  celltype = rownames(pca$rotation),
  PC = pca$rotation[,2]
)
ggbarplot(to.plot, x="celltype", y="PC", fill="celltype", sort.val = "asc")   +
  scale_fill_manual(values=opts$celltype.colors) +
  guides(x = guide_axis(angle = 90)) +
  theme(
    axis.text.x = element_text(color="black", size=rel(0.75)),
    # axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    # legend.title = element_blank(),
    legend.position = "none"
    # axis.title = element_blank(),
  )

to.plot <- pca$x %>% as.data.table %>% .[,c("archR_cell"):=rownames(mtx)] %>%
  merge(foo %>% dcast(cell~celltype_marker) %>% setnames("cell","archR_cell"), by="archR_cell") %>%
  merge(sample_metadata[,c("archR_cell","celltype.predicted")])

ggscatter(to.plot, x="PC2", y="PC3", color="Surface_ectoderm") +
  viridis::scale_colour_viridis()

ggscatter(to.plot, x="PC2", y="PC3", color="celltype.predicted") +
  scale_color_manual(values=opts$celltype.colors) 

for (i in unique(foo$cell)) {
  p <- ggbarplot(foo[cell==i], x="celltype_marker", y="value", fill="celltype_marker", outlier.shape=NA) +
    scale_fill_manual(values=opts$celltype.colors) +
    # guides(fill = guide_legend(override.aes = list(size=4))) +
    # coord_cartesian(ylim=c(0,0.5)) +
    guides(x = guide_axis(angle = 90)) +
    theme_classic() +
    theme(
      axis.text.x = element_text(color="black", size=rel(0.75)),
      # axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      # legend.title = element_blank(),
      legend.position = "none"
      # axis.title = element_blank(),
    )
  
  # pdf(sprintf("%s/archr_umap_celltype_%s_LSIiter%s_nfeatures%s_neighb%s_mindist%s.pdf",io$outdir,opts$matrix,opts$lsi.iterations, opts$lsi.varFeatures, opts$umap.neighbours,opts$umap.minDist))
  print(p)
  # dev.off()
}
#####################
## Plot pseudobulk ##
#####################

# Create data.table
average.peak.celltype.dt <- assay(pseudobulk.atac.peak.se) %>% as.data.table(keep.rownames = T) %>%
  setnames("rn","idx") %>% melt(id.vars="idx", variable.name="celltype")

to.plot <- merge(
  markersPeaks.dt.filt[,c("idx","MeanDiff","celltype")] %>% setnames("celltype","celltype_marker"), 
  average.peak.celltype.dt, 
  by = "idx", 
  allow.cartesian = TRUE
)

p <- ggboxplot(to.plot, x="celltype_marker", y="value", fill="celltype_marker", outlier.shape=NA) +
  scale_fill_manual(values=opts$celltype.colors) +
  # guides(fill = guide_legend(override.aes = list(size=4))) +
  # coord_cartesian(ylim=c(0,0.5)) +
  guides(x = guide_axis(angle = 90)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(color="black", size=rel(0.75)),
    # axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    # legend.title = element_blank(),
    legend.position = "none"
    # axis.title = element_blank(),
  )

# pdf(sprintf("%s/archr_umap_celltype_%s_LSIiter%s_nfeatures%s_neighb%s_mindist%s.pdf",io$outdir,opts$matrix,opts$lsi.iterations, opts$lsi.varFeatures, opts$umap.neighbours,opts$umap.minDist))
print(p)
 # dev.off()
