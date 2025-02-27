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

# io$metadata <- paste0(io$basedir,"/sample_metadata.txt.gz")
io$outdir <- paste0(io$basedir,"/results/atac/archR/celltype_hierarchies")
io$atac.pseudobulk.peaks <- paste0(io$archR.directory,"/pseudobulk/pseudobulk_PeakMatrix_summarized_experiment.rds")
io$atac.marker.peaks <- paste0(io$basedir, "/results/atac/archR/marker_peaks/all_peaks/markersPeaks.tsv.gz")


opts$celltypes = c(
  "Caudal_Mesoderm",
  "Epiblast",
  "Primitive_Streak",
  "Caudal_epiblast",
  "Somitic_mesoderm",
  "NMP",
  "Forebrain_Midbrain_Hindbrain",
  "Spinal_cord"
)

########################
## Load cell metadata ##
########################

sample_metadata <- fread(io$metadata) %>%
  .[pass_atacQC==TRUE &] %>%
  .[celltype.predicted%in%opts$celltypes]
table(sample_metadata$celltype.predicted)

##################ç
## Subset ArchR ##
##################

ArchRProject.filt <- ArchRProject[sample_metadata$archR_cell]
ArchRProject.filt@sampleColData <- ArchRProject.filt@sampleColData[opts$samples,,drop=F]

###################################
## Load precomputed marker peaks ##
###################################

# markersPeaks.dt <- fread(io$atac.marker.peaks)
# 
# markersPeaks.dt.filt <- markersPeaks.dt %>%
#   .[celltype%in%c("Notochord")] %>%
#   .[,sig:=FDR<0.1 & abs(Log2FC)>1] %>% 
#   .[,sign:="Up in Notochord"] %>% .[Log2FC<0,sign:=c("Down in Notochord")] %>%
#   .[,idx:=sprintf("%s:%s-%s",chr,start,end)]
# 
# markersPeaks.dt.filt[,mean(sig),by="sign"]
# sum(markersPeaks.dt.filt$sig)

#########################
## Define marker peaks ##
#########################

markersPeaks <- getMarkerFeatures(
  ArchRProject.filt, 
  useMatrix = "PeakMatrix",
  groupBy = "celltype.predicted",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = c("Spinal_cord"),
  bgdGroups = c("Somitic_mesoderm")
)

# markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= Inf & Log2FC >= -Inf")

markersPeaks.dt.filt <- names(markerList) %>%
  map(function(i) as.data.table(markerList[[i]]) %>% .[,celltype:=i]) %>% 
  rbindlist %>%
  setnames("seqnames","chr") %>%
  .[,c("idx"):=NULL]
head(markersPeaks.dt.filt)

markersPeaks.dt.filt %>%
  .[,idx:=sprintf("%s:%s-%s",chr,start,end)] %>%
  .[,sig:=FDR<0.05 & abs(MeanDiff)>0.15] %>% 
  .[,sign:="Up in Spinal_cord"] %>% 
  .[Log2FC<0,sign:=c("Up in Somitic_mesoderm")]

#########################
## Pseudobulk analysis ##
#########################

pseudobulk.atac.peak.se <- readRDS(io$atac.pseudobulk.peaks)

row.ranges.dt <- rowData(pseudobulk.atac.peak.se) %>% as.data.table %>% 
  setnames("seqnames","chr") %>%
  .[,idx:=sprintf("%s:%s-%s",chr,start,end)]
rownames(pseudobulk.atac.peak.se) <- row.ranges.dt$idx

# Subset
pseudobulk.atac.peak.se.filt <- pseudobulk.atac.peak.se %>% .[markersPeaks.dt.filt[sig==T,idx]] %>% .[,opts$celltypes]
dim(pseudobulk.atac.peak.se.filt)

# Create data.table
average.peak.celltype.dt <- assay(pseudobulk.atac.peak.se.filt) %>% as.data.table(keep.rownames = T) %>%
  setnames("rn","idx") %>% melt(id.vars="idx", variable.name="celltype")

# Plot pseudobulk 
to.plot <- merge(
  markersPeaks.dt.filt[sig==T,c("idx","Log2FC","FDR","sig","sign")], 
  average.peak.celltype.dt, 
  by = "idx", 
  allow.cartesian = TRUE
)
length(unique(to.plot$celltype))
length(unique(to.plot$idx))

# to.plot[,mean(value),by=c("celltype","sign")]

order.celltypes <- to.plot %>% 
  .[,.(value=mean(value)),by=c("celltype","sign")] %>% 
  dcast(celltype~sign) %>% 
  .[,diff:=`Up in Somitic_mesoderm` - `Up in Spinal_cord`] %>%
  setorder(-diff) %>% .$celltype

to.plot[,celltype:=factor(celltype,levels=order.celltypes)]

p <- ggboxplot(to.plot, x="celltype", y="value", fill="sign", outlier.shape=NA) +
  # scale_fill_manual(values=opts$celltype.colors) +
  # guides(fill = guide_legend(override.aes = list(size=4))) +
  coord_cartesian(ylim=c(0,1)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(color="black", angle=30, hjust=1),
    # legend.title = element_blank(),
    # legend.position = "none",
    # axis.text = element_blank(),
    # axis.title = element_blank(),
    # axis.ticks = element_blank()
  )

# pdf(sprintf("%s/archr_umap_celltype_%s_LSIiter%s_nfeatures%s_neighb%s_mindist%s.pdf",io$outdir,opts$matrix,opts$lsi.iterations, opts$lsi.varFeatures, opts$umap.neighbours,opts$umap.minDist))
print(p)
 # dev.off()

#########################
## Single cell analysis ##
#########################

markersPeaks.dt.filt[sig==T,idx]

atac.peak.se <- getMatrixFromProject(ArchRProject.filt, useMatrix="PeakMatrix")

peak_names <- rowRanges(atac.peak.se) %>% as.data.table %>% .[,id:=sprintf("%s:%s-%s",seqnames,start,end)] %>% .$id
rownames(atac.peak.se) <- peak_names

atac.peak.se.filt <- atac.peak.se[markersPeaks.dt.filt[sig==T,idx],]
dim(atac.peak.se.filt)


foo <- assay(atac.peak.se.filt) %>% as.matrix %>% as.data.table(keep.rownames = T) %>%
  setnames("rn","idx") %>% 
  melt(id.vars="idx", variable.name="cell") %>%
  merge(
    markersPeaks.dt.filt[,c("idx","sign")] %>% setnames("sign","celltype"), 
    by = "idx", 
    allow.cartesian = TRUE
  ) %>%
  .[,.(value=mean(value)),by=c("cell","celltype")]

# PCA
to.plot <- foo %>% dcast(cell~celltype) %>%
  merge(sample_metadata[,c("archR_cell","celltype.predicted")] %>% setnames(c("cell","celltype")))

ggscatter(to.plot, x="Up in Somitic_mesoderm", y="Up in Spinal_cord", color="celltype", size=0.45)   +
  labs(x="Peaks up in Somitic", y="Peaks up in Spinal") +
  scale_color_manual(values=opts$celltype.colors) +
  guides(color = guide_legend(override.aes = list(size=4))) +
  theme(
    axis.text = element_text(color="black", size=rel(0.75)),
    legend.position = "right"
  )


