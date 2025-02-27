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

# io$metadata <- paste0(io$basedir,"/processed/atac/archR/sample_metadata_after_archR.txt.gz")
# io$metadata <- paste0(io$basedir,"/sample_metadata_old.txt.gz")
# io$metadata <- paste0(io$basedir,"/results/atac/archR/celltype_assignment/sample_metadata_after_archR.txt.gz")
io$outdir <- paste0(io$basedir,"/results/atac/archR/celltype_hierarchies")

opts$samples <- c(
  "E7.5_rep1",
  "E7.5_rep2"
  # "E8.5_rep1",
  # "E8.5_rep2"
)

opts$celltypes <- c(
  "Epiblast",
  "Primitive_Streak",
  # "Caudal_epiblast",
  # "PGC",
  # "Anterior_Primitive_Streak",
  # "Notochord",
  "Def._endoderm",
  # "Gut",
  "Nascent_mesoderm",
  # "Mixed_mesoderm"
  # "Intermediate_mesoderm",
  # "Caudal_Mesoderm",
  # "Paraxial_mesoderm",
  "Somitic_mesoderm",
  # "Pharyngeal_mesoderm",
  # "Cardiomyocytes",
  # "Allantois",
  # "ExE_mesoderm",
  # "Mesenchyme",
  # "Haematoendothelial_progenitors",
  # "Endothelium",
  # "Blood_progenitors_1",
  # "Blood_progenitors_2",
  "Erythroid1",
  # "Erythroid2",
  # "Erythroid3"
  # "NMP",
  # "Rostral_neurectoderm",
  # "Caudal_neurectoderm",
  # "Neural_crest",
  # "Forebrain_Midbrain_Hindbrain"
  # "Spinal_cord",
  # "Surface_ectoderm",
  "Visceral_endoderm"
  # "ExE_endoderm",
  # "ExE_ectoderm",
  # "Parietal_endoderm"
)

########################
## Load cell metadata ##
########################

sample_metadata <- fread(io$metadata) %>%
  .[pass_atacQC==TRUE & doublet_call==FALSE] %>%
  .[sample%in%opts$samples & celltype.mapped%in%opts$celltypes]

table(sample_metadata$celltype.mapped)

stopifnot(sample_metadata$cell %in% rownames(ArchRProject))

##################
## Subset ArchR ##
##################

ArchRProject.filt <- ArchRProject[sample_metadata$cell]
# ArchRProject.filt@sampleColData <- ArchRProject.filt@sampleColData[opts$samples,,drop=F]
table(getCellColData(ArchRProject.filt,"Sample")[[1]])
table(getCellColData(ArchRProject.filt,"celltype.mapped")[[1]])

###################################
## Load precomputed marker peaks ##
###################################

# io$atac.marker.peaks <- paste0(io$basedir, "/results/atac/archR/marker_peaks/all_peaks/markersPeaks.tsv.gz")
# markersPeaks.dt <- fread(io$atac.marker.peaks)
# markersPeaks.dt.filt <- markersPeaks.dt %>%
#   .[celltype%in%c("Haematoendothelial_progenitors")] %>%
#   .[,sig:=FDR<0.1 & abs(Log2FC)>1] %>% 
#   .[,sign:="Up in Haematoendothelial_progenitors"] %>% .[Log2FC<0,sign:=c("Down in Haematoendothelial_progenitors")] %>%
#   .[,idx:=sprintf("%s:%s-%s",chr,start,end)]
# 
# markersPeaks.dt.filt[,mean(sig),by="sign"]
# sum(markersPeaks.dt.filt$sig)

##################
## Define peaks ##
##################

# markersPeaks <- getMarkerFeatures(
#   ArchRProject.filt, 
#   useMatrix = "PeakMatrix",
#   groupBy = "celltype.mapped",
#   testMethod = "wilcoxon",
#   bias = c("TSSEnrichment", "log10(nFrags)"),
#   useGroups = c("Nascent_mesoderm"),
#   bgdGroups = c("Def._endoderm")
# )
# 
# # markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")
# markerList <- getMarkers(markersPeaks, cutOff = "FDR <= Inf & Log2FC >= -Inf")
# 
# markersPeaks.dt.filt <- names(markerList) %>%
#   map(function(i) as.data.table(markerList[[i]]) %>% .[,celltype:=i]) %>% 
#   rbindlist %>%
#   setnames("seqnames","chr") %>%
#   .[,c("idx"):=NULL] %>%
#   .[,idx:=sprintf("%s:%s-%s",chr,start,end)] %>%
#   .[,sig:=FDR<0.01 & abs(MeanDiff)>0.1] %>% 
#   .[,sign:="Up in Mesoderm"] %>% .[Log2FC<0,sign:=c("Up in Endoderm")]
# 
# markersPeaks.dt.filt[,sum(sig),by="sign"]

#######################
## Precomputed peaks ##
#######################

opts$celltype1 <- "Nascent_mesoderm"
opts$celltype2 <- "Primitive_Streak"

markersPeaks.dt.filt <- fread(sprintf("%s/PeakMatrix_%s_vs_%s.txt.gz",io$archR.peak.differential.dir,opts$celltype1,opts$celltype2)) %>%
  .[,sig:=FDR<0.01 & abs(MeanDiff)>0.1] %>%
  .[,sign:=sprintf("Up in %s",opts$celltype1)] %>% .[MeanDiff<0,sign:=sprintf("Up in %s",opts$celltype2)]

markersPeaks.dt.filt[,sum(sig),by="sign"]

###################################################
## Load average peak accessibility per cell type ##
###################################################

# io$archR.pseudobulk.peakMatrix.se
# io$archR.pseudobulk.peakMatrix.se <- "/Users/ricard/data/gastrulation_multiome_10x/processed/atac/archR/pseudobulk/pseudobulk_PeakMatrix_summarized_experiment.rds"
average.peak.celltype.se <- readRDS(io$archR.pseudobulk.peakMatrix.se)
rownames(average.peak.celltype.se) <- rowData(average.peak.celltype.se) %>% as.data.table %>% 
  setnames("seqnames","chr") %>% .[,idx:=sprintf("%s:%s-%s",chr,start,end)] %>% .$idx

# Subset
average.peak.celltype.se.filt <- average.peak.celltype.se %>% .[markersPeaks.dt.filt[sig==T,idx]] %>% .[,opts$celltypes]

# Create data.table
average.peak.celltype.dt <- assay(average.peak.celltype.se.filt) %>% as.data.table(keep.rownames = T) %>%
  setnames("rn","idx") %>% melt(id.vars="idx", variable.name="celltype")

##########
## Plot ##
##########

# Are all peaks that appear at the PS stage inherited by mesoderm and endodermal lineages?

to.plot <- merge(
  markersPeaks.dt.filt[,c("idx","MeanDiff","FDR","sig","sign")], 
  average.peak.celltype.dt, 
  by = "idx", 
  allow.cartesian = TRUE
)
length(unique(to.plot$celltype))
length(unique(to.plot$idx))

# to.plot[,mean(value),by=c("celltype","sign")]

order.celltypes <- to.plot %>% .[,.(value=mean(value)),by=c("celltype","sign")] %>% 
  dcast(celltype~sign) %>% 
  .[,diff:=`Up in Nascent_mesoderm` - `Up in Primitive_Streak`] %>%
  setorder(-diff) %>% .$celltype

to.plot[,celltype:=factor(celltype,levels=order.celltypes)]

p <- ggboxplot(to.plot, x="celltype", y="value", fill="sign", outlier.shape=NA) +
  # scale_fill_manual(values=opts$celltype.colors) +
  # guides(fill = guide_legend(override.aes = list(size=4))) +
  coord_cartesian(ylim=c(0,0.6)) +
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

