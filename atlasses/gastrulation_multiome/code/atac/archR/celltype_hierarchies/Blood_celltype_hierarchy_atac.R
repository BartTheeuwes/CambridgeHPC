
#####################
## Define settings ##
#####################

# Load default settings
if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
  source("/Users/ricard/gastrulation_multiome_10x/utils.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
  source("/homes/ricard/gastrulation_multiome_10x/utils.R")
} else {
  stop("Computer not recognised")
}

# Options
opts$motif_annotation <- "Motif_cisbp"

opts$celltypes <- c(
  "Cardiomyocytes",
  "Allantois",
  # "ExE_mesoderm",
  # "Mesenchyme",
  "Haematoendothelial_progenitors",
  "Endothelium",
  "Blood_progenitors_1",
  "Blood_progenitors_2",
  "Erythroid1",
  "Erythroid2",
  "Erythroid3"
)

# I/O
io$outdir <- paste0(io$basedir,"/results/atac/archR/celltype_hierarchies")

#######################################
## Load pseudobulk RNA and ATAC data ##
#######################################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/rna_atac/load_rna_atac_pseudobulk.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/rna_atac/load_rna_atac_pseudobulk.R")
} else {
  stop("Computer not recognised")
}

atac.chromvar.se <- atac.chromvar.se[,opts$celltypes]
atac.peakMatrix.se <- atac.peakMatrix.se[,opts$celltypes]

#########################
## Define marker peaks ##
#########################

# atac.diff <- fread(paste0(io$archR.peak.differential.dir,"/PeakMatrix_Endothelium_vs_Erythroid3.txt.gz")) %>%
atac.diff <- fread(paste0(io$archR.peak.differential.dir,"/PeakMatrix_Blood_progenitors_1_vs_Endothelium.txt.gz")) %>%
  .[,sig:=FDR<0.01 & abs(MeanDiff)>0.15] %>% 
  .[,sign:="Up in Endothelium"] %>% 
  .[MeanDiff<0,sign:=c("Up in Erythroid")]

peaks.to.plot <- atac.diff[sig==T,idx]

# Create data.table
atac.peakMatrix.dt <- assay(atac.peakMatrix.se[peaks.to.plot,]) %>% as.data.table(keep.rownames = T) %>%
  setnames("rn","idx") %>% melt(id.vars="idx", variable.name="celltype")

##########
## Plot ##
##########

# Plot pseudobulk 
to.plot <- merge(
  atac.diff[idx%in%peaks.to.plot,c("idx","sign")], 
  atac.peakMatrix.dt, 
  by = "idx", 
  allow.cartesian = TRUE
)
length(unique(to.plot$celltype))
length(unique(to.plot$idx))

# to.plot[,mean(value),by=c("celltype","sign")]

order.celltypes <- to.plot %>% 
  .[,.(value=mean(value)),by=c("celltype","sign")] %>% 
  dcast(celltype~sign) %>% 
  .[,diff:=`Up in Endothelium` - `Up in Erythroid`] %>%
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
