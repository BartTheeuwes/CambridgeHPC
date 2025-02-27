#####################
## define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
  source("/Users/ricard/gastrulation_multiome_10x/utils.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
  source("/homes/ricard/gastrulation_multiome_10x/utils.R")
} else {
  stop("Computer not recognised")
}


# I/O
io$outdir <- paste0(io$basedir,"/results/atac/archR/marker_peaks/pseudobulk")

# Options
opts$celltypes <- c(
  "Epiblast",
  "Primitive_Streak",
  "Caudal_epiblast",
  # "PGC",
  # "Anterior_Primitive_Streak",
  "Notochord",
  "Def._endoderm",
  "Gut",
  "Nascent_mesoderm",
  # "Mixed_mesoderm",
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
  # "Caudal_neurectoderm",
  "Neural_crest",
  "Forebrain_Midbrain_Hindbrain",
  "Spinal_cord",
  "Surface_ectoderm"
  # "Visceral_endoderm"
  # "ExE_endoderm",
  # "ExE_ectoderm",
  # "Parietal_endoderm"
)

opts$aggregate.celltypes <- c(
  "Erythroid1" = "Erythroid1",
  "Erythroid2" = "Erythroid1",
  "Erythroid3" = "Erythroid1",
  "Blood_progenitors_1" = "Blood_progenitors_1",
  "Blood_progenitors_2" = "Blood_progenitors_1",
  "Rostral_neurectoderm" = "Neurectoderm",
  "Caudal_neurectoderm" = "Neurectoderm",
  "Anterior_Primitive_Streak" = "Primitive_Streak"
)

#######################
## Load marker peaks ##
#######################

marker_peaks.dt <- fread(io$archR.markers_peaks) %>%
  .[celltype%in%opts$celltypes] %>%
  .[,celltype:=stringr::str_replace_all(celltype,opts$aggregate.celltypes)] %>%
  .[,.(score=mean(score),N=mean(N)), by=c("celltype","idx")]

##########################
## Load pseudobulk ATAC ##
##########################

# Load SummarizedExperiment
atac.peakMatrix.se <- readRDS(io$archR.pseudobulk.peakMatrix.se)[,opts$celltypes]

# Merge cell types

# Load peak metadata
peak_metadata.dt <- fread(io$archR.peak.metadata) %>% 
  .[,idx:=sprintf("%s:%s-%s",chr,start,end)]

# Define peak names
peak_names <- rowData(atac.peakMatrix.se) %>% as.data.table %>% .[,idx:=sprintf("%s:%s-%s",seqnames,start,end)] %>% .$id
rownames(atac.peakMatrix.se) <- peak_names

# Subset marker peaks
peaks.to.plot <- unique(marker_peaks.dt$idx)
atac.peakMatrix.se <- atac.peakMatrix.se[peaks.to.plot,]
peak_metadata.dt <-peak_metadata.dt[idx%in%peaks.to.plot]

# Convert to data.table
atac.dt <- assay(atac.peakMatrix.se) %>% t %>%
  as.data.table(keep.rownames = T) %>%
  setnames("rn","celltype") %>%
  melt(id.vars=c("celltype"), variable.name="idx", value.name="acc") %>%
  .[,celltype:=stringr::str_replace_all(celltype,opts$aggregate.celltypes)] %>%
  .[,.(acc=mean(acc)), by=c("celltype","idx")]

###############################
## Scatterplot per cell type ##
###############################

# opts$min.atac <- 0
# opts$max.atac <- 3
# atac.dt[atac>opts$max.atac,atac:=opts$max.atac]

to.plot <- atac.dt %>%
  merge(marker_peaks.dt[,c("celltype","idx")] %>% setnames("celltype","class") %>% .[,class:=sprintf("%s markers",class)], by=c("idx"), allow.cartesian=T) %>%
  .[,.(acc=mean(acc)),by=c("class","celltype")] %>%
  .[,celltype:=factor(celltype,levels=opts$celltypes)]
  
p <- ggbarplot(to.plot, x="celltype", y="acc", fill="celltype", stat="identity") +
  facet_wrap(~class) +
  scale_fill_manual(values=opts$celltype.colors) +
  labs(x="", y="Chromatin accessibility") +
  theme(
    axis.text = element_text(color="black"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  )
  
pdf(sprintf("%s/barplot_marker_peaks_acc_pseudobulk.pdf",io$outdir), width = 12, height = 10)
print(p)
dev.off()


################################
## Exploration of hierarchies ##
################################

# Load precomputed PCA
pca.mtx <- readRDS("/Users/ricard/data/gastrulation_multiome_10x/results/rna/pseudobulk/rna_pca_pseudobulk.rds")

# Define celltype distances
celltype_distance.dt <- dist(pca.mtx) %>% as.matrix %>% as.data.table(keep.rownames = T) %>%
  setnames("rn","celltypeA") %>%
  melt(id.vars="celltypeA", variable.name="celltypeB", value.name="distance")

# Plot all celltypes
to.plot2 <- to.plot %>%
  # .[class%in%c("Epiblast markers")] %>% 
  merge(
    celltype_distance.dt[celltypeA=="Epiblast"] %>% .[,celltypeA:=NULL] %>% setnames("celltypeB","celltype"), by="celltype"
  ) %>% melt(id.vars=c("celltype","class","distance"), variable.name="modality")

p <- ggscatter(to.plot2, x="distance", y="value", fill="celltype", size=4, shape=21,
          add="reg.line", add.params = list(color="blue", fill="lightgray"), conf.int=TRUE) +
  stat_cor(method = "pearson") +
  facet_wrap(~class, scales="fixed") +
  scale_fill_manual(values=opts$celltype.colors) +
  labs(x="Distance between celltypes", y="Chromatin accessibility") +
  theme(
    strip.text = element_text(size=rel(0.85)),
    axis.text = element_text(color="black", size=rel(0.60)),
    axis.title = element_text(color="black", size=rel(0.80)),
    legend.position = "none"
  )


pdf(sprintf("%s/marker_peaks_acc_vs_celltype_distance_pseudobulk.pdf",io$outdir), width = 12, height = 10)
print(p)
dev.off()

# Plot epiblast
to.plot2 <- to.plot[class=="Epiblast markers"] %>%
  merge(
    celltype_distance.dt[celltypeA=="Epiblast"] %>% .[,celltypeA:=NULL] %>% setnames("celltypeB","celltype"), by="celltype"
  ) %>% melt(id.vars=c("celltype","class","distance"), variable.name="modality")

p <- ggscatter(to.plot2, x="distance", y="value", fill="celltype", size=4, shape=21,
               add="reg.line", add.params = list(color="black", fill="lightgray"), conf.int=TRUE) +
  stat_cor(method = "pearson") +
  facet_wrap(~class, scales="fixed") +
  scale_fill_manual(values=opts$celltype.colors) +
  labs(x="Distance between celltypes", y="Chromatin accessibility") +
  theme(
    strip.text = element_text(size=rel(0.85)),
    axis.text = element_text(color="black", size=rel(0.60)),
    axis.title = element_text(color="black", size=rel(0.80)),
    legend.position = "none"
  )


pdf(sprintf("%s/marker_peaks_acc_vs_celltype_distance_pseudobulk_Epiblast.pdf",io$outdir), width = 5, height = 5)
print(p)
dev.off()