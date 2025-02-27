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
# io$metadata <- paste0(io$basedir,"/processed/atac/archR/sample_metadata_after_archR.txt.gz")
io$outdir <- paste0(io$basedir,"/results/atac/archR/chromvar/JASPAR")

# Options
opts$samples <- c(
  "E7.5_rep1",
  "E7.5_rep2",
  "E8.5_rep1",
  "E8.5_rep2"
)

opts$deviation.matrix <- "DeviationMatrix_JASPAR"

#####################
## Update metadata ##
#####################

sample_metadata <- fread(io$metadata) %>%
  .[pass_atacQC==TRUE] %>%
  .[sample%in%opts$samples]
stopifnot(sample_metadata$archR_cell %in% rownames(ArchRProject))

# subset celltypes with sufficient number of cells
opts$min.cells <- 100
sample_metadata <- sample_metadata %>%
  .[,N:=.N,by=c("celltype.predicted")] %>% .[N>opts$min.cells] %>% .[,N:=NULL]
table(sample_metadata$celltype.predicted)

##################
## Subset ArchR ##
##################

ArchRProject.filt <- ArchRProject[sample_metadata$archR_cell,]
table(getCellColData(ArchRProject.filt,"Sample")[[1]])
table(getCellColData(ArchRProject.filt,"celltype.predicted")[[1]])

#############################################
## Boxplots of motif z-score per cell type ##
#############################################

deviations.se <- getMatrixFromProject(ArchRProject.filt, opts$deviation.matrix)

markerMotifs <- getFeatures(ArchRProject.filt, useMatrix = opts$deviation.matrix)
markerMotifs <- grep("z:", markerMotifs, value = TRUE)
# markerMotifs <- getFeatures(ArchRProject.filt, select = paste(motifs, collapse="|"), useMatrix = opts$deviation.matrix)

z.dt <- assay(deviations.se,"z") %>% as.matrix %>% as.data.frame %>%
  as.data.table(keep.rownames = T) %>% setnames("rn","motif") %>%
  .[,motif:=stringr::str_split(motif,"_") %>% map_chr(1) %>% gsub("z:","",.)] %>%
  melt(id.vars="motif", variable.name="archR_cell") %>%
  merge(sample_metadata[,c("archR_cell","celltype.predicted")], by="archR_cell") %>%
  setnames("celltype.predicted","celltype")# %>%
  # .[,.(value=mean(value)),by=c("celltype","motif")]

opts$celltype.colors <- opts$celltype.colors[names(opts$celltype.colors)%in%unique(z.dt$celltype)]

for (i in unique(z.dt$motif)) {
  to.plot <- z.dt[motif==i] %>% .[,celltype:=factor(celltype,levels=names(opts$celltype.colors))]
  
  p <- ggboxplot(to.plot, x="celltype", y="value", fill="celltype", outlier.shape=NA) +
    scale_fill_manual(values=opts$celltype.colors) +
    geom_hline(yintercept=0, linetype="dashed") +
    labs(x="", y=sprintf("%s z-score",i)) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(color="black", angle=40, hjust=1, size=rel(0.75)),
      axis.text.y = element_text(color="black", size=rel(0.8))
    )
  pdf(sprintf("%s/pdf/boxplots/chromvar_boxplot_%s.pdf",io$outdir,i), width=7, height=4)
  print(p)
  dev.off()
}

###########################################
## UMAP, cells coloured by motif z-score ##
###########################################

io$umap <- "/Users/ricard/data/gastrulation_multiome_10x/results/atac/archR/dimensionality_reduction/archr_umap_coordinates_PeakMatrix_LSIiter1_nfeatures50000_neighb25_mindist0.3.txt"
umap.dt <- fread(io$umap)

for (i in unique(z.dt$motif)) {
  to.plot <- z.dt[motif==i] %>% merge(umap.dt, by="archR_cell")
  
  p <- ggplot(to.plot, aes(x=V1, y=V2)) +
    # geom_point(aes_string(fill=i), shape=21, color="black", stroke=0.1) +
    # scale_fill_gradient(low = "gray80", high = "purple") +
    geom_point(aes(color=value), size=0.75) +
    labs(title=i) +
    scale_color_gradient(low = "gray80", high = "purple") +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.title = element_blank(),
      legend.position = "none",
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank()
    )
  pdf(sprintf("%s/pdf/umaps/chromvar_umap_%s.pdf",io$outdir,i), width=6, height=5)
  print(p)
  dev.off()
}

##########
## Test ##
##########

# plotVarDev <- getVarDeviations(ArchRProject.filt, name = "MotifMatrix", plot = TRUE)
# plotVarDev
