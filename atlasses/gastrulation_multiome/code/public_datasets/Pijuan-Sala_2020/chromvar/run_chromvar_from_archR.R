
########################
## Load ArchR project ##
########################

source("/Users/ricard/gastrulation_multiome_10x/public_datasets/Pijuan-Sala_2020/archR/load_archR_project.R")

#####################
## Define settings ##
#####################

io$outdir <- paste0(io$basedir,"/results/chromVAR/archR")

#####################
## Update metadata ##
#####################

sample_metadata <-sample_metadata %>% 
  .[,celltype:=factor(celltype,levels=names(opts$celltype.colors))]

##########################
## Add motif annotation ##
##########################

opts$motifSet <- "cisbp"     # [JASPAR2016, JASPAR2018, JASPAR2020, cisbp, encode, homer]
opts$collection <- "CORE"    # only for JASPAR motif sets.
opts$motif.pvalue.cutoff <- 5e-05  # default is 5e-05

# add motif set 
ArchRProject <- addMotifAnnotations(
  ArchRProject, 
  motifSet = opts$motifSet,      
  collection = opts$collection,  
  cutOff = opts$motif.pvalue.cutoff,   
  name = "Motif"
)
names(getPeakAnnotation(ArchRProject))

##########################
## Add background peaks ##
##########################

# Background peaks are chosen by sampling peaks based on similarity in GC content and # of fragments across samples using the Mahalanobis distance. 
# The w paramter controls how similar background peaks should be. The bs parameter controls the precision with which the similarity is computed; 
# increasing bs will make the function run slower.
# Returns a matrix with one row per peak and one column per iteration. values in a row represent indices of background peaks for the peak with that index

# see chromVAR::getBackgroundPeaks()
ArchRProject <- addBgdPeaks(
  ArchRProject,
  nIterations = 50,
  w = 0.1,
  binSize = 50,
  method = "chromVAR",
  seed = 42
)
getBgdPeaks(ArchRProject, method = "chromVAR")

###################################
## Compute deviations for motifs ##
###################################

# The function computeDeviations returns a SummarizedExperiment with two "assays":  
# - The first matrix (accessible via `deviations(dev)` or `assays(dev)$deviations)` will give the bias corrected deviation in accessibility for each set of peaks (rows) for each cell or sample (columns). This metric represent how accessible the set of peaks is relative to the expectation based on equal chromatin accessibility profiles across cells/samples, normalized by a set of background peak sets matched for GC and average accessibility. 
# - The second matrix `deviationScores(dev)` or `assays(deviations)$z` gives the deviation Z-score, which takes into account how likely such a score would occur if randomly sampling sets of beaks with similar GC content and average accessibility.

ArchRProject <- addDeviationsMatrix(ArchRProject, 
  bgdPeaks = getBgdPeaks(ArchRProject, method = "chromVAR"),
  matrixName = "DeviationMatrix",
  peakAnnotation = "Motif",
  out = c("z", "deviations"),
  binarize = FALSE,
  force = TRUE
)

deviations.se <- getMatrixFromProject(ArchRProject, "DeviationMatrix")

#############################################
## Boxplots of motif z-score per cell type ##
#############################################

markerMotifs <- getFeatures(ArchRProject, useMatrix = "DeviationMatrix")
markerMotifs <- grep("z:", markerMotifs, value = TRUE)
# markerMotifs <- getFeatures(ArchRProject, select = paste(motifs, collapse="|"), useMatrix = "DeviationMatrix")

z.dt <- assay(deviations.se,"z") %>% as.data.frame %>%
  as.data.table(keep.rownames = T) %>% setnames("rn","motif") %>%
  .[,motif:=stringr::str_split(motif,"_") %>% map_chr(1) %>% gsub("z:","",.)] %>%
  melt(id.vars="motif", variable.name="cell") %>%
  .[,cell:=gsub("E8.5_gastrulation#","",cell)] %>%
  merge(sample_metadata[,c("cell","celltype","umap_X","umap_Y")])# %>%
  # .[,.(value=mean(value)),by=c("celltype","motif")]

for (i in unique(z.dt$motif)) {
  to.plot <- z.dt[motif==i]
  
  p <- ggboxplot(to.plot, x="celltype", y="value", fill="celltype") +
    scale_fill_manual(values=opts$celltype.colors) +
    geom_hline(yintercept=0, linetype="dashed") +
    labs(x="", y=sprintf("%s z-score",i)) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(color="black", angle=40, hjust=1, size=rel(0.75)),
      axis.text.y = element_text(color="black", size=rel(0.8))
    )
  pdf(sprintf("%s/pdf/chromvar_boxplot_%s.pdf",io$outdir,i), width=7, height=4)
  print(p)
  dev.off()
}

###########################################
## UMAP, cells coloured by motif z-score ##
###########################################

for (i in unique(z.dt$motif)) {
  to.plot <- z.dt[motif==i]
  
  p <- ggplot(to.plot, aes(x=umap_X, y=umap_Y)) +
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
  pdf(sprintf("%s/pdf/chromvar_umap_%s.pdf",io$outdir,i), width=6, height=5)
  print(p)
  dev.off()
}

##########
## Test ##
##########

# plotVarDev <- getVarDeviations(ArchRProject, name = "MotifMatrix", plot = TRUE)
# plotVarDev

##########
## Save ##
##########

# Save summarised experiment object
saveRDS(deviations.se, sprintf("%s/deviations_summarized_experiment.rds",io$outdir))

# Save deviations matrix and z-score matrix
# for (i in c("z","deviations")) {
#   saveRDS(assay(deviations.se,i), sprintf("%s/%s.rds",io$outdir,i))
# }

# VarDev <- getDev(ArchRProject, name = "MotifMatrix", plot=F)
