# https://www.ArchRProject.com/bookdown/identifying-marker-peaks-with-archr.html

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
io$metadata <- paste0(io$basedir,"/processed/atac/archR/sample_metadata_after_archR.txt.gz")
io$outdir <- paste0(io$basedir,"/results/atac/archR/marker_peaks")

# Options
opts$samples <- c(
  "E7.5_rep1",
  "E7.5_rep2",
  "E8.5_rep1",
  "E8.5_rep2"
)

opts$celltypes <- c(
  # "Epiblast",
  # "Primitive_Streak",
  # "Caudal_epiblast",
  # "PGC",
  # "Anterior_Primitive_Streak",
  # "Notochord",
  # "Def._endoderm",
  # "Gut",
  # "Nascent_mesoderm",
  # "Mixed_mesoderm",
  # "Intermediate_mesoderm",
  # "Caudal_Mesoderm",
  # "Paraxial_mesoderm",
  # "Somitic_mesoderm",
  # "Pharyngeal_mesoderm",
  # "Cardiomyocytes",
  # "Allantois",
  # "ExE_mesoderm",
  # "Mesenchyme",
  # "Haematoendothelial_progenitors",
  # "Endothelium",
  # "Blood_progenitors_1",
  # "Blood_progenitors_2",
  # "Erythroid1",
  # "Erythroid2",
  "Erythroid3",
  "NMP"
  # "Rostral_neurectoderm",
  # "Caudal_neurectoderm",
  # "Neural_crest",
  # "Forebrain_Midbrain_Hindbrain",
  # "Spinal_cord",
  # "Surface_ectoderm",
  # "Visceral_endoderm",
  # "ExE_endoderm",
  # "ExE_ectoderm",
  # "Parietal_endoderm"
)

########################
## Load cell metadata ##
########################

sample_metadata <- fread(io$metadata) %>%
  .[pass_atacQC==TRUE] %>%
  .[,archR_cell:=sprintf("%s#%s",sample,barcode)] %>%
  .[sample%in%opts$samples & celltype.mapped%in%opts$celltypes]

stopifnot(sample_metadata$archR_cell %in% rownames(ArchRProject))

##################
## Subset ArchR ##
##################

ArchRProject.filt <- ArchRProject[sample_metadata$archR_cell]
table(getCellColData(ArchRProject.filt,"Sample")[[1]])
table(getCellColData(ArchRProject.filt,"celltype.mapped")[[1]])

######################################
## Identify celltype-specific peaks ##
######################################

markersPeaks <- getMarkerFeatures(
  ArchRProj = ArchRProject.filt, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "celltype.mapped",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

# returns a SummarizedExperiment object
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")
names(markerList)
head(markerList[[1]])

# returns a GRanges object
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)
dt <- names(markerList) %>% 
  map(function(i) as.data.table(markerList[[i]]) %>% .[,celltype:=i]) %>% 
  rbindlist %>%
  setnames("seqnames","chr") %>%
  .[,c("strand","width"):=NULL]
head(dt)

# save
# outfile <- paste0(io$basedir,"/results/markers/markers.tsv.gz")
# fwrite(dt, outfile, sep="\t")

#######################
## Plot Marker genes ##
#######################

# Plot heatmap
plotMarkerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1",
  scaleRows = TRUE,
  clusterCols = TRUE,
  transpose = TRUE,
  labelRows = FALSE
)

# MA plot
markerPlot(
  seMarker = markersPeaks, 
  name = "C1", 
  cutOff = "FDR <= 0.1 & Log2FC >= 1", 
  plotAs = "MA"
)

# Volcano plot
markerPlot(
  seMarker = markersPeaks, 
  name = "C1", 
  cutOff = "FDR <= 0.1 & Log2FC >= 1", 
  plotAs = "Volcano"
)

