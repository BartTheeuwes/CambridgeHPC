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
# io$metadata <- paste0(io$basedir,"/processed/atac/archR/sample_metadata_after_archR.txt.gz")
io$outdir <- paste0(io$basedir,"/results/atac/archR/marker_peaks/genes")

# Options
opts$samples <- c(
  "E7.5_rep1",
  "E7.5_rep2",
  "E8.5_rep1",
  "E8.5_rep2"
)

########################
## Load cell metadata ##
########################

sample_metadata <- fread(io$metadata) %>%
  .[pass_atacQC==TRUE] %>%
  .[sample%in%opts$samples]
stopifnot(sample_metadata$archR_cell %in% rownames(ArchRProject))

# subset celltypes with sufficient number of cells
opts$min.cells <- 100
sample_metadata <- sample_metadata %>%
  .[,N:=.N,by=c("celltype.predicted")] %>% .[N>opts$min.cells] %>% .[,N:=NULL]

# to.save <- sample_metadata %>% .[,.(N=.N),by=c("celltype.predicted")]
# fwrite(to.save, paste0(io$outdir,"/celltype_numbers.txt.gz"))
# sample_metadata[celltype.predicted!=celltype.mapped,c("sample","cell","celltype.predicted","celltype.mapped")] %>% View

##################
## Subset ArchR ##
##################

ArchRProject.filt <- ArchRProject[sample_metadata$archR_cell,]
table(getCellColData(ArchRProject.filt,"Sample")[[1]])
table(getCellColData(ArchRProject.filt,"celltype.predicted")[[1]])

######################################
## Identify celltype-specific peaks ##
######################################

# Load pre-computed
# markersPeaks <- readRDS(paste0(io$outdir,"/genes/markersPeaks_summarized_experiment.rds"))

markersPeaks <- getMarkerFeatures(
  ArchRProj = ArchRProject.filt, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "celltype.predicted",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
saveRDS(markersPeaks, paste0(io$outdir,"/markersPeaks_summarized_experiment.rds"))

#################
## Get Markers ##
#################

# returns a SummarizedExperiment object
# markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= Inf & Log2FC >= -Inf")
# markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)

dt <- names(markerList) %>% 
  map(function(i) as.data.table(markerList[[i]]) %>% .[,celltype:=i]) %>% 
  rbindlist %>%
  setnames("seqnames","chr") %>%
  setnames("name","gene") %>%
  .[,c("strand","width","idx"):=NULL]
head(dt)

# save
fwrite(dt, paste0(io$outdir,"/markers_all.tsv.gz"), sep="\t")



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

