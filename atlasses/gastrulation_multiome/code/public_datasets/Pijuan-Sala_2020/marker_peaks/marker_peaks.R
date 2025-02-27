######################################
## Identify celltype-specific peaks ##
######################################

# https://www.ArchRProject.com/bookdown/identifying-marker-peaks-with-archr.html

markersPeaks <- getMarkerFeatures(
  ArchRProj = ArchRProject, 
  useMatrix = "PeakMatrix", 
  groupBy = "celltype",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

# returns a SummarizedExperiment object
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")
names(markerList)
head(markerList[["C1"]])

# returns a GRanges object
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)
dt <- names(markerList) %>% 
  map(function(i) as.data.table(markerList[[i]]) %>% .[,celltype:=i]) %>% 
  rbindlist %>%
  setnames("seqnames","chr") %>%
  .[,c("strand","width"):=NULL]
head(dt)

# save
outfile <- paste0(io$basedir,"/results/markers/markers.tsv.gz")
fwrite(dt, outfile, sep="\t")

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

