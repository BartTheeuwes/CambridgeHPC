################
## Clustering ##
################

ArchRProject.filt <- addClusters(
  input = ArchRProject.filt,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters"
)
head(ArchRProject.filt@cellColData$Clusters)