R.utils::sourceDirectory("/Users/ricard/git/ArchR/R/", verbose=T, modifiedOnly=FALSE)
R.utils::sourceDirectory("/homes/ricard/git/ArchR/R/", verbose=T, modifiedOnly=FALSE)

# peaks.bed <- fread(io$archR.peakSet.bed) %>% .[,foo:=sprintf("%s_%s_%s",V1,V2,V3)]
# foo <- sprintf("%s_%s_%s",seqnames(ArchRProject.filt@peakSet),start(ArchRProject.filt@peakSet),end(ArchRProject.filt@peakSet))

ArchRProj = ArchRProject.filt
groupBy = "celltype.predicted"
useGroups = NULL
bgdGroups = NULL
useMatrix = "PeakMatrix"
bias = c("TSSEnrichment","log10(nFrags)")
normBy = NULL
testMethod = "wilcoxon"
maxCells = 500
scaleTo = 10^4
threads = 1
k = 100
bufferRatio = 0.8
binarize = FALSE
useSeqnames = NULL
verbose = TRUE
logFile = createLogFile("getMarkerFeatures")
