#############################
## Peak calling with MACS2 ##
#############################

# https://www.ArchRProject.filt.com/bookdown/calling-peaks-with-archr.html

pathToMacs2 <- findMacs2()

# This function will get insertions from coverage files, call peaks, and merge peaks to get a "Union Reproducible Peak Set".
ArchRProject.filt <- addReproduciblePeakSet(
  ArchRProj = ArchRProject.filt, 
  groupBy = "celltype", 
  peakMethod = "Macs2",
  excludeChr = c("chrM", "chrY"),
  pathToMacs2 = pathToMacs2,
  cutOff = 0.05,
  extendSummits = 300,
  plot = TRUE
)


########################
## Create peak matrix ##
########################

ArchRProject.filt <- addPeakMatrix(ArchRProject.filt)
getAvailableMatrices(ArchRProject.filt)

foo <- getMatrixFromProject(ArchRProject.filt, useMatrix = "GeneScoreMatrix")
mtx <- getMatrixFromProject(ArchRProject.filt, useMatrix = "PeakMatrix")@assays@data[[1]]

# Save matrix.mtx.gz
outfile <- "/Users/ricard/data/gastrulation_multiome_10x/public_datasets/Pijuan-Sala_2020/data/peak_matrix/archR/matrix.mtx.gz"
Matrix::writeMM(mtx@assays@data$PeakMatrix, file=outfile)

# Save barcodes.tsv.gz
outfile <- "/Users/ricard/data/gastrulation_multiome_10x/public_datasets/Pijuan-Sala_2020/data/peak_matrix/archR/barcodes.tsv.gz"
fwrite(as.data.table(gsub("E8.5_gastrulation#","", colnames(mtx))), outfile, col.names=F)

# Save features.tsv.gz and peak metadata
dt <- getPeakSet(ArchRProject.filt) %>% as.data.table() %>% setnames(c("seqnames"),c("chr")) %>%
  .[,c("strand","idx","nearestTSS","distToTSS","GroupReplicate","replicateScoreQuantile","groupScoreQuantile","Reproducibility"):=NULL]
outfile <- "/Users/ricard/data/gastrulation_multiome_10x/public_datasets/Pijuan-Sala_2020/data/peak_matrix/archR/peak_metadata.tsv.gz"
fwrite(dt, outfile, sep="\t")

foo <- dt[,c("chr","start","end")] %>% .[,foo:=sprintf("%s_%s_%s",chr,start,end)] %>% .[,"foo"]
outfile <- "/Users/ricard/data/gastrulation_multiome_10x/public_datasets/Pijuan-Sala_2020/data/peak_matrix/archR/features.tsv.gz"
fwrite(foo, outfile, sep="\t", col.names = F)

# save peaks in bed format
fwrite(dt[,c("chr","start","end")], paste0(io$archR.directory,"/PeakCalls/peaks_archR_macs2.bed"), sep="\t", col.names = F)

