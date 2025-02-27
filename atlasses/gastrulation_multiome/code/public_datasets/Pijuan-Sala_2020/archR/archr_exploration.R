library(ArchR)
library(ggpubr)
library(Matrix)

#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/public_datasets/Pijuan-Sala_2020/settings.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/public_datasets/Pijuan-Sala_2020/settings.R")
} else {
  stop("Computer not recognised")
}
setwd(io$archR.directory)

################
## Define I/O ##
################

io$outdir <- paste0(io$archR.directory,"/pdf")

####################
## Define options ##
####################

addArchRThreads(threads = 2) 
addArchRGenome("mm10")

########################
## Load ArchR project ##
########################

ArchRProject <- loadArchRProject(io$archR.directory)

# rename cells and samples (not allowed)
# ArchRProject$cellNames <- ArchRProject$cellNames %>% stringr::str_replace_all("#","_") %>% stringr::str_replace_all("-1","")
# ArchRProject$Sample <- stringr::str_replace_all(ArchRProject$Sample,"multiome","E8.5_rep")

#####################
## Update metadata ##
#####################

sample_metadata[,archR_cell:=paste0("E8.5_gastrulation#",cell)]
mean(sample_metadata$archR_cell%in%ArchRProject$cellNames)
mean(ArchRProject$cellNames%in%sample_metadata$archR_cell)

# subset cells
cells <- ArchRProject$cellNames[ArchRProject$cellNames %in% sample_metadata$archR_cell]
ArchRProject.filt <- ArchRProject[cells,]

sample_metadata.to.archR <- sample_metadata %>%
  .[archR_cell%in%cells] %>% 
  setkey(archR_cell) %>% .[cells] %>%
  as.data.frame() %>%
  tibble::column_to_rownames("archR_cell")

ArchRProject.filt
dim(sample_metadata.to.archR)

colnames(getCellColData(ArchRProject))

for (i in colnames(sample_metadata.to.archR)) {
  ArchRProject.filt <- addCellColData(ArchRProject.filt,
    data = sample_metadata.to.archR[[i]], 
    name = i,
    cells = rownames(sample_metadata.to.archR)
  )
}

#####################
## Quality control ##
#####################

# Number of unique nuclear fragments (i.e. not mapping to mitochondrial DNA).
plotGroups(ArchRProject.filt, groupBy = "celltype", colorBy = "cellColData", name = "TSSEnrichment", plotAs = "violin", alpha = 0.4, addBoxPlot = TRUE)
plotGroups(ArchRProject.filt, groupBy = "celltype", colorBy = "cellColData", name = "NucleosomeRatio", plotAs = "violin", alpha = 0.4, addBoxPlot = TRUE)
plotGroups(ArchRProject.filt, groupBy = "celltype", colorBy = "cellColData", name = "BlacklistRatio", plotAs = "violin", alpha = 0.4, addBoxPlot = TRUE)
plotGroups(ArchRProject.filt, groupBy = "celltype", colorBy = "cellColData", name = "nFrags", plotAs = "violin", alpha = 0.4, addBoxPlot = TRUE)

# TSS enrichment score. Low signal-to-background ratio is often attributed to dead or dying cells which have de-chromatinzed DNA which allows for random transposition genome-wide.
plotFragmentSizes(ArchRProject.filt)
to.plot <- plotFragmentSizes(ArchRProject.filt, returnDF=T) %>% as.data.table() %>%
  merge(sample_metadata.to.archR,by="cell")
ggline(to.plot, x="fragmentSize", y="fragmentPercent", plot_type="l")

# Fragment size distribution. Due to nucleosomal periodicity, we expect to see depletion of fragments that are the length of DNA wrapped around a nucleosome (approximately 147 bp).
plotTSSEnrichment(ArchRProject.filt)
to.plot <- plotTSSEnrichment(ArchRProject.filt, returnDF=T) %>% as.data.table() %>%
  melt(id.vars=c("sampleName","x"))
ggline(to.plot, x="x", y="value", color="sampleName", plot_type="l") +
  facet_wrap(~variable, scales="free_y") +
  labs(x="") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.title = element_blank()
  )

# Mitochondrial coverage

#####################
## Filter doublets ##
#####################

# To predict which cells are actually doublets:
# (1) synthesize in silico doublets from the data by mixing the reads from thousands of combinations of individual cells. 
# (2) project these synthetic doublets into the UMAP embedding 
# (3) Identify their nearest neighbor. 
# By iterating this procedure thousands of times, we can identify cells in our data whose signal looks very similar to synthetic doublets.

# Adding doublet scores will create plots in the QualityControl directory.
doubScores <- addDoubletScores(ArrowFiles,
  useMatrix = "TileMatrix",
  k = 10, # refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP",
  LSIMethod = 1
)
class(doubScores)
names(doubScores[[1]])
head(doubScores[[1]][[1]])

# ArchRProject.filt <- filterDoublets(ArchRProject.filt, cutEnrich = 1, cutScore = -Inf, filterRatio = 1)




#################
## Gene scores ##
#################

# AchR model: 
# (1) Accessibility within the entire gene body contributes to the gene score.
# (2) An exponential weighting function that accounts for the activity of putative distal regulatory elements in a distance-dependent fashion.
# (3) Imposed gene boundaries that minimizes the contribution of unrelated regulatory elements to the gene score.

# Gene scores are calculated for each Arrow file at the time of creation if the parameter addGeneScoreMat is set to TRUE - this is the default behavior
# returns a SummarizedExperiment class
gene.score.matrix <- getMatrixFromProject(ArchRProject.filt, useMatrix = "GeneScoreMatrix")
assay(gene.score.matrix)[1:3,1:3]
head(colData(gene.score.matrix))
head(rowData(gene.score.matrix))

##################
## Gene markers ##
##################

table(ArchRProject.filt$celltype)

# It is important to note that not all genes behave well with gene scores. 
# In particular, genes that reside in very gene-dense areas can be problematic. 
# Thus, it is always best to sanity check all gene score analyses by looking at sequencing tracks which is described in a later chapter.
markersGS <- getMarkerFeatures(
  ArchRProj = ArchRProject.filt, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "celltype",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
dim(markersGS@assays@data$Log2FC)


################
## Pseudobulk ##
################

# https://www.ArchRProject.filt.com/bookdown/how-does-archr-make-pseudo-bulk-replicates.html

# This function will merge cells within each designated cell group for the generation of pseudo-bulk replicates 
# and then merge these replicates into a single insertion coverage file.
# Output: creates files in archR/GroupCoverages/celltype: [X]._.Rep[Y].insertions.coverage.h5
ArchRProject.filt <- addGroupCoverages(ArchRProject.filt, groupBy = "celltype")


##########
## Save ##
##########

getOutputDirectory(ArchRProject.filt)
saveArchRProject(ArchRProject.filt)

##########
## Load ##
##########

io$outdir  <- getOutputDirectory(ArchRProject)
ArchRProject <- loadArchRProject(io$outdir)

# getAvailableMatrices(ArchRProject)