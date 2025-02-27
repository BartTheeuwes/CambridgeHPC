library(ArchR)

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

################
## Define I/O ##
################

io$archR.directory <- paste0(io$basedir,"/data/processed/archR")

setwd(io$archR.directory)

####################
## Define options ##
####################

addArchRThreads(threads = 2) 

# The precompiled version of the hg19 genome in ArchR uses 
# BSgenome.Hsapiens.UCSC.hg19, TxDb.Hsapiens.UCSC.hg19.knownGene, org.Hs.eg.db, 
# and a blacklist that was merged using ArchR::mergeGR() from the hg19 v2 blacklist regions 
# and from mitochondrial regions that show high mappability to the hg19 nuclear genome

# Add Genome built
addArchRGenome("mm10")

########################
## create Arrow Files ##
########################

# Steps:
# (1) Read accessible fragments from the provided input files.
# (2) Calculate quality control information for each cell (i.e. TSS enrichment scores and nucleosome info).
# (3) Filter cells based on quality control parameters.
# (4) Create a genome-wide TileMatrix using 500-bp bins.
# (5) Create a GeneScoreMatrix using the custom geneAnnotation that was defined when we called addArchRGenome().

ArrowFiles <- createArrowFiles(
  inputFiles = io$fragments,
  sampleNames = "E8.5_gastrulation",
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  excludeChr = c("chrM", "chrY")
)

############################
## create an ArchRProject ##
############################

ArchRProject <- ArchRProject(ArrowFiles, 
  outputDirectory = io$archR.directory,
  copyArrows = FALSE # This is recommended so that if you modify the Arrow files you have an original copy for later usage.
)

#####################
## Update metadata ##
#####################

# rename cells and samples (not allowed)
# ArchRProject$cellNames <- ArchRProject$cellNames %>% stringr::str_replace_all("#","_") %>% stringr::str_replace_all("-1","")
# ArchRProject$Sample <- stringr::str_replace_all(ArchRProject$Sample,"multiome","E8.5_rep")


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

##########
## Save ##
##########

saveArchRProject(ArchRProject)
