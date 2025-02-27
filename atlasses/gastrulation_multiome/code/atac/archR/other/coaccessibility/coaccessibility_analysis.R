# https://www.archrproject.com/bookdown/co-accessibility-with-archr.html

########################
## Load ArchR project ##
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

io$metadata <- paste0(io$basedir,"/processed/atac/archR/sample_metadata_after_archR.txt.gz")
# io$metadata <- paste0(io$basedir,"/results/atac/archR/celltype_assignment/sample_metadata_after_archR.txt.gz")
io$outdir <- paste0(io$basedir,"/results/atac/archR/coaccessibility")

opts$samples <- c(
  # "E7.5_rep1",
  # "E7.5_rep2",
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

##################
## Subset ArchR ##
##################

# subset celltypes with sufficient number of cells
opts$min.cells <- 100
sample_metadata <- sample_metadata %>%
  .[,N:=.N,by=c("celltype.predicted")] %>% .[N>opts$min.cells] %>% .[,N:=NULL]
table(sample_metadata$celltype.predicted)

ArchRProject.filt <- ArchRProject[sample_metadata$archR_cell]
ArchRProject.filt@sampleColData <- ArchRProject.filt@sampleColData[opts$samples,,drop=F]
table(getCellColData(ArchRProject.filt,"Sample")[[1]])
table(getCellColData(ArchRProject.filt,"celltype.predicted")[[1]])

###################################
## Load dimensionality reduction ##
###################################

getReducedDims(ArchRProject.filt)

# LSI parameters
opts$lsi.iterations <- 1
opts$lsi.cluster.resolution <- 2
opts$lsi.varFeatures <- 25000
opts$lsi.dims <- 30
opts$matrix <- "PeakMatrix"


# Iterative LSI: two iterations
ArchRProject.filt <- addIterativeLSI(
   ArchRProject.filt,
   useMatrix = opts$matrix, 
   name = "IterativeLSI", 
   iterations = opts$lsi.iterations, 
   clusterParams = list(
     resolution = opts$lsi.cluster.resolution, 
     sampleCells = 10000, 
     n.start = 10
   ), 
   varFeatures = opts$lsi.varFeatures, 
   dimsToUse = 1:opts$lsi.dims,
   force = TRUE
)

#####################
## Define Peak Set ##
#####################

io$peaks <- "/Users/ricard/data/gastrulation_multiome_10x/processed/atac/archR/PeakCalls/E8.5_peaks_archR_macs2.bed"

# Load Peak Set
peaks.dt <- fread(io$peaks) %>%
  setnames(c("chr","start","end")) %>%
  .[,id:=sprintf("%s:%s-%s",chr,start,end)]

peaks.granges <- makeGRangesFromDataFrame(peaks.dt, keep.extra.columns = T)
ArchRProject.filt <- addPeakSet(ArchRProject.filt, peaks.granges)

head(getPeakSet(ArchRProject.filt))

##################################
## Run coaccessibility analysis ##
##################################

ArchRProject.filt <- addCoAccessibility(
    ArchRProj = ArchRProject.filt,
    reducedDims = "IterativeLSI",
    dimsToUse = 1:opts$lsi.dims,
    corCutOff = 0.75
)

##############
## Analysis ##
##############

# Columns
# - "queryHits" and "subjectHits" columns denote the index of the two peaks that were found to be correlated. 
# - "correlation": numeric correlation of the accessibility between those two peaks.
coacc.df <- getCoAccessibility(
    ArchRProj = ArchRProject.filt,
    corCutOff = 0.5,
    resolution = 1,      # the bp resolution to return loops as. This helps with overplotting of correlated regions.
    returnLoops = FALSE  # return the co-accessibility signal as a GRanges "loops" object designed for use with the ArchRBrowser()
    
    
)
head(coacc.df)

coacc.df2 <- getCoAccessibility(
    ArchRProj = ArchRProject.filt,
    corCutOff = 0.5,
    resolution = 10000,
    returnLoops = TRUE
)

head(coacc.df2[[1]])

##########
## Save ##
##########

to.save <- coacc.df %>% as.data.table
fwrite(to.save, sprintf("%s/archr_coaccessibility_scores.txt",io$outdir))

#############################
## Plotting browser tracks ##
#############################

p <- plotBrowserTrack(
    ArchRProj = ArchRProject.filt, 
    groupBy = "celltype.predicted", 
    geneSymbol = "T", 
    upstream = 50000,
    downstream = 50000,
    loops = coacc.df2
)

grid::grid.newpage()
grid::grid.draw(p[[1]])
