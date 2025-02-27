
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

# io$metadata <- paste0(io$basedir,"/processed/atac/archR/sample_metadata_after_archR.txt.gz")
io$metadata <- paste0(io$basedir,"/sample_metadata.txt.gz")
# io$metadata <- paste0(io$basedir,"/results/atac/archR/celltype_assignment/sample_metadata_after_archR.txt.gz")
io$outdir <- paste0(io$basedir,"/results/atac/archR/denoising")

opts$samples <- c(
  "E7.5_rep1",
  "E7.5_rep2",
  "E8.5_rep1",
  "E8.5_rep2"
)

# LSI parameters
opts$lsi.iterations = 2
opts$lsi.cluster.resolution = 2
opts$lsi.varFeatures <- 25000
opts$lsi.dims <- 30
opts$matrix <- "PeakMatrix"

# UMAP parameters
opts$umap.seed <- 42
opts$umap.neighbours <- 25
opts$umap.minDist <- 0.3

opts$batch.correction <- FALSE

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

ArchRProject.filt <- ArchRProject[sample_metadata$archR_cell]

###########################
## Latent Semantic Index ##
###########################

# Iterative LSI: two iterations
ArchRProject.filt <- addIterativeLSI(ArchRProject.filt,
  useMatrix = opts$matrix, 
  name = "IterativeLSI", 
  firstSelection = "Top", # "Top" or "Var"
  depthCol = "nFrags",
  iterations = opts$lsi.iterations, 
  clusterParams = list(
    resolution = opts$lsi.cluster.resolution, 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = opts$lsi.varFeatures, 
  force = TRUE
)

# Correlation between latent dimensions and number of peaks
cor(getReducedDims(ArchRProject.filt, "IterativeLSI"),ArchRProject.filt$nFrags)[,1] %>% abs %>% sort(decreasing = T)


##########
## Save ##
##########
  
# save umap coordinates
# to.save <- getEmbedding(ArchRProject.filt,"UMAP") %>%
#   as.data.table(keep.rownames = T) %>%
#   setnames(c("archR_cell","V1","V2"))
# fwrite(to.save, sprintf("%s/archr_umap_coordinates_%s_LSIiter%s_nfeatures%s_neighb%s_mindist%s.txt",io$outdir,opts$matrix,opts$lsi.iterations, opts$lsi.varFeatures, opts$umap.neighbours,opts$umap.minDist))

##########
## TEST ##
##########

# getAvailableMatrices(ArchRProject.filt)
# deviations.se <- getMatrixFromProject(ArchRProject.filt, "DeviationMatrix")
