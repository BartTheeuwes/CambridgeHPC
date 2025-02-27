# https://www.archrproject.com/bookdown/calculating-gene-scores-in-archr.html

#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
  source("/Users/ricard/gastrulation_multiome_10x/utils.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
  source("/homes/ricard/gastrulation_multiome_10x/utils.R")
}

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

# I/O
# io$metadata <- paste0(io$basedir,"/sample_metadata.txt.gz")
# io$outdir <- paste0(io$basedir,"/results/atac/archR/gene_scores")

#####################
## Define gene set ##
#####################

# getGeneAnnotation(ArchRProject)
genes.gr <- getGenes(ArchRProject)

getAvailableMatrices(ArchRProject)

#########################################
## Add Gene Scores using default model ##
#########################################

# Note that this will add the matrices to the arrowFiles
addGeneScoreMatrix(
  input = ArchRProject,
  genes = genes.gr,
  geneModel = "exp(-abs(x)/5000) + exp(-1)", # string should be a function of x, where x is the distance from the TSS.
  matrixName = "GeneScoreMatrix",
  extendUpstream = c(1000, 1e+05),
  extendDownstream = c(1000, 1e+05),  # The minimum and maximum number of bp downstream of the transcription termination site to consider for gene activity score calculation.
  geneUpstream = 5000,                # Number of bp upstream the gene to extend the gene body.
  geneDownstream = 0,                 # Number of bp downstream the gene to extend the gene body.
  tileSize = 500,                     # The size of the tiles used for binning counts prior to gene activity score calculation.
  geneScaleFactor = 5,                # A numeric scaling factor to weight genes based on the inverse of their length 
  scaleTo = 10000,                    # Each column in the calculated gene score matrix will be normalized
  excludeChr = c("chrY", "chrM")
)

#############################################
## Add Gene Scores ignoring distal regions ##
#############################################

# Note that this will add the matrices to the arrowFiles
addGeneScoreMatrix(
  input = ArchRProject,
  genes = genes.gr,
  geneModel = "1",                    # string should be a function of x, where x is the distance from the TSS.
  matrixName = "GeneScoreMatrix_nodistal",
  extendUpstream = c(0,0),
  extendDownstream = c(0,0),          # The minimum and maximum number of bp downstream of the transcription termination site to consider for gene activity score calculation.
  geneUpstream = 2000,                # Number of bp upstream the gene to extend the gene body.
  geneDownstream = 0,                 # Number of bp downstream the gene to extend the gene body.
  tileSize = 100,                     # The size of the tiles used for binning counts prior to gene activity score calculation.
  geneScaleFactor = 5,                # A numeric scaling factor to weight genes based on the inverse of their length 
  scaleTo = 10000,                    # Each column in the calculated gene score matrix will be normalized
  excludeChr = c("chrY", "chrM")
)


addGeneScoreMatrix(
  input = ArchRProject,
  genes = genes.gr,
  useTSS = TRUE,
  extendTSS = TRUE,
  geneModel = "1",                    # string should be a function of x, where x is the distance from the TSS.
  matrixName = "GeneScoreMatrix_nodistal_v2",
  extendUpstream = c(0,0),
  extendDownstream = c(0,0),          # The minimum and maximum number of bp downstream of the transcription termination site to consider for gene activity score calculation.
  geneUpstream = 500,                # Number of bp upstream the gene to extend the gene body.
  geneDownstream = 100,                 # Number of bp downstream the gene to extend the gene body.
  tileSize = 100,                     # The size of the tiles used for binning counts prior to gene activity score calculation.
  geneScaleFactor = 1,                # A numeric scaling factor to weight genes based on the inverse of their length 
  scaleTo = 10000,                    # Each column in the calculated gene score matrix will be normalized
  excludeChr = c("chrY", "chrM"),
  force = TRUE
)


##########
## Test ##
##########

# atac.GeneScoreMatrix.se <- getMatrixFromProject(ArchRProject, binarize = FALSE, useMatrix = "GeneScoreMatrix_nodistal_v2")
# rownames(atac.GeneScoreMatrix.se) <- rowData(atac.GeneScoreMatrix.se)$name
# tmp <- rowMeans(assay(atac.GeneScoreMatrix.se)) %>% sort(decreasing = T)
# head(tmp)
# tmp[c("Auts2","Nav2")]
# tmp[c("Actb","Polr2f")]
