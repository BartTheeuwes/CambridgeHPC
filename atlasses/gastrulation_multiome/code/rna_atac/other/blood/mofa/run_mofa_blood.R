
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(MOFA2))


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
io$metadata <- paste0(io$basedir, "/sample_metadata.txt.gz")
# io$outdir <- paste0(io$basedir, "/results/atac/signac/dimensionality_reduction")
io$mofa.output <- paste0(io$basedir, "/results/rna_atac/mofa/mofa_model_blood_rna.rds")

# Options
opts$celltypes <- c(
  # "Cardiomyocytes",
  "Allantois",
  # "ExE_mesoderm",
  # "Mesenchyme",
  "Haematoendothelial_progenitors",
  "Endothelium",
  "Blood_progenitors_1",
  "Blood_progenitors_2",
  "Erythroid1",
  "Erythroid2",
  "Erythroid3"
)
opts$npeaks <- 5e3

########################
## Load cell metadata ##
########################

sample_metadata <- fread(io$metadata) %>%
  .[pass_atacQC==TRUE & pass_rnaQC==TRUE] %>%
  .[celltype.predicted%in%opts$celltypes]

##################
## Subset ArchR ##
##################

stopifnot(sample_metadata$archR_cell %in% rownames(ArchRProject))
ArchRProject.filt <- ArchRProject[sample_metadata$archR_cell]

######################
## Load peak matrix ##
#####################

# Get Peak Matrix from ArchR
atac.peak.se <- getMatrixFromProject(ArchRProject.filt, useMatrix="PeakMatrix")
dim(atac.peak.se)

# Define peak names
peak_names <- rowRanges(atac.peak.se) %>% as.data.table %>% .[,id:=sprintf("%s_%s_%s",seqnames,start,end)] %>% .$id
rownames(atac.peak.se) <- peak_names

# Rename cells 
colnames(atac.peak.se) <- sample_metadata %>% 
  .[archR_cell%in%colnames(atac.peak.se)] %>%
  setkey(archR_cell) %>% .[colnames(atac.peak.se)] %>% .$cell

#####################################
## Feature selection on ATAC peaks ##
#####################################

# Load peak variability estimates
peak.variability.dt <- fread(io$archR.peak.variability)

# Define highly variable peaks
peaks <- peak.variability.dt %>% 
  setorder(-variance_pseudobulk) %>% 
  head(n=opts$npeaks) %>% 
  .$peak

# Subset
atac.peak.se.filt <- atac.peak.se[peaks,]
dim(atac.peak.se.filt)

# Get matrix
peak.mtx <- assay(atac.peak.se.filt, "PeakMatrix")

# TFIDF normalisation
tfidf.peak.mtx <- tfidf(peak.mtx)
dim(tfidf.peak.mtx)

# hist(tfidf.peak.mtx[1:10000,1:1000] %>% as.matrix)

##########################
## Load cisTopic matrix ##
##########################

cistopic.mtx <- fread(io$cistopic) %>%
  .[,c("UMAP1","UMAP2"):=NULL] %>%
  matrix.please
dim(cistopic.mtx)

##############################
## Load chromVAR deviations ##
##############################

# Get Deviation Matrix from ArchR
# atac.deviation.mtx <- getMatrixFromProject(ArchRProject.filt, useMatrix="DeviationMatrix")

# Load Deviation Matrix as a SummarizedExperiment object
atac.deviation.mtx <- readRDS(io$archR.deviations.se) %>% assay(.,"z")

# Rename features
rownames(atac.deviation.mtx) <- paste0("chromVAR_",rownames(atac.deviation.mtx))

# Subset cells
atac.deviation.mtx <- atac.deviation.mtx[,colnames(atac.deviation.mtx) %in% sample_metadata$archR_cell]

# Rename cells 
colnames(atac.deviation.mtx) <- sample_metadata %>% 
  .[archR_cell%in%colnames(atac.deviation.mtx)] %>%
  setkey(archR_cell) %>% .[colnames(atac.deviation.mtx)] %>% .$cell

##############################
## Load RNA expression data ##
##############################

sce <- load_SingleCellExperiment(io$sce, normalise = TRUE, cells = sample_metadata$cell)
dim(sce)

# Add sample metadata to the colData of the SingleCellExperiment
colData(sce) <- sample_metadata %>% as.data.frame %>% tibble::column_to_rownames("cell") %>%
  .[colnames(sce),] %>% DataFrame()

# Feature selection
decomp <- modelGeneVar(sce)
decomp <- decomp[decomp$mean > 0.001,]
hvgs <- rownames(decomp)[decomp$p.value <= 0.01]
sce_filt <- sce[hvgs,]

rna.mtx <- assay(sce_filt,"logcounts")
dim(rna.mtx)

###########################
## Prepare data for MOFA ##
###########################

cells <- Reduce(intersect,list(colnames(rna.mtx),colnames(tfidf.peak.mtx)))
# cells <- Reduce(intersect,list(colnames(rna.mtx),colnames(atac.deviation.mtx),colnames(tfidf.peak.mtx)))

rna.mtx <- rna.mtx[,cells]
# atac.deviation.mtx <- atac.deviation.mtx[,cells]
tfidf.peak.mtx <- tfidf.peak.mtx[,cells]

# hist(rna.mtx[1:500,] %>% as.matrix)
# hist(atac.deviation.mtx[1:500,] %>% as.matrix)


########################
## Create MOFA object ##
########################

# Sanity check
stopifnot(colnames(rna.mtx) == colnames(atac.deviation.mtx))
stopifnot(colnames(atac.deviation.mtx) == colnames(tfidf.peak.mtx))

# MOFAobject <- create_mofa_from_SingleCellExperiment(sce_filt, extract_metadata = T)
# views_names(MOFAobject) <- "RNA"

MOFAobject <- create_mofa_from_matrix(
  list(
  "RNA" = rna.mtx, 
  # "ATAC_chromVAR" = atac.deviation.mtx, 
  "ATAC_peaks" = tfidf.peak.mtx
  )
)
MOFAobject

# Visualise data structure
# plot_data_overview(MOFAobject)

####################
## Define options ##
####################

# Data options
data_opts <- get_default_data_options(MOFAobject)
data_opts$use_float32 <- TRUE

# Model options
model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 10
model_opts$spikeslab_weights <- FALSE

# Training options
train_opts <- get_default_training_options(MOFAobject)
train_opts$convergence_mode <- "fast"
# train_opts$maxiter <- 5

#########################
## Prepare MOFA object ##
#########################

MOFAobject <- prepare_mofa(
  MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

#####################
## Train the model ##
#####################

mofa <- run_mofa(MOFAobject)

#########################
# Add samples metadata ##
#########################

metadata.to.mofa <- sample_metadata %>%
  setnames("sample","batch") %>% setnames("cell","sample") %>%
  .[sample%in%unlist(samples_names(mofa))] %>%
  setkey(sample) %>% .[unlist(samples_names(mofa))]
samples_metadata(mofa) <- metadata.to.mofa

##########
## Save ##
##########

saveRDS(mofa, io$mofa.output)
