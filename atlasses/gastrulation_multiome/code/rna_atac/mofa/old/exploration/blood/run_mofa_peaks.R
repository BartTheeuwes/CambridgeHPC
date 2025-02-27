
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

################
## Define I/O ##
################

io$mofa.output <- paste0(io$basedir, "/results/rna_atac/mofa/NMPs/mofa_model.rds")

# Options
opts$samples <- c(
  "E7.5_rep1",
  "E7.5_rep2",
  "E8.5_rep1",
  "E8.5_rep2"
)
opts$celltypes = c(
  # "Intermediate_mesoderm",
  # "Caudal_epiblast",
  # "Paraxial_mesoderm",
  "Spinal_cord",
  "NMP",
  "Somitic_mesoderm"
  # "Forebrain_Midbrain_Hindbrain",
)
opts$npeaks <- 5e3

########################
## Load cell metadata ##
########################

io$metadata <- paste0(io$basedir, "/sample_metadata.txt.gz")

sample_metadata <- fread(io$metadata) %>%
  .[pass_atacQC==TRUE & pass_rnaQC==TRUE] %>%
  .[sample%in%opts$samples & celltype.mapped%in%opts$celltypes]

##########################
## Subset ArchR project ##
##########################

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
## Feature selection on ATAc peaks ##
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
hvgs <- rownames(decomp)[decomp$p.value <= 0.10]
sce_filt <- sce[hvgs,]

rna.mtx <- assay(sce_filt,"logcounts")
dim(rna.mtx)

###########################
## Prepare data for MOFA ##
###########################

cells <- intersect(colnames(rna.mtx),colnames(tfidf.peak.mtx))

rna.mtx <- rna.mtx[,cells]
tfidf.peak.mtx <- tfidf.peak.mtx[,cells]

########################
## Create MOFA object ##
########################

MOFAobject <- create_mofa_from_matrix(
  list(
  "RNA" = rna.mtx,
  "ATAC_peaks" = tfidf.peak.mtx
  )
)
MOFAobject

####################
## Define options ##
####################

# Data options
data_opts <- get_default_data_options(MOFAobject)
data_opts$use_float32 <- TRUE

# Model options
model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 5
model_opts$spikeslab_weights <- TRUE

# Training options
train_opts <- get_default_training_options(MOFAobject)
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

MOFAobject <- run_mofa(MOFAobject)

#########################
# Add samples metadata ##
#########################

metadata.to.mofa <- sample_metadata %>% copy %>%
  setnames("sample","batch") %>% setnames("cell","sample") %>%
  .[sample%in%unlist(samples_names(MOFAobject))] %>%
  setkey(sample) %>% .[unlist(samples_names(MOFAobject))]
samples_metadata(MOFAobject) <- metadata.to.mofa

##########
## Save ##
##########

saveRDS(MOFAobject, io$mofa.output)