
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(MOFA2))

#####################
## Define settings ##
#####################

# Load default settings
if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
  source("/Users/ricard/gastrulation_multiome_10x/utils.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
  source("/homes/ricard/gastrulation_multiome_10x/utils.R")
} else {
  stop("Computer not recognised")
}


# I/O
io$mofa.output <- paste0(io$basedir, "/results/rna_atac/mofa/all_cells/mofa_allcells_lsi_noExE.rds")

# Options
opts$samples <- c(
  "E7.5_rep1",
  "E7.5_rep2",
  "E8.0_rep1",
  "E8.0_rep2",
  "E8.5_rep1",
  "E8.5_rep2"
)

opts$celltypes = c(
  "Epiblast",
  "Primitive_Streak",
  "Caudal_epiblast",
  # "PGC",
  "Anterior_Primitive_Streak",
  "Notochord",
  "Def._endoderm",
  "Gut",
  "Nascent_mesoderm",
  "Mixed_mesoderm",
  "Intermediate_mesoderm",
  "Caudal_Mesoderm",
  "Paraxial_mesoderm",
  "Somitic_mesoderm",
  "Pharyngeal_mesoderm",
  "Cardiomyocytes",
  "Allantois",
  "ExE_mesoderm",
  "Mesenchyme",
  "Haematoendothelial_progenitors",
  "Endothelium",
  "Blood_progenitors_1",
  "Blood_progenitors_2",
  "Erythroid1",
  "Erythroid2",
  "Erythroid3",
  "NMP",
  "Rostral_neurectoderm",
  "Caudal_neurectoderm",
  "Neural_crest",
  "Forebrain_Midbrain_Hindbrain",
  "Spinal_cord",
  "Surface_ectoderm"
  # "Visceral_endoderm",
  # "ExE_endoderm",
  # "ExE_ectoderm",
  # "Parietal_endoderm"
)

########################
## Load cell metadata ##
########################

sample_metadata <- fread(io$metadata) %>%
  .[(pass_atacQC==TRUE | pass_rnaQC==TRUE) & doublet_call==FALSE] %>%
  .[sample%in%opts$samples & celltype.predicted%in%opts$celltypes]

#########################
## Load LSI dimensions ##
#########################

io$lsi.atac <- paste0(io$basedir,"/results/atac/archR/dimensionality_reduction/PeakMatrix/all_cells/E7.5_rep1-E7.5_rep2-E8.0_rep1-E8.0_rep2-E8.5_rep1-E8.5_rep2_lsi_features50000_ndims30.txt.gz")
lsi.mtx <- fread(io$lsi.atac) %>% matrix.please %>% .[sample_metadata[pass_atacQC==T,cell]] %>% t

##############################
## Load RNA expression data ##
##############################

sce <- load_SingleCellExperiment(io$rna.sce, normalise = TRUE, cells = sample_metadata$cell)
dim(sce)

# Add sample metadata to the colData of the SingleCellExperiment
colData(sce) <- sample_metadata %>% as.data.frame %>% tibble::column_to_rownames("cell") %>%
  .[colnames(sce),] %>% DataFrame()

# Feature selection
decomp <- modelGeneVar(sce)
decomp <- decomp[decomp$mean > 0.001,]
hvgs <- rownames(decomp)[decomp$p.value <= 0.10]
sce_filt <- sce[hvgs,]

rna.mtx <- assay(sce_filt,"logcounts") %>% as.matrix
dim(rna.mtx)

###########################
## Prepare data for MOFA ##
###########################

cells <- intersect(colnames(rna.mtx),colnames(lsi.mtx))

rna.mtx <- rna.mtx[,cells]
lsi.mtx <- lsi.mtx[,cells]

########################
## Create MOFA object ##
########################

MOFAobject <- create_mofa_from_matrix(
  list(
  "RNA" = rna.mtx,
  "ATAC" = lsi.mtx
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
model_opts$num_factors <- 30
model_opts$spikeslab_weights <- FALSE
model_opts$ard_weights <- FALSE

# Training options
train_opts <- get_default_training_options(MOFAobject)
train_opts$convergence_mode <- "medium"
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