suppressPackageStartupMessages(library(MOFA2))

#####################
## Define settings ##
#####################

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))


# I/O
io$outdir <- file.path(io$basedir, "results_new/rna_atac/mofa/all_cells")
io$mofa.output <- file.path(io$basedir, "results_new/rna_atac/mofa/all_cells/mofa_allcells_pca_lsi.rds")

# Options
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
  "Surface_ectoderm",
  "Visceral_endoderm",
  "ExE_endoderm",
  "ExE_ectoderm",
  "Parietal_endoderm"
)

########################
## Load cell metadata ##
########################

sample_metadata <- fread(io$metadata) %>%
  .[(pass_atacQC==TRUE & pass_rnaQC==TRUE) & doublet_call==FALSE] %>%
  .[celltype.mapped_mnn%in%opts$celltypes]

opts$rna.cells <- sample_metadata[pass_rnaQC==TRUE,cell]
opts$atac.cells <- sample_metadata[pass_atacQC==TRUE,cell]

###############
## Load data ##
###############

# RNA (PCA)
io$pca.rna <- paste0(io$basedir,"/results_new/rna/dimensionality_reduction/sce/pca_features2500_pcs50.txt.gz")
tmp <- fread(io$pca.rna)
opts$rna.cells <- intersect(tmp$cell,opts$rna.cells)
rna.mtx <- tmp %>% matrix.please %>% .[opts$rna.cells,] %>% t
rm(tmp)

# ATAC (LSI)
io$lsi.atac <- paste0(io$basedir,"/results_new/atac/archR/dimensionality_reduction/lsi_PeakMatrix_nfeatures25000_ndims50.txt.gz")
tmp <- fread(io$lsi.atac)
opts$atac.cells <- intersect(tmp$cell,opts$atac.cells)
atac.mtx <- tmp %>% matrix.please %>% .[opts$atac.cells,] %>% t
rm(tmp)

###########################
## Prepare data for MOFA ##
###########################

rna.mtx <- .augment_matrix(rna.mtx, unique(rna.cells,atac.cells))
atac.mtx <- .augment_matrix(atac.mtx, unique(rna.cells,atac.cells))

########################
## Create MOFA object ##
########################

MOFAobject <- create_mofa_from_matrix(
  list(
  "RNA" = rna.mtx,
  "ATAC" = atac.mtx
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
model_opts$num_factors <- 25
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


fwrite(metadata.to.mofa, file.path(io$outdir,"sample_metadata.txt.gz"), quote=F, sep="\t", na="NA")
