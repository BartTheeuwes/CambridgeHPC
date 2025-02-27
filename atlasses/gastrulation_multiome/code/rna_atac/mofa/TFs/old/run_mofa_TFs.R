suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(MOFA2))

########################
## Define settings ##
########################

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


io$mofa.output <- paste0(io$basedir, "/results/rna_atac/mofa/TFs/mofa_TFs_forebrain.rds")

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
  "Visceral_endoderm"
  # "ExE_endoderm"
  # "ExE_ectoderm",
  # "Parietal_endoderm"
)
opts$celltypes <- "Forebrain_Midbrain_Hindbrain"

########################
## Load cell metadata ##
########################

sample_metadata <- fread(io$metadata) %>%
  .[(pass_atacQC==TRUE & pass_rnaQC==TRUE) & doublet_call==FALSE] %>%
  # .[(pass_atacQC==TRUE | pass_rnaQC==TRUE) & doublet_call==FALSE] %>%
  .[sample%in%opts$samples & celltype.predicted%in%opts$celltypes]

#########################
## Load RNA expression ##
#########################

source("/Users/ricard/gastrulation_multiome_10x/rna_atac/load_rna_atac_single_cells.R")

rm(list=c("atac.peakMatrix.se","rna.sce")); gc()

#######################
## Feature selection ##
#######################

# RNA
tmp <- sparseMatrixStats::rowVars(logcounts(rna.sce.tf))
rna.sce.tf <- rna.sce.tf[tmp>0.01,]

# rna.mtx <- assay(sce_filt,"logcounts")
# dim(rna.mtx)

TFs <- intersect(rownames(rna.sce.tf),rownames(atac.chromvar.se))
rna.sce.tf <- rna.sce.tf[TFs,]
atac.chromvar.se <- atac.chromvar.se[TFs,]

###########################
## Prepare data for MOFA ##
###########################

# cells <- intersect(colnames(rna.sce.tf),colnames(atac.chromvar.se))
# rna.sce.tf <- rna.sce.tf[,cells]
# atac.chromvar.se <- atac.chromvar.se[,cells]

########################
## Create MOFA object ##
########################

MOFAobject <- create_mofa_from_matrix(list(
  "TF_RNA" = as.matrix(logcounts(rna.sce.tf)),
  "TF_chromVAR" = assay(atac.chromvar.se,"z")
  ))
MOFAobject

####################
## Define options ##
####################

# Data options
data_opts <- get_default_data_options(MOFAobject)
data_opts$use_float32 <- TRUE
# data_opts$scale_views <- FALSE

# Model options
model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 3# 30
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
