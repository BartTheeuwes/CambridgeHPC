
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(MOFA2))


#####################
## Define settings ##
#####################

# load default settings
source("/Users/ricard/gastrulation_multiome_10x/settings.R")

# I/O
# io$metadata <- paste0(io$basedir, "/sample_metadata.txt.gz")
# io$outdir <- paste0(io$basedir, "/results/atac/signac/dimensionality_reduction")
io$mofa.output <- paste0(io$basedir, "/results/rna/mofa/mofa_model_blood_rna.rds")

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

########################
## Load cell metadata ##
########################

sample_metadata <- fread(io$metadata) %>%
  .[pass_rnaQC==TRUE & doublet_call==FALSE] %>%
  .[celltype.predicted%in%opts$celltypes]
table(sample_metadata$celltype.predicted)

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

########################
## Create MOFA object ##
########################

# Sanity check
MOFAobject <- create_mofa_from_SingleCellExperiment(sce_filt, extract_metadata = TRUE)
views_names(MOFAobject) <- "RNA"

MOFAobject

####################
## Define options ##
####################

# Data options
data_opts <- get_default_data_options(MOFAobject)
# data_opts$use_float32 <- FALSE

# Model options
model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 15
model_opts$spikeslab_weights <- FALSE
model_opts$ard_weights <- FALSE

# Training options
train_opts <- get_default_training_options(MOFAobject)
# train_opts$convergence_mode <- "fast"
train_opts$maxiter <- 150

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
