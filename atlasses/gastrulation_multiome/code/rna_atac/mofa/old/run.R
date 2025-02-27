
suppressPackageStartupMessages(library(MOFA2))
suppressPackageStartupMessages(library(argparse))

here::i_am("rna_atac/mofa/run.R")

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--python',          type="character",   help='Python binary')
p$add_argument('--nfactors',           type="integer",    default=30,                  help='Number of MOFA factors')
p$add_argument('--seed',            type="integer",    default=42,                  help='Random seed')
p$add_argument('--cores',            type="integer",    default=1,                  help='Number of cores')
p$add_argument('--outfile',          type="character",                               help='Output file')
p$add_argument('--test',  action="store_true",  help='Test mode?')

args <- p$parse_args(commandArgs(TRUE))



##################################
## Set up reticulate connection ##
##################################

if (!is.null(args$python)) {
  reticulate::use_python(args$python, required=T)
}

########################
## Create MOFA object ##
########################

MOFAobject <- create_mofa_from_matrix(list("RNA" = rna.mtx, "ATAC" = atac_tfidf.mtx))
MOFAobject

# Clear memory
rm(list=c("rna.mtx","atac_tfidf.mtx","rna.sce","atac.se"))
gc()


####################
## Define options ##
####################

# Data options
data_opts <- get_default_data_options(MOFAobject)
data_opts$use_float32 <- TRUE
# data_opts$scale_views <- FALSE

# Model options
model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- args$nfactors
model_opts$spikeslab_weights <- FALSE
# model_opts$ard_weights <- FALSE

# Training options
train_opts <- get_default_training_options(MOFAobject)
train_opts$convergence_mode <- "fast"
train_opts$seed <- args$seed

if (args$test) {
  train_opts$maxiter <- 5  
}

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

saveRDS(MOFAobject, args$outfile)