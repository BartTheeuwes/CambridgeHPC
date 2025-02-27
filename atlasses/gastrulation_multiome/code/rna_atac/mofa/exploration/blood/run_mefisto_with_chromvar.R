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

io$mofa.output <- paste0(io$basedir, "/results/rna_atac/mofa/blood/mofa_model.rds")

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
  "Mixed_mesoderm",
  "Haematoendothelial_progenitors",
  "Endothelium",
  "Blood_progenitors_1",
  "Blood_progenitors_2",
  "Erythroid1"
  # "Erythroid2",
  # "Erythroid3"
)

########################
## Load cell metadata ##
########################

# io$metadata <- paste0(io$basedir, "/sample_metadata.txt.gz")
sample_metadata <- fread(io$metadata) %>%
  .[pass_atacQC==TRUE & pass_rnaQC==TRUE & doublet_call==FALSE] %>%
  .[sample%in%opts$samples & celltype.mapped%in%opts$celltypes] %>%
  .[,celltype.mapped:=factor(celltype.mapped,levels=opts$celltypes)]

# Select N cells per cell type
sample_metadata <- sample_metadata %>% 
  split(.$celltype.mapped) %>% map(~ head(.,n=50)) %>%
  rbindlist
table(sample_metadata$celltype.mapped)

#####################
## Load trajectory ##
#####################

io$pseudotime <- "/Users/ricard/data/gastrulation_multiome_10x/results/rna/trajectories/haematoendothelium_trajectory/haematoendothelium_trajectory.txt.gz"
trajectory.dt <- fread(io$pseudotime) %>% 
  .[,V1:=minmax.normalisation(V1)] %>%
  .[,V2:=minmax.normalisation(V2)] %>%
  .[cell%in%sample_metadata$cell]

ggplot(trajectory.dt, aes(x=V1, y=V2)) +
  geom_point() +
  theme_classic()

##########################
## Subset ArchR project ##
##########################

stopifnot(sample_metadata$archR_cell %in% rownames(ArchRProject))
ArchRProject.filt <- ArchRProject[sample_metadata$archR_cell]

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

################################
## Feature selection chromVAR ##
################################

chromvar.features <- sort(apply(atac.deviation.mtx,1,var), decreasing = T) %>% head(n=250) %>% names

atac.deviation.mtx <- atac.deviation.mtx[chromvar.features,]
dim(atac.deviation.mtx)

##############################
## Load RNA expression data ##
##############################

io$rna.sce <- "/Users/ricard/data/gastrulation_multiome_10x/processed/rna/old/SingleCellExperiment.rds"
sce <- load_SingleCellExperiment(io$rna.sce, normalise = TRUE, cells = sample_metadata$cell)
dim(sce)

# Add sample metadata to the colData of the SingleCellExperiment
colData(sce) <- sample_metadata %>% as.data.frame %>% tibble::column_to_rownames("cell") %>%
  .[colnames(sce),] %>% DataFrame()

# Feature selection
sce <- sce[!grepl("Rik|Gm",rownames(sce)),]

decomp <- modelGeneVar(sce)
decomp <- decomp[decomp$mean > 0.01,]
hvgs <- rownames(decomp)[decomp$p.value <= 0.05]
sce_filt <- sce[hvgs,]

rna.mtx <- assay(sce_filt,"logcounts")
dim(rna.mtx)

##############################
## Prepare data for MEFISTO ##
##############################

cells <- Reduce("intersect",list(colnames(rna.mtx),colnames(atac.deviation.mtx),trajectory.dt$cell))

rna.mtx <- rna.mtx[,cells]
atac.deviation.mtx <- atac.deviation.mtx[,cells]
pseudotime.values <- trajectory.dt %>% matrix.please %>% t %>% .[,cells,drop=F]

########################
## Create MOFA object ##
########################

MOFAobject <- create_mofa_from_matrix(
  list(
    "RNA" = rna.mtx
    # "ATAC_chromVAR" = atac.deviation.mtx
  )
)
MOFAobject

# Set pseudotime as a covariate
all(samples_names(MOFAobject)[[1]]==names(pseudotime.values))
MOFAobject <- set_covariates(MOFAobject, covariates=pseudotime.values)

####################
## Define options ##
####################

# Data options
data_opts <- get_default_data_options(MOFAobject)
data_opts$use_float32 <- TRUE

# Model options
model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 5
model_opts$spikeslab_weights <- FALSE

# Training options
train_opts <- get_default_training_options(MOFAobject)
# train_opts$maxiter <- 5

mefisto_opts <- get_default_mefisto_options(MOFAobject)
# mefisto_opts$sparseGP <- TRUE
# mefisto_opts$frac_inducing <- 0.50

#########################
## Prepare MOFA object ##
#########################

MOFAobject <- prepare_mofa(
  MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts,
  mefisto_options = mefisto_opts
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
