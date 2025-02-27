suppressPackageStartupMessages(library(cisTopic))
suppressPackageStartupMessages(library(ggpubr))

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

io$metadata <- paste0(io$basedir, "/sample_metadata.txt.gz")
io$outdir <- paste0(io$basedir,"/results/atac/cistopic/blood")

####################
## Define options ##
####################

opts$samples <- c(
  "E7.5_rep1",
  "E7.5_rep2",
  "E8.5_rep1",
  "E8.5_rep2"
)

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

opts$test <- FALSE

opts$npeaks <- 5e3

########################
## Load cell metadata ##
########################

sample_metadata <- fread(io$metadata) %>%
  .[pass_atacQC==TRUE] %>%
  .[sample%in%opts$samples & celltype.predicted%in%opts$celltypes]

if (opts$test) sample_metadata <- head(sample_metadata,n=100)

##################
## Subset ArchR ##
##################

stopifnot(sample_metadata$archR_cell %in% rownames(ArchRProject))
ArchRProject.filt <- ArchRProject[sample_metadata$archR_cell]

######################
## Load peak matrix ##
######################

# Get Peak Matrix from ArchR
atac.peak.se <- getMatrixFromProject(ArchRProject.filt, useMatrix="PeakMatrix")
dim(atac.peak.se)

# Define peak names
peak_names <- rowRanges(atac.peak.se) %>% as.data.table %>% .[,id:=sprintf("%s:%s-%s",seqnames,start,end)] %>% .$id
rownames(atac.peak.se) <- peak_names

#######################
## Feature selection ##
#######################

# Load peak variability estimates
peak.variability.dt <- fread(io$archR.peak.variability)

# Define highly variable peaks
peaks <- peak.variability.dt %>% 
  setorder(-variance_pseudobulk) %>% 
  head(n=opts$npeaks) %>% 
  .[,peak:=stringr::str_replace(peak,"_",":")] %>%
  .[,peak:=stringr::str_replace(peak,"_","-")] %>%
  .$peak

# Subset
atac.peak.se.filt <- atac.peak.se[peaks,]
dim(atac.peak.se.filt)
rm(atac.peak.se)

# Extract sparse mtx
peak.mtx <- assay(atac.peak.se.filt,"PeakMatrix")

##############
## CisTopic ##
##############

# Create cisTopicObject
cisTopicObject <- createcisTopicObject(peak.mtx, project.name='GastrulationMultiome10x', keepCountsMatrix = FALSE)

# Define number of topics
if (opts$test) {
  ntopics <- c(10,20)
} else {
  ntopics <- c(30,40,50,60,70,80,90)
}

# Run Latent Dirichlet Allocation
# cisTopicObject <- runCGSModels(cisTopicObject, topic, burnin = 125, iterations = 250, nCores = 2, seed = 123)

# Run Latent Dirichlet Allocation with WarpLDA
cisTopicObject <- runWarpLDAModels(
  object = cisTopicObject, 
  topic = ntopics, 
  iterations = 500, 
  nCores = 2, 
  seed = 42,
  returnType = "allModels"
)

#####################
## Model selection ##
#####################

# Returns a cisTopic object (when the input is a cisTopic object) with the selected model 
# stored in object@selected.model, and the log likelihoods of the models in object@log.lik
cisTopicObject.best_model <- selectModel(
  cisTopicObject, 
  type = "maximum", 
  keepBinaryMatrix = F, 
  keepModels = F
)

# Create data.table 
lik <- cisTopicObject@models %>% map_dbl(function(x) x$log.likelihoods[2,ncol(x$log.likelihoods)])
lik.dt <- data.table(
  ntopics = names(lik),
  lik = lik
)

p <- ggline(lik.dt, x="ntopics", y="lik")

pdf(paste0(io$outdir,"/model_selection.pdf"), width=7, height=3)
print(p)
dev.off()

##########
## Save ##
##########

saveRDS(cisTopicObject, paste0(io$outdir,"/cistopic_warpLDA_allModels.rds"))
saveRDS(cisTopicObject.best_model, paste0(io$outdir,"/cistopic_warpLDA_bestModel.rds"))

