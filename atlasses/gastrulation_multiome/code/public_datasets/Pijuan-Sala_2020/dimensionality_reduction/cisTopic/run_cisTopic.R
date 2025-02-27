library(Seurat)
library(cisTopic)

#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/public_datasets/Pijuan-Sala_2020/settings.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/public_datasets/Pijuan-Sala_2020/settings.R")
} else {
  stop("Computer not recognised")
}

################
## Define I/O ##
################

io$outdir <- paste0(io$basedir,"/results/cisTopic")

####################
## Define options ##
####################

opts$test <- FALSE

#####################
## Update metadata ##
#####################

# sample_metadata <- sample_metadata
  # .[batch%in%opts$batches & celltype.mapped%in%opts$celltypes]

if (opts$test) sample_metadata <- head(sample_metadata,n=100)

################################
## Load data in Seurat format ##
################################

# Load Seurat object
seurat <- readRDS(io$seurat)#[,sample_metadata$cell]
dim(seurat)

# Add metadata to Seurat object
# seurat@meta.data <- sample_metadata %>% as.data.frame %>% 
#   tibble::column_to_rownames("cell") %>%
#   .[cell%in%colnames(seurat)] %>% .[colnames(seurat),]

#######################################
## Load data in sparse matrix format ##
#######################################

# Load sparse matrix
# m <- Matrix::readMM(io$matrix)
# barcodes <- fread(io$barcodes, header=F)[[1]]
# features <- fread(io$features, header=F)[[1]]
# if (opts$test) {
#   n <- 500
#   m <- m[1:n,1:n]
#   features <- features[1:n]
#   barcodes <- barcodes[1:n]
# }
# features <- features %>% stringr::str_replace(.,"_",":") %>% stringr::str_replace(.,"_","-")
# rownames(m) <- features
# colnames(m) <- barcodes

##############
## CisTopic ##
##############

# Create cisTopicObject
m <- seurat@assays[["peaks"]]@counts
rownames(m) <- rownames(m) %>% stringr::str_replace(.,"-",":")
cisTopicObject <- createcisTopicObject(m, project.name='PijuanSala2020', keepCountsMatrix = FALSE)

# Define number of topics
topic <- c(20,30,40,50,60,70)
# topic <- 10

# Run Latent Dirichlet Allocation
cisTopicObject <- runCGSModels(cisTopicObject, topic, burnin = 125, iterations = 250, nCores = 3, seed = 123)
# cisTopicObject <- selectModel(cisTopicObject)

# Run Latent Dirichlet Allocation with WarpLDA
# cisTopicObject <- runWarpLDAModels(cisTopicObject, topic, iterations = 500, nCores=2, seed = 42)

##########
## Save ##
##########

saveRDS(cisTopicObject, paste0(io$outdir,"/cistopic.rds"))
