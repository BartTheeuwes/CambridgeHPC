
####################
## Load libraries ##
####################
library(SingleCellExperiment)
library(scran)
#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/public_datasets/settings.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/public_datasets/settings.R")
} else {
  stop("Computer not recognised")
}

# Define I/O
io$outdir <- paste0(io$basedir,"/results/GRN")
io$sce <- paste0(io$basedir,'/processed/rna/pseudobulk/SingleCellExperiment.rds')
io$tf <- paste0(io$basedir,'/mm_tfs.txt')
io$expr <- paste0(io$outdir,'/expr.csv')
io$path.to.pyscenic <- '/home/lijingyu/anaconda3/bin/pyscenic'
io$adj <- paste0(io$outdir,'/adj.csv')
io$adj_cor <- paste0(io$outdir,'/adj_cor.csv')
# Define options
opts$test <- TRUE

#####################
## Update metadata ##
#####################

if (opts$test) sample_metadata <- head(sample_metadata,n=100)

###############
## Load data ##
###############

# Load singlecellexpresssion object
sce <- readRDS(io$sce)

#################
## Filter data ##
#################

dec <- modelGeneVar(sce)
plot(dec$mean, dec$total, xlab="Mean log-expression", ylab="Variance")
curve(metadata(dec)$trend(x), col="blue", add=TRUE)
top.hvgs <- getTopHVGs(dec, n=3000)
tf= read.table(io$tf)
tf=tf$V1[tf$V1 %in% rownames(sce)]
genes=union(tf,top.hvgs)
sce<- sce[genes,]


#################
## Run grnboost##
#################

write.csv(t(sce@assays@data$logcounts),file = io$expr)
system(paste(io$path.to.pyscenic,'grn',io$expr,io$tf,'-o',io$adj))
system(paste(io$path.to.pyscenic,'add_cor',io$adj,io$expr,'-o',io$adj_cor))
regulon<-read.csv(io$adj_cor)
