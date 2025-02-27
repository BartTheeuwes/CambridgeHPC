here::i_am("atac/archR/dimensionality_reduction/cells/archR_dimensionality_reduction.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

suppressPackageStartupMessages(library(ArchR))

######################
## Define arguments ##
######################

## START TEST ##
args <- list()
args$metadata <- "/rds/project/rds-SDzz0CATGms/users/bt392/09_Eomes_invitro_blood/results/atac/archR/qc/sample_metadata_after_qc.txt.gz"
args$nfeatures <- 25000
args$matrix <- "PeakMatrix"
args$feature_stats = '/rds/project/rds-SDzz0CATGms/users/bt392/09_Eomes_invitro_blood/results/atac/archR/feature_stats/PeakMatrix/PeakMatrix_Clusters_PeakMatrix2_stats.txt.gz'
args$peakmatrix = '/rds/project/rds-SDzz0CATGms/users/bt392/09_Eomes_invitro_blood/processed/atac/archR/Matrices/PeakMatrix_summarized_experiment.rds'
args$vars_to_regress <- c('TSSEnrichment_atac')
args$regression_out = file.path(io$basedir, 'results/atac/archR/regression/')

dir.create(args$regression_out, recursive=TRUE, showWarnings =FALSE)
## END TEST ##

# Load ArchR project 

source(here::here("atac/archR/load_archR_project.R"))

# get top features
stats = fread(args$feature_stats)
features = stats[order(-var_pseudobulk)] %>% head(., args$nfeatures) %>% .$feature

# load peakmatrix
PeakMatrix = readRDS(args$peakmatrix)


# Load sample metadata 
sample_metadata <- fread(args$metadata) %>%
  .[,log_nFrags_atac:=log10(nFrags_atac)]

sample_metadata = sample_metadata[match(colnames(PeakMatrix), cell),]
summary(sample_metadata$cell == colnames(PeakMatrix))

# TF-IDF normalisation
peak.mtx = PeakMatrix[features,]@assays@data$PeakMatrix
peak_norm.mtx = tfidf(peak.mtx)

a = Sys.time()
# Regress out variables 
peak_norm_regressed.mtx <- RegressOutMatrix(
    mtx = peak_norm.mtx,
    covariates = sample_metadata[,args$vars_to_regress,with=F]
  )

b = Sys.time()
print('time for regression:')
b - a

a = Sys.time()
fwrite(peak_norm_regressed.mtx, file.path(args$regression_out,"peak_norm_regressed_mtx.txt.gz"))
b = Sys.time()
print('time for saving:')
b - a
