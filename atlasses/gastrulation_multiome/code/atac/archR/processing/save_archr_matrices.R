here::i_am("atac/archR/processing/save_archr_matrices.R")

#####################
## Define settings ##
#####################

source(here::here("settings.R"))
source(here::here("utils.R"))
  
# I/O

###################
## Load metadata ##
###################

sample_metadata <- fread(io$metadata) %>%
  .[pass_atacQC==TRUE & doublet_call==FALSE]

#########################
## Load ATAC-based PCA ##
#########################

pca.atac <- fread(io$pca.atac) %>% matrix.please

########################
## Load ArchR Project ##
########################

source(here::here("atac/archR/load_archR_project.R"))

# Subset
ArchRProject.filt <- ArchRProject[sample_metadata$cell,]

getAvailableMatrices(ArchRProject.filt)

################
## PeakMatrix ##
################

atac.peakMatrix.se <- getMatrixFromProject(ArchRProject.filt, binarize = TRUE, useMatrix = "PeakMatrix")

# Define peak names
row_ranges.dt <- rowRanges(atac.peakMatrix.se) %>% as.data.table %>% 
  setnames("seqnames","chr") %>%
  .[,c("chr","start","end")] %>%
  .[,idx:=sprintf("%s:%s-%s",chr,start,end)]
rownames(atac.peakMatrix.se) <- row_ranges.dt$idx

# Load peak metadata
peak_metadata.dt <- fread(io$archR.peak.metadata) %>%
  .[,idx:=sprintf("%s:%s-%s",chr,start,end)]
atac.peakMatrix.se <- atac.peakMatrix.se[rownames(atac.peakMatrix.se) %in% peak_metadata.dt$idx]

# Denoise

# Save
io$outfile <- paste0(io$archR.directory,"/PeakCalls/PeakMatrix_summarized_experiment.rds")
saveRDS(atac.peakMatrix.se, io$outfile)

#####################
## GeneScoreMatrix ##
#####################

# atac.GeneScoreMatrix.se <- getMatrixFromProject(ArchRProject.filt, binarize = TRUE, useMatrix = "GeneScoreMatrix_nodistal")
atac.GeneScoreMatrix.se <- getMatrixFromProject(ArchRProject.filt, binarize = FALSE, useMatrix = "GeneScoreMatrix_nodistal")
rownames(atac.GeneScoreMatrix.se) <- rowData(atac.GeneScoreMatrix.se)$name
stopifnot(sum(duplicated(rownames(atac.GeneScoreMatrix.se)))==0)

# Filter genes
atac.GeneScoreMatrix.se <- atac.GeneScoreMatrix.se[grep("^Rik|Rik$|^mt-|^Rps-|^Rpl-|^Gm|^Mir|^Olfr",rownames(atac.GeneScoreMatrix.se),invert=T),]

# Remove lowly variable genes
cells <- intersect(rownames(pca.atac),colnames(atac.GeneScoreMatrix.se))
tmp <- sparseMatrixStats::rowVars(assay(atac.GeneScoreMatrix.se))
names(tmp) <- rownames(atac.GeneScoreMatrix.se)
atac.GeneScoreMatrix.se <- atac.GeneScoreMatrix.se[tmp>0.05,]

# Denoise
atac.GeneScoreMatrix.mtx <- smoother_aggregate_nearest_nb(mat=as.matrix(assay(atac.GeneScoreMatrix.se[,cells])), D=pdist(pca.atac[cells,]), k=25)
assay(atac.GeneScoreMatrix.se) <- atac.GeneScoreMatrix.mtx

# Save
# io$outfile <- paste0(io$archR.directory,"/GeneScoreMatrix_no_distal_summarized_experiment_denoised.rds")
io$outfile <- paste0(io$archR.directory,"/GeneScoreMatrix_no_distal_summarized_experiment.rds")
saveRDS(atac.GeneScoreMatrix.se, io$outfile)

##############
## chromVAR ##
##############


