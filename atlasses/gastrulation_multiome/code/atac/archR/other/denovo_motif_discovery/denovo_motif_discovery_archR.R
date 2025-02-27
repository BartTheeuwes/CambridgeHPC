# https://www.archrproject.com/bookdown/chromvar-deviatons-enrichment-with-archr.html
# https://greenleaflab.github.io/chromVAR/articles/Articles/Applications.html#motif-kmer-similarity

library(SummarizedExperiment)
library(chromVAR)
library(BSgenome.Mmusculus.UCSC.mm10)

########################
## Load ArchR Project ##
########################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/atac/archR/load_archR_project.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/atac/archR/load_archR_project.R")
} else {
  stop("Computer not recognised")
}

#####################
## Define settings ##
#####################

# I/O
# io$metadata <- paste0(io$basedir,"/processed/atac/archR/sample_metadata_after_archR.txt.gz")
io$outdir <- paste0(io$basedir,"/results/atac/archR/denovo_motif_discovery")

# Options
opts$samples <- c(
  "E7.5_rep1",
  "E7.5_rep2",
  "E8.5_rep1",
  "E8.5_rep2"
)

#####################
## Update metadata ##
#####################

sample_metadata <- fread(io$metadata) %>%
  .[pass_atacQC==TRUE] %>%
  .[sample%in%opts$samples]
stopifnot(sample_metadata$archR_cell %in% rownames(ArchRProject))

# subset celltypes with sufficient number of cells
opts$min.cells <- 100
sample_metadata <- sample_metadata %>%
  .[,N:=.N,by=c("celltype.predicted")] %>% .[N>opts$min.cells] %>% .[,N:=NULL]
table(sample_metadata$celltype.predicted)

##################
## Subset ArchR ##
##################

ArchRProject.filt <- ArchRProject[sample_metadata$archR_cell,]
table(getCellColData(ArchRProject.filt,"Sample")[[1]])
table(getCellColData(ArchRProject.filt,"celltype.predicted")[[1]])

##################
## Add Peak Set ##
##################

# Load Peak Set
# peaks.dt <- fread(io$archR.peakSet.bed) %>%
#   setnames(c("chr","start","end")) %>%
#   .[,id:=sprintf("%s:%s-%s",chr,start,end)]
# peaks.granges <- makeGRangesFromDataFrame(peaks.dt, keep.extra.columns = T)
# ArchRProject.filt <- addPeakSet(ArchRProject.filt, peaks.granges)

###################
## Sanity checks ##
###################

# Motif annotations require a PeakSet
peaks.se <- getPeakSet(ArchRProject.filt) %>% head(n=1000)
mcols(peaks.se)$id <- paste(seqnames(peaks.se), mcols(peaks.se)$idx, sep="_")
stopifnot(sum(duplicated(mcols(peaks.se)$id))==0)

metadata(ArchRProject.filt@peakSet)$bgdPeaks <- file.path(io$archR.directory, "Background-Peaks.rds")
background_peaks <- getBgdPeaks(ArchRProject.filt, method = "chromVAR") %>% head(n=1000)

###################
## Define k-mers ##
###################

kmer_ix <- matchKmers(6, peaks.se, genome = BSgenome.Mmusculus.UCSC.mm10)

# mtx <- getMatrixFromProject(ArchRProject.filt, useMatrix="PeakMatrix")
# assayNames(mtx) <- "counts"
# rowRanges(mtx)$id <- paste(seqnames(rowRanges(mtx)), mcols(rowRanges(mtx))$idx, sep="_")
# stopifnot(sum(duplicated(mcols(rowRanges(mtx))$id))==0)
# mtx <- mtx[rowRanges(mtx)$id %in% peaks.se$id,]
# saveRDS(mtx, "~/Downloads/mtx.rds")
mtx <- readRDS("~/Downloads/mtx.rds")

mtx <- addGCBias(mtx, genome = BSgenome.Mmusculus.UCSC.mm10)


# reteurns chromVARDeviations-class, which inherits from SummarizedExperiment, and has two assays: deviations and deviation scores.
kmer_dev <- computeDeviations(mtx, kmer_ix)

# Save summarised experiment object
# saveRDS(deviations.se, sprintf("%s/deviations_summarized_experiment.rds",io$outdir))

###########################################################
## Create de novo motifs from k-mers based on deviations ##
###########################################################

# returns list with 
# - (1) motifs: de novo motif matrices
# - (2) seed: seed kmer for de novo motif
de_novos <- assembleKmers(kmer_dev, threshold = 1.5, p = 0.01, progress = F)

saveRDS(de_novos, sprintf("%s/deviations_summarized_experiment.rds",io$outdir))

#############################
## Compare to known motifs ##
#############################

# dist_to_known <- pwmDistance(de_novos, motifs)
# closest_match1 <- which.min(dist_to_known$dist[1,])
# dist_to_known$strand[1,closest_match1]

################
## START TEST ##
################

# R.utils::sourceDirectory("/Users/ricard/git/chromVAR/R/", verbose=T, modifiedOnly=FALSE)
# library(BiocParallel)
# register(SerialParam())
# 
# object <- mtx
# annotations <- kmer_ix
# background_peaks <- getBackgroundPeaks(object)
# expectation <- computeExpectations(object)
# 
# object <- counts_check(object)
# annotations <- matches_check(annotations)
# peak_indices <- convert_to_ix_list(annotationMatches(annotations))
# dev <- compute_deviations_core(counts(object), peak_indices, background_peaks, expectation, colData = colData(object), rowData = colData(annotations))
  
##############
## END TEST ##
##############
