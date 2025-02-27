here::i_am("atac/archR/chromvar_chip/run_chromvar_chip.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

suppressPackageStartupMessages(library(GenomicRanges))

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--metadata',  type="character",              help='Cell metadata') 
p$add_argument('--atac_peak_matrix',  type="character",              help='ATAC Peak matrix')
p$add_argument('--virtual_chip_mtx',             type="character",        help='Virtual ChIP matrix')
p$add_argument('--peak_metadata',  type="character",              help='Peak metadata')
p$add_argument('--motifmatcher',  type="character",              help='Motif annotation') 
p$add_argument('--motif_annotation',  type="character",              help='Motif annotation') 
p$add_argument('--min_number_peaks',     type="integer",    default=30,    help='Minimum number of peaks per TF')
p$add_argument('--min_chip_score',     type="integer",    default=0.15,    help='Minimum ChIP score')
p$add_argument('--ignore_negative_values',  action="store_true",  help='Ignote negative ChIP-seq scores')
p$add_argument('--outdir',          type="character",                help='Output directory')
p$add_argument('--test',  action="store_true",  help='Test mode?')
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

## START TEST ##
io$basedir <- file.path(io$basedir,"test")
args <- list()
args$metadata <- file.path(io$basedir,"results/atac/archR/celltype_assignment/sample_metadata_after_celltype_assignment.txt.gz")
args$virtual_chip_mtx <- file.path(io$basedir,"results/rna_atac/virtual_chipseq/pseudobulk/CISBP/virtual_chip_matrix.rds")
args$atac_peak_matrix <- file.path(io$basedir,"processed/atac/archR/Matrices/PeakMatrix_summarized_experiment.rds")
args$atac_peak_metadata <- file.path(io$basedir,"processed/atac/archR/PeakCalls/peak_metadata.tsv.gz")
args$background_peaks <- file.path(io$basedir,"processed/atac/archR/PeakCalls/peak_metadata.tsv.gz")
args$motif_annotation <- "CISBP"
args$min_number_peaks <- 30
args$min_chip_score <- 0.15
args$motifmatcher <- file.path(io$basedir,"results/rna_atac/virtual_chipseq/pseudobulk/CISBP/motifmatchr_virtual_chip.rds")
args$ignore_negative_values <- TRUE
args$test <- TRUE
args$outdir <- file.path(io$basedir,"results/atac/archR/chromvar_chip/cells")
## END TEST ##

# I/O

###################
## Load metadata ##
###################

sample_metadata <- fread(args$metadata) %>%
  .[pass_atacQC==TRUE & doublet_call==FALSE]

if (args$test) {
  sample_metadata <- sample_metadata %>% head(n=250)
}

##########################
## Load ATAC PeakMatrix ##
##########################

atac_peakMatrix.se <- readRDS(args$atac_peak_matrix)

# Load peak metadata
peak_metadata.dt <- fread(args$atac_peak_metadata) %>%
  .[,idx:=sprintf("%s:%s-%s",chr,start,end)] %>%
  .[idx%in%rownames(atac_peakMatrix.se)] %>% setkey(idx) %>%
  .[rownames(atac_peakMatrix.se)]

stopifnot(peak_metadata.dt$idx==rownames(atac_peakMatrix.se))

# temporary sanity check
if (any(!c("start","strand")%in%colnames(rowData(atac_peakMatrix.se)))) {
  tmp <- rownames(atac_peakMatrix.se) %>% strsplit(":") %>% map_chr(2)
  rowData(atac_peakMatrix.se)$start <- tmp %>% strsplit("-") %>% map_chr(1)
  rowData(atac_peakMatrix.se)$end <- tmp %>% strsplit("-") %>% map_chr(2)
}

###################################
## Load virtual ChIP-seq library ##
###################################

# Load virtual chip matrix
virtual_chip.mtx <- readRDS(args$virtual_chip_mtx)

# Define TFs
if (args$test) {
  virtual_chip.mtx <- virtual_chip.mtx[,c("FOXA2","MIXL1","TAL1")]
}

# Filter peaks with no TFs
virtual_chip.mtx <- virtual_chip.mtx[rowSums(virtual_chip.mtx!=0)>0,]

###############################
## Load motifmatcher results ##
###############################

# source(here::here("load_motifmatchR.R"))
motifmatcher.se <- readRDS(args$motifmatcher)

############################
## Select overlapping TFs ##
############################

TFs <- intersect(colnames(motifmatcher.se),colnames(virtual_chip.mtx))

motifmatcher.se <- motifmatcher.se[,TFs]
virtual_chip.mtx <- virtual_chip.mtx[,TFs]

########################################################################
## Update motifmatcher results according to the minimum binding score ##
########################################################################

if (args$ignore_negative_values){
  print("Ignoring negative TF binding scores...")
  virtual_chip.mtx[virtual_chip.mtx<0] <- 0
}

for (i in TFs) {
  # peaks <- which(virtual_chip.mtx[,i]>=args$min_chip_score) %>% names
  peaks <- which(abs(virtual_chip.mtx[,i])>=args$min_chip_score) %>% names
  
  assay(motifmatcher.se)[,i] <- FALSE
  assay(motifmatcher.se)[,i][rownames(motifmatcher.se)%in%peaks] <- TRUE
  
  cat(sprintf("%s (%s/%s):\n(-) %s (+) %s\tbefore filtering\n(-) %s (+) %s\tafter filtering\n", i,match(i,TFs),length(TFs),sum(virtual_chip.mtx[,i]<0),sum(virtual_chip.mtx[,i]>0), sum(virtual_chip.mtx[peaks,i]<0),sum(virtual_chip.mtx[peaks,i]>0)))
}


###########################
## Load background peaks ##
###########################

bgdPeaks.se <- readRDS(io$archR.bgdPeaks)
tmp <- rowRanges(bgdPeaks.se)
rownames(bgdPeaks.se) <- sprintf("%s:%s-%s",seqnames(tmp), start(tmp), end(tmp))
bgdPeaks.se <- bgdPeaks.se[rownames(atac_peakMatrix.se),]


###################
## Sanity checks ##
###################

stopifnot(rownames(atac_peakMatrix.se)==rownames(motifmatcher.se))
stopifnot(rownames(bgdPeaks.se)==rownames(atac_peakMatrix.se))

#####################################################################
## Use chromVAR default implementation (slow and memory intensive) ##
#####################################################################

# # prepare data for chromvar
# assayNames(atac_peakMatrix.se) <- "counts"
# 
# # Compute deviations
# chromvar_deviations.se <- chromVAR::computeDeviations(
#   object = atac_peakMatrix.se,
#   annotations = motifmatcher.se,
#   background_peaks = assay(bgdPeaks.se)
# )
# 
# # Save
# if (!args$test) {
#   saveRDS(chromvar_deviations.se, sprintf("%s/chromVAR_deviations_%s_chromvar_chip.rds",args$outdir,args$motif_annotation))
# }


#########################################
## Use ArchR's chromVAR implementation ##
#########################################

source(here::here("atac/archR/chromvar/utils.R"))

featureDF <- data.frame(
  rowSums = rowSums(assay(atac_peakMatrix.se)),
  start = start(rowRanges(atac_peakMatrix.se)),
  end = end(rowRanges(atac_peakMatrix.se)),
  GC = peak_metadata.dt$GC
)

archr_deviations.se <- .customDeviations(
  countsMatrix = assay(atac_peakMatrix.se),
  annotationsMatrix = as(assay(motifmatcher.se,"motifMatches"),"dgCMatrix"),
  backgroudPeaks = assay(bgdPeaks.se),
  expectation = featureDF$rowSums/sum(featureDF$rowSums),
  prefix = "",
  out = c("deviations", "z"),
  threads = 1,
  verbose = TRUE
)

# Save
saveRDS(archr_deviations.se, file.path(args$outdir,sprintf("chromVAR_chip_deviations_%s_archr.rds",args$motif_annotation)))

##########
## Test ##
##########

# hist(assay(chromvar_deviations.se,"deviations"))

# countsMatrix = assay(atac_peakMatrix.se)
# annotationsMatrix = as(assay(motifmatcher.se),"dgCMatrix")
# backgroudPeaks = assay(bgdPeaks.se)
# expectation = featureDF$rowSums/sum(featureDF$rowSums)
# prefix = ""
# out = c("z", "deviations")
# threads = 1
# verbose = TRUE


##########
## Test ##
##########

# chromvar_deviations_nofilt.se <- readRDS("/Users/ricard/data/gastrulation_multiome_10x/results/atac/archR/chromvar/deviations_summarized_experiment_Motif_cisbp.rds")[,TFs]
# rownames(chromvar_deviations_nofilt.se) <- rowData(chromvar_deviations_nofilt.se)$name %>% toupper %>% stringr::str_split(.,"_") %>% map_chr(1)
# chromvar_deviations_nofilt.se <- chromvar_deviations_nofilt.se[TFs,]
# 
# to.plot <- assay(chromvar_deviations_nofilt.se,"z") %>% as.matrix %>% as.data.table(keep.rownames = T) %>% setnames("rn","TF") %>% 
#   melt(id.vars="TF", variable.name="cell") %>%
#   merge(sample_metadata[,c("cell","celltype.predicted")])
# 
# p <- ggboxplot(to.plot, x="celltype.predicted", y="value", fill="celltype.predicted", outlier.shape=NA) +
#   scale_fill_manual(values=opts$celltype.colors) +
#   facet_wrap(~TF, ncol=1) +
#   geom_hline(yintercept=0, linetype="dashed") +
#   # coord_cartesian(ylim=c(-1,1)) +
#   labs(x="", y="") +
#   theme(
#     legend.position = "none",
#     axis.text.y = element_text(size=rel(0.7)),
#     # axis.text.x = element_text(size=rel(0.7), hjust=1, angle=40)
#     axis.text.x = element_blank(),
#     legend.title = element_blank(),
#     axis.ticks.x = element_blank()
#   )
#   
# print(p)