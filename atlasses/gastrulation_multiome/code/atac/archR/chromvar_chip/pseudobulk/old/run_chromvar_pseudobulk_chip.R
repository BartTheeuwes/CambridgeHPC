here::i_am("atac/archR/chromvar_chip/pseudobulk/run_chromvar_pseudobulk_chip.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(chromVAR))

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--motif_annotation',  type="character",              help='Motif annotation') 
p$add_argument('--atac_peak_matrix',  type="character",              help='ATAC Peak matrix (pseudobulk)') 
p$add_argument('--motifmatcher',  type="character",              help='Motif annotation') 
p$add_argument('--peak_metadata',  type="character",              help='') 
p$add_argument('--background_peaks',  type="character",              help='') 
p$add_argument('--min_number_peaks',     type="integer",    default=30,    help='Minimum number of peaks per TF')
p$add_argument('--min_chip_score',     type="double",    default=0.15,    help='Minimum ChIP score')
p$add_argument('--ignore_negative_values',  action="store_true",  help='Ignote negative ChIP-seq scores')
p$add_argument('--outdir',          type="character",                help='Output directory')
p$add_argument('--test',  action="store_true",  help='Test mode?')
args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# io$basedir <- file.path(io$basedir,"test")
# args <- list()
# args$motif_annotation <- "CISBP"
# args$atac_peak_matrix <- file.path(io$basedir,"results/atac/archR/pseudobulk/celltype/PeakMatrix/pseudobulk_PeakMatrix_summarized_experiment.rds")
# args$motifmatcher <- file.path(io$basedir,sprintf("results/rna_atac/virtual_chipseq/pseudobulk/%s/motifmatchr_virtual_chip.rds",args$motif_annotation))
# args$background_peaks <- file.path(io$basedir,"/processed/atac/archR/Background-Peaks.rds")
# args$peak_metadata <- file.path(io$basedir,"/processed/atac/archR/PeakCalls/peak_metadata.tsv.gz")
# args$min_number_peaks <- 30
# args$min_chip_score <- 0.10
# args$outdir <- file.path(io$basedir,"results/atac/archR/chromvar_chip/pseudobulk")
# args$ignore_negative_values <- TRUE
# args$test <- TRUE
## END TEST ##

dir.create(args$outdir, showWarnings=F, recursive=T)

################################
## Load pseudobulk PeakMatrix ##
################################

atac_peakMatrix_pseudobulk.se <- readRDS(args$atac_peak_matrix)#[,opts$celltypes]

# Load peak metadata
peak_metadata.dt <- fread(args$peak_metadata) %>% 
  .[,idx:=sprintf("%s:%s-%s",chr,start,end)] %>%
  .[,c("idx","score","GC")] %>% 
  setkey(idx) %>% .[rownames(atac_peakMatrix_pseudobulk.se)]
stopifnot(peak_metadata.dt$idx==rownames(atac_peakMatrix_pseudobulk.se))

###############################
## Load motifmatcher results ##
###############################

print("Loading motifmatcher results...")

motifmatcher.se <- readRDS(args$motifmatcher)

stopifnot(c("motifMatches","VirtualChipScores")%in%names(assays(motifmatcher.se)))

###################################################################
## Update motifmatcher results using the virtual ChIP-seq scores ##
###################################################################

if (args$ignore_negative_values){
  print(sprintf("Number of matches before filtering negative TF binding values: %d",sum(assay(motifmatcher.se,"motifMatches"))))
  # assays(motifmatcher.se) <- assays(motifmatcher.se)[""]
  assay(motifmatcher.se,"motifMatches")[assay(motifmatcher.se,"VirtualChipScores")<0] <- F
  print(sprintf("Number of matches after filtering negative TF binding values: %d",sum(assay(motifmatcher.se,"motifMatches"))))
}

print(sprintf("Number of matches before filtering based on minimum ChIP-seq score: %d",sum(assay(motifmatcher.se,"motifMatches"))))
# assays(motifmatcher.se) <- assays(motifmatcher.se)[""]
assay(motifmatcher.se,"motifMatches")[abs(assay(motifmatcher.se,"VirtualChipScores"))<=args$min_chip_score] <- F
print(sprintf("Number of matches after filtering based on minimum ChIP-seq score: %d",sum(assay(motifmatcher.se,"motifMatches"))))

assays(motifmatcher.se) <- assays(motifmatcher.se)["motifMatches"]

################
## Filter TFs ##
################

# Filter TFs with too few peaks
TFs <- which(colSums(assay(motifmatcher.se,"motifMatches"))>=args$min_number_peaks) %>% names
TFs.removed <- which(colSums(assay(motifmatcher.se,"motifMatches"))<args$min_number_peaks) %>% names

cat(sprintf("%s TFs removed because they don't have enough binding sites: %s", length(TFs.removed), paste(TFs.removed, collapse=" ")))

if (args$test) {
  TFs <- c("FOXA2","MIXL1","GATA1","EOMES","BCL11B","DLX2","FOXC1")
}
motifmatcher.se <- motifmatcher.se[,TFs]

###########################
## Load background peaks ##
###########################

print("Loading background peaks...")

bgdPeaks.se <- readRDS(args$background_peaks)
tmp <- rowRanges(bgdPeaks.se)
rownames(bgdPeaks.se) <- sprintf("%s:%s-%s",seqnames(tmp), start(tmp), end(tmp))
stopifnot(sort(rownames(bgdPeaks.se))==sort(rownames(atac_peakMatrix_pseudobulk.se)))
bgdPeaks.se <- bgdPeaks.se[rownames(atac_peakMatrix_pseudobulk.se),]

###################
## Sanity checks ##
###################

stopifnot(rownames(atac_peakMatrix_pseudobulk.se)==rownames(motifmatcher.se))

#########################################
## Use chromVAR default implementation ##
#########################################

print("Running chromVAR's default implementation...")

# prepare data for chromvar
assayNames(atac_peakMatrix_pseudobulk.se) <- "counts"

# Compute deviations
chromvar_deviations_chromvar.se <- computeDeviations(
  object = atac_peakMatrix_pseudobulk.se,
  annotations = motifmatcher.se,
  background_peaks = assay(bgdPeaks.se)
)

# Save
saveRDS(chromvar_deviations_chromvar.se, file.path(args$outdir,sprintf("chromVAR_chip_%s_chromvar.rds",args$motif_annotation)))

#########################################
## Use ArchR's chromVAR implementation ##
#########################################

print("Running ArchR's chromVAR implementation...")

source(here::here("atac/archR/chromvar/utils.R"))

featureDF <- data.frame(
  rowSums = rowSums(assay(atac_peakMatrix_pseudobulk.se)),
  start = rowData(atac_peakMatrix_pseudobulk.se)$start,
  end = rowData(atac_peakMatrix_pseudobulk.se)$end,
  GC = peak_metadata.dt$GC
)

chromvar_deviations_archr.se <- .customDeviations(
  countsMatrix = assay(atac_peakMatrix_pseudobulk.se),
  annotationsMatrix = as(assay(motifmatcher.se),"dgCMatrix"),
  backgroudPeaks = assay(bgdPeaks.se),
  expectation = featureDF$rowSums/sum(featureDF$rowSums),
  prefix = "",
  out = c("deviations", "z"),
  threads = 1,
  verbose = TRUE
)

# Save
saveRDS(chromvar_deviations_archr.se, file.path(args$outdir,sprintf("chromVAR_chip_%s_archr.rds",args$motif_annotation)))

###################################################################
## Compare default chromVAR with ArchR's chromVAR implementation ##
###################################################################

# chromvar.deviations.dt_A <- assay(chromvar_deviations_chromvar.se,"z") %>% as.matrix %>% as.data.frame %>%
#   as.data.table(keep.rownames = T) %>% setnames("rn","motif") %>%
#   melt(id.vars="motif", variable.name="celltype") %>%
#   .[,class:=factor("chromvar")]
# 
# chromvar.deviations.dt_B <- assay(chromvar_deviations_archr.se) %>% as.matrix %>% as.data.frame %>%
#   as.data.table(keep.rownames = T) %>% setnames("rn","motif") %>%
#   melt(id.vars="motif", variable.name="celltype") %>%
#   .[,class:=factor("archr")]
# 
# to.plot <- rbind(chromvar.deviations.dt_A, chromvar.deviations.dt_B)
# 
# opts$celltype.colors <- opts$celltype.colors[names(opts$celltype.colors)%in%unique(to.plot$celltype)]
# 
# motifs.to.plot <- unique(to.plot$motif)
# # motifs.to.plot <- c("FOXA2","GATA2","TAL1")
# 
# for (i in motifs.to.plot) {
#   
#   to.plot.barplot <- to.plot[motif==i] %>% 
#     .[,celltype:=factor(celltype,levels=names(opts$celltype.colors))]
#   
#   p1 <- ggbarplot(to.plot.barplot, x="celltype", y="value", fill="celltype") +
#     scale_fill_manual(values=opts$celltype.colors) +
#     facet_wrap(~class, scales="free_y") +
#     geom_hline(yintercept=0, linetype="dashed") +
#     guides(x = guide_axis(angle = 90)) +
#     labs(x="", y=sprintf("%s chromVAR z-score",i)) +
#     theme(
#       legend.position = "none",
#       # axis.text.x = element_text(color="black", angle=40, hjust=1, size=rel(0.75)),
#       axis.text = element_text(color="black", size=rel(0.6)),
#       axis.title = element_text(color="black", size=rel(0.8))
#     )
#   
#   to.plot.scatter <- to.plot.barplot %>% dcast(motif+celltype~class, value.var="value")
#   p2 <- ggscatter(to.plot.scatter, x="chromvar", y="archr", fill="celltype", shape=21, stroke=0.1, size=3.5, 
#                   add="reg.line", add.params = list(color="black", fill="lightgray"), conf.int=TRUE) +
#     stat_cor(method = "pearson", label.x.npc = "middle", label.y.npc = "bottom") +
#     scale_fill_manual(values=opts$celltype.colors) +
#     labs(x="z-score (chromVAR)", y="z-score (ArchR)") +
#     theme(
#       legend.position = "none",
#       axis.text = element_text(size=rel(0.7)),
#       axis.title = element_text(size=rel(0.85))
#     )
#   
#   p <- cowplot::plot_grid(plotlist=list(p1,p2), nrow=1, rel_widths = c(2/3,1/3))
#   
#   pdf(file.path(args$outdir,sprintf("pdf/%s_chromvar_vs_archr.pdf",i)), width=13, height=5)
#   print(p)
#   dev.off()
# }

################
## START TEST ##
################

# library(BiocParallel)
# R.utils::sourceDirectory("/Users/ricard/git/chromVAR/R/", verbose=T, modifiedOnly=FALSE)
# R.utils::sourceDirectory("/homes/ricard/git/chromVAR/R/", verbose=T, modifiedOnly=FALSE)

# object <- atac_peakMatrix_pseudobulk.se
# background_peaks = assay(bgdPeaks.se)
# annotations = motifmatcher.se
# rowData = NULL
# colData = NULL
# 
# counts_mat <- counts(atac_peakMatrix_pseudobulk.se)
# expectation = computeExpectations(counts_mat)
# peak_indices <- convert_to_ix_list(annotationMatches(annotations))
# 
# peak_set = peak_indices[[2]]

