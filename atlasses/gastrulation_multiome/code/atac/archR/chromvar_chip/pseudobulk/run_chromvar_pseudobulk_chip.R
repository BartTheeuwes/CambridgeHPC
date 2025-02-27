suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(chromVAR))
suppressPackageStartupMessages(library(argparse))

here::i_am("atac/archR/chromvar/pseudobulk/run_chromvar_pseudobulk_chip.R")

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--atac_peak_matrix',  type="character",              help='ATAC Peak matrix (pseudobulk)') 
p$add_argument('--virtual_chip_mtx',             type="character",        help='Virtual ChIP matrix')
# p$add_argument('--motif2gene',  type="character",              help='Motif annotation')
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

source(here::here("settings.R"))
source(here::here("utils.R"))

## START TEST ##
# args <- list()
# args$atac_peak_matrix <- file.path(io$basedir,"results_new/atac/archR/pseudobulk/celltype.mapped_mnn/pseudobulk_PeakMatrix_summarized_experiment.rds")
# args$virtual_chip_mtx <- file.path(io$basedir,"results_new/rna_atac/virtual_chipseq/CISBP/virtual_chip.mtx")
# args$motif_annotation <- "CISBP"
# args$min_number_peaks <- 30
# args$min_chip_score <- 0.15
# args$motifmatcher <- file.path(io$basedir,"results_new/rna_atac/virtual_chipseq/CISBP/motifmatchr_virtual_chip.rds")
# args$outdir <- file.path(io$basedir,"results_new/atac/archR/chromvar_chip/pseudobulk")
# args$ignore_negative_values <- TRUE
# args$test <- TRUE
## END TEST ##

################################
## Load pseudobulk PeakMatrix ##
################################

atac_peakMatrix_pseudobulk.se <- readRDS(args$atac_peak_matrix)#[,opts$celltypes]

# Load peak metadata
peak_metadata.dt <- fread(io$archR.peak.metadata) %>% 
  .[,idx:=sprintf("%s:%s-%s",chr,start,end)]

# Define peak names
# peak_names <- rowData(atac_peakMatrix_pseudobulk.se) %>% as.data.table %>% .[,idx:=sprintf("%s:%s-%s",seqnames,start,end)] %>% .$id
# rownames(atac_peakMatrix_pseudobulk.se) <- peak_names

# Subset peaks
# atac_peakMatrix_pseudobulk.se <- atac_peakMatrix_pseudobulk.se[rownames(atac_peakMatrix_pseudobulk.se) %in% peak_metadata.dt$idx]
# peak_metadata.dt <- peak_metadata.dt %>%
#   .[idx%in%rownames(atac_peakMatrix_pseudobulk.se)] %>% setkey(idx) %>%
#   .[rownames(atac_peakMatrix_pseudobulk.se)]

###################################
## Load virtual ChIP-seq library ##
###################################

# virtual_chip.dt <- opts$TFs %>% map(function(i) {
#   fread(sprintf("%s/%s.bed.gz",io$virtual_chip.dir,i)) %>%
#     setnames(c("chr","start","end","score")) %>%
#     .[,idx:=sprintf("%s:%s-%s",chr,start,end)] %>%
#     .[score>=args$min_chip_score] %>%
#     .[,tf:=i]
# }) %>% rbindlist
virtual_chip.mtx <- readRDS(args$virtual_chip_mtx)

# Define TFs
if (args$test) {
  virtual_chip.mtx <- virtual_chip.mtx[,c("FOXA2","MIXL1","TAL1")]
}

# Fix peak names (temporary)
if (!all(grepl(":",rownames(virtual_chip.mtx)))) {
  rownames(virtual_chip.mtx) <- sub("-",":",rownames(virtual_chip.mtx))
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

# Filter peaks with no motif matches (DON'T, WE NEED ALL PEAKS TO MATCH BACKGROUND PEAKS)
# motifmatcher.se <- motifmatcher.se[rowSums(assay(motifmatcher.se))>0,]
# atac.peakMatrix.se <- atac.peakMatrix.se[rownames(motifmatcher.se),]

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

################
## Filter TFs ##
################

# Filter TFs with too few peaks
TFs <- which(colSums(assay(motifmatcher.se,"motifMatches"))>=args$min_number_peaks) %>% names
TFs.removed <- which(colSums(assay(motifmatcher.se,"motifMatches"))<args$min_number_peaks) %>% names

cat(sprintf("%s TFs removed because they don't have enough binding sites: %s", length(TFs.removed), paste(TFs.removed, collapse=" ")))

motifmatcher.se <- motifmatcher.se[,TFs]
virtual_chip.mtx <- virtual_chip.mtx[,TFs]

################
## Rename TFs ##
################

# source(here::here("atac/archR/load_motif_annotation.R"))
# 
# motifs <- intersect(colnames(motifmatcher.se),motif2gene.dt$motif)
# genes <- intersect(colnames(virtual_chip.mtx),motif2gene.dt$gene)
# 
# motif2gene_filt.dt <- motif2gene.dt[motif%in%motifs & gene%in%genes]
# motifs <- motif2gene_filt.dt$motif
# genes <- motif2gene_filt.dt$gene
# 
# tmp <- genes; names(tmp) <- motifs
# 
# stopifnot(motif2gene_filt.dt$motif%in%colnames(motifmatcher.se))
# stopifnot(motif2gene_filt.dt$gene%in%colnames(virtual_chip.mtx))
# 
# motifmatcher.se <- motifmatcher.se[,motifs]
# colnames(motifmatcher.se) <- tmp[colnames(motifmatcher.se)]
# virtual_chip.mtx <- virtual_chip.mtx[,genes]
# 
# stopifnot(colnames(motifmatcher.se)==colnames(virtual_chip.mtx))

###########################
## Load background peaks ##
###########################

bgdPeaks.se <- readRDS(io$archR.bgdPeaks)
tmp <- rowRanges(bgdPeaks.se)
rownames(bgdPeaks.se) <- sprintf("%s:%s-%s",seqnames(tmp), start(tmp), end(tmp))
bgdPeaks.se <- bgdPeaks.se[rownames(atac_peakMatrix_pseudobulk.se),]

###################
## Sanity checks ##
###################

stopifnot(rownames(atac_peakMatrix_pseudobulk.se)==rownames(motifmatcher.se))

#########################################
## Use chromVAR default implementation ##
#########################################

# prepare data for chromvar
assayNames(atac_peakMatrix_pseudobulk.se) <- "counts"

# Compute deviations
chromvar_deviations_chromvar.se <- computeDeviations(
  object = atac_peakMatrix_pseudobulk.se,
  annotations = motifmatcher.se,
  background_peaks = assay(bgdPeaks.se)
)

# Save
saveRDS(chromvar_deviations_chromvar.se, file.path(args$outdir,sprintf("chromVAR_deviations_%s_chromvar_chip.rds",args$motif_annotation)))

#########################################
## Use ArchR's chromVAR implementation ##
#########################################

source(here::here("atac/archR/chromvar/utils.R"))

featureDF <- data.frame(
  rowSums = rowSums(assay(atac_peakMatrix_pseudobulk.se)),
  start = rowData(atac_peakMatrix_pseudobulk.se)$start,
  end = rowData(atac_peakMatrix_pseudobulk.se)$end,
  GC = peak_metadata.dt$GC
)

chromvar_deviations_archr.se <- .customDeviations(
  countsMatrix = assay(atac_peakMatrix_pseudobulk.se),
  annotationsMatrix = as(assay(motifmatcher.se,"motifMatches"),"dgCMatrix"),
  backgroudPeaks = assay(bgdPeaks.se),
  expectation = featureDF$rowSums/sum(featureDF$rowSums),
  prefix = "",
  out = c("deviations", "z"),
  threads = 1,
  verbose = TRUE
)

# Save
saveRDS(chromvar_deviations_archr.se, file.path(args$outdir,sprintf("chromVAR_deviations_%s_archr_chip.rds",args$motif_annotation)))

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
