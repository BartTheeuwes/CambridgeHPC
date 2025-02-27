suppressPackageStartupMessages(library(argparse))

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--groupA',    type="character",    help='group A')
p$add_argument('--groupB',    type="character",    help='group B')
p$add_argument('--test',      type="character",    help='Statistical test')
p$add_argument('--motif_annotation',    type="character",    help='Motif annotation')
p$add_argument('--outfile',   type="character",    help='Output file')
args <- p$parse_args(commandArgs(TRUE))

## START TEST
# args <- list()
# args$motif_annotation <- "Motif_cisbp"
# args$test <- "wilcoxon"
# args$groupA <- "Epiblast"
# args$groupB <- "Cardiomyocytes"
# args$outfile <- tempfile()
## END TEST

# Sanity checks
stopifnot(args$test%in%c("ttest","wilcoxon"))

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

########################
## Load cell metadata ##
########################

# io$metadata <- "/Users/ricard/data/gastrulation_multiome_10x/processed/atac/archR/sample_metadata_after_archR.txt.gz"
sample_metadata <- fread(io$metadata) %>%
  .[pass_atacQC==TRUE] %>%
  .[sample%in%opts$samples & celltype.predicted%in%c(args$groupA,args$groupB)]
table(sample_metadata$celltype.predicted)

#########################
## Subset ArchR object ##
#########################

ArchRProject.filt <- ArchRProject[sample_metadata$cell]

# add celltype.predicted to ArchR's CellColData
foo <- sample_metadata %>% 
  .[cell%in%rownames(ArchRProject.filt)] %>% setkey(cell) %>% .[rownames(ArchRProject.filt)] %>%
  as.data.frame() %>% tibble::column_to_rownames("cell")
stopifnot(all(foo$TSSEnrichment_atac == getCellColData(ArchRProject.filt, "TSSEnrichment")[[1]]))
ArchRProject.filt <- addCellColData(
  ArchRProject.filt,
  data = foo$celltype.predicted, 
  name = "celltype.predicted",
  cells = rownames(foo),
  force = TRUE
)

table(ArchRProject.filt$celltype.predicted)

###################
## Sanity checks ##
###################

stopifnot(args$motif_annotation%in%names(ArchRProject.filt@peakAnnotation))
stopifnot(c(args$groupA,args$groupB)%in%ArchRProject.filt$celltype.predicted)

#######################
## Differential test ##
#######################

markerTest <- getMarkerFeatures(
  ArchRProject.filt, 
  # useMatrix = sprintf("DeviationMatrix_%s",args$motif_annotation),
  useMatrix = sprintf("DeviationMatrix_%s",args$motif_annotation),
  groupBy = "celltype.predicted",
  testMethod = args$test,
  useSeqnames = "z",
  useGroups = args$groupA,
  bgdGroups = args$groupB
)

####################
## Prepare output ##
####################

chromvar.diff.dt <- data.table(
  # Log2FC = round(assay(markerTest,"Log2FC")[,1],3),
  name =  rowData(markerTest)$name,
  MeanDiff = round(assay(markerTest,"MeanDiff")[,1],3),
  FDR = assay(markerTest,"FDR")[,1]
  # AUC = round(assay(markerTest,"AUC")[,1],3)
) %>% setorder(FDR)

##########
## Save ##
##########

fwrite(chromvar.diff.dt, args$outfile, sep="\t")


################
## START TEST ##
################

# R.utils::sourceDirectory("/Users/ricard/git/ArchR/R/", verbose=T, modifiedOnly=FALSE)

# ArchRProj = ArchRProject.filt
# groupBy = "celltype.predicted"
# useGroups = args$groupA
# bgdGroups = args$groupB
# useMatrix = "PeakMatrix"
# bias = c("TSSEnrichment","log10(nFrags)")
# normBy = NULL
# testMethod = "wilcoxon"
# maxCells = 500
# scaleTo = 10^4
# threads = 1
# k = 100
# bufferRatio = 0.8
# binarize = FALSE
# useSeqnames = NULL
# verbose = TRUE
# logFile = createLogFile("getMarkerFeatures")

##############
## END TEST ##
##############