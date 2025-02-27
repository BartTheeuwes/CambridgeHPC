
suppressPackageStartupMessages(library(argparse))

## Initialize argument parser ##
p <- ArgumentParser(description='')
p$add_argument('--groupA',    type="character",    help='group A')
p$add_argument('--groupB',    type="character",    help='group B')
p$add_argument('--test',      type="character",    help='Statistical test')
p$add_argument('--matrix',    type="character",    help='Matrix to use, see getAvailableMatrices')
# p$add_argument('--test_mode', action="store_true", help='Test mode? subset number of cells')
p$add_argument('--outfile',   type="character",    help='Output file')
args <- p$parse_args(commandArgs(TRUE))

## START TEST
# args <- list()
# args$matrix <- "GeneScoreMatrix"  # PeakMatrix
# args$test <- "wilcoxon"
# args$groupA <- "Erythroid"
# args$groupB <- "Cardiomyocytes"
# # args$outfile <- sprintf("%s/results/differential/archR/%s_%s_vs_%s.tsv.gz",io$basedir,args$matrix,args$celltypeA,args$celltypeB)
# args$outfile <- tempfile()
# args$test_mode <- FALSE
## END TEST

# Sanity checks
stopifnot(args$test%in%c("binomial","ttest","wilcoxon"))

########################
## Load ArchR project ##
########################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/public_datasets/Pijuan-Sala_2020/archR/load_archR_project.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/public_datasets/Pijuan-Sala_2020/archR/load_archR_project.R")
} else {
  stop("Computer not recognised")
}

stopifnot(args$matrix%in%getAvailableMatrices(ArchRProject))

#######################
## Differential test ##
#######################

markerTest <- getMarkerFeatures(
  ArchRProject, 
  useMatrix = args$matrix,
  groupBy = "celltype",
  testMethod = args$test,
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = args$groupA,
  bgdGroups = args$groupB
)

##########
## Plot ##
##########

# plotMarkers(markerTest, name = "Erythroid", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "MA")
# plotMarkers(markerTest, name = "Erythroid", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "Volcano")

##########
## Save ##
##########

cols <- c("seqnames","start","end")
if ("name" %in% colnames(rowData(markerTest))) {
  cols <- c(cols,"name")
}
dt.1 <- rowData(markerTest)[,cols] %>% as.data.table %>%
  setnames("seqnames","chr")# %>%
  # .[,id:=sprintf("%s:%s-%s",chr,start,end)]

dt.2 <- data.table(
  Log2FC = round(assay(markerTest,"Log2FC")[,1],3),
  MeanDiff = round(assay(markerTest,"MeanDiff")[,1],3),
  FDR = assay(markerTest,"FDR")[,1],
  AUC = round(assay(markerTest,"AUC")[,1],3)
)
dt <- cbind(dt.1,dt.2) %>% setorder(-AUC)

fwrite(dt, args$outfile, sep="\t")
