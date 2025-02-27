
suppressPackageStartupMessages(library(argparse))

here::i_am("atac/archR/differential/archr_differential_accessibility.R")


## Initialize argument parser ##
p <- ArgumentParser(description='')
p$add_argument('--matrix',    type="character",    help='Matrix to use, see getAvailableMatrices')
p$add_argument('--groupA',    type="character",    help='group A')
p$add_argument('--groupB',    type="character",    help='group B')
p$add_argument('--celltype_label',    type="character",    help='Celltype label')
p$add_argument('--test',      type="character",    help='Statistical test')
p$add_argument('--outfile',   type="character",    help='Output file')
args <- p$parse_args(commandArgs(TRUE))

## START TEST
# args <- list()
# args$matrix <- "GeneScoreMatrix_TSS" # "PeakMatrix"
# args$groupA <- "Epiblast"
# args$groupB <- "Cardiomyocytes"
# args$celltype_label <- "celltype.mapped_mnn"
# args$test <- "wilcoxon"
# args$outfile <- tempfile()
## END TEST

# Sanity checks
stopifnot(args$test%in%c("binomial","ttest","wilcoxon"))

###################
## Load settings ##
###################

source(here::here("settings.R"))
source(here::here("utils.R"))

########################
## Load cell metadata ##
########################

sample_metadata <- fread(io$metadata) %>%
  .[pass_atacQC==TRUE & doublet_call==FALSE]

stopifnot(args$celltype_label%in%colnames(sample_metadata))

sample_metadata <- sample_metadata %>%
  .[,celltype:=eval(as.name(args$celltype_label))] %>%
  .[celltype%in%c(args$groupA,args$groupB)]

########################
## Load ArchR project ##
########################

source(here::here("atac/archR/load_archR_project.R"))

# Subset ArchR object
ArchRProject.filt <- ArchRProject[sample_metadata$cell]

# add celltype to ArchR's CellColData
foo <- sample_metadata %>% 
  .[cell%in%rownames(ArchRProject.filt)] %>% setkey(cell) %>% .[rownames(ArchRProject.filt)] %>%
  as.data.frame() %>% tibble::column_to_rownames("cell")
stopifnot(all(foo$nFrags_atac == getCellColData(ArchRProject.filt,"nFrags_atac")[[1]]))

ArchRProject.filt <- addCellColData(
  ArchRProject.filt,
  data = foo$celltype, 
  name = "celltype",
  cells = rownames(foo),
  force = TRUE
)

###################
## Sanity checks ##
###################

stopifnot(args$matrix%in%getAvailableMatrices(ArchRProject.filt))
stopifnot(c(args$groupA,args$groupB)%in%ArchRProject.filt$celltype)
stopifnot(sample_metadata$cell %in% rownames(ArchRProject.filt))

#######################
## Differential test ##
#######################

markerTest <- getMarkerFeatures(
  ArchRProject.filt, 
  useMatrix = args$matrix,
  groupBy = "celltype",
  testMethod = args$test,
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = args$groupA,
  bgdGroups = args$groupB
)

####################
## Prepare output ##
####################

dt.1 <- rowData(markerTest) %>% as.data.table %>%
  setnames("seqnames","chr") %>%
  .[,idx:=sprintf("%s:%s-%s",chr,start,end)] %>%
  .[,c("name","idx")]
  # .[,c("chr","start","end"):=NULL]

dt.2 <- data.table(
  # Log2FC = round(assay(markerTest,"Log2FC")[,1],3),
  MeanDiff = round(assay(markerTest,"MeanDiff")[,1],3),
  FDR = assay(markerTest,"FDR")[,1],
  AUC = round(assay(markerTest,"AUC")[,1],3)
)

dt <- cbind(dt.1,dt.2) %>% setorder(-AUC)

##########
## Save ##
##########

fwrite(dt, args$outfile, sep="\t")


################
## START TEST ##
################

# R.utils::sourceDirectory("/Users/ricard/git/ArchR/R/", verbose=T, modifiedOnly=FALSE)

# ArchRProj = ArchRProject.filt
# groupBy = "celltype"
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