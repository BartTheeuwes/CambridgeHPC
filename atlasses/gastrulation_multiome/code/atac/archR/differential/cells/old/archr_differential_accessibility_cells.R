here::i_am("atac/archR/differential/cells/archr_differential_accessibility_cells.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

suppressPackageStartupMessages(library(ArchR))

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--archr_directory',    type="character",    help='ArchR directory')
p$add_argument('--metadata',    type="character",    help='Cell metadata')
p$add_argument('--matrix',    type="character",    help='Matrix to use, see getAvailableMatrices')
p$add_argument('--groupA',    type="character",    help='group A')
p$add_argument('--celltypes',    type="character",    default="all", nargs="+", help='Celltypes to use')
p$add_argument('--samples',    type="character",    default="all", nargs="+", help='Samples to use')
p$add_argument('--groupB',    type="character",    help='group B')
p$add_argument('--group_variable',    type="character",    help='Group label')
p$add_argument('--statistical_test',      type="character",    help='Statistical test')
p$add_argument('--outfile',   type="character",    help='Output file')
args <- p$parse_args(commandArgs(TRUE))

## START TEST
# args <- list()
# args$archr_directory <- file.path(io$basedir,"test/processed/atac/archR")
# args$matrix <- "PeakMatrix" 
# # args$groupA <- "Epiblast"
# # args$groupB <- "Cardiomyocytes"
# # args$group_variable <- "celltype"
# args$metadata <- file.path(io$basedir,"results/atac/archR/qc/cells_metadata.dt_after_qc.txt.gz")
# args$groupA <- "WT"
# args$groupB <- "T_KO"
# args$samples <- opts$samples
# args$group_variable <- "genotype"
# args$celltypes <- "NMP"
# args$statistical_test <- "wilcoxon"
# args$outfile <- tempfile()
## END TEST

#####################
## Parse arguments ##
#####################

# Sanity checks
stopifnot(args$statistical_test%in%c("binomial","ttest","wilcoxon"))

# Define cell types
if (args$celltypes[1]=="all") {
  args$celltypes <- opts$celltypes
} else {
  # stopifnot(args$celltypes%in%c(opts$celltypes,"Erythroid","Blood_progenitors"))
  stopifnot(args$celltypes%in%opts$celltypes)
}

# Define samples
if (args$samples[1]=="all") {
  args$samples <- opts$samples
} else {
  stopifnot(args$samples%in%opts$samples)
}

# Define groups
opts$groups <- c(args$groupA,args$groupB)

########################
## Load cell metadata ##
########################

cells_metadata.dt <- fread(args$metadata) %>%
  .[pass_atacQC==TRUE & pass_rnaQC==TRUE & doublet_call==FALSE & celltype%in%args$celltypes & sample%in%args$samples] %>%
  .[,celltype_genotype:=sprintf("%s_%s",celltype,genotype)]

stopifnot(args$group_variable%in%colnames(cells_metadata.dt))

cells_metadata.dt <- cells_metadata.dt %>%
  .[,group:=eval(as.name(args$group_variable))] %>%
  .[group%in%opts$groups]

print(table(cells_metadata.dt$group))

########################
## Load ArchR project ##
########################

# source(here::here("atac/archR/load_archR_project.R"))

setwd(args$archr_directory)

addArchRGenome("mm10")
addArchRThreads(threads = 1)

ArchRProject <- loadArchRProject(args$archr_directory)[cells_metadata.dt$cell]

# Sanity checks
stopifnot(args$matrix%in%getAvailableMatrices(ArchRProject))

# add group to ArchR's CellColData
foo <- cells_metadata.dt %>% 
  .[cell%in%rownames(ArchRProject)] %>% setkey(cell) %>% .[rownames(ArchRProject)] %>%
  as.data.frame() %>% tibble::column_to_rownames("cell")

# sanity check
stopifnot(all(foo$nFrags_atac == getCellColData(ArchRProject,"nFrags_atac")[[1]]))

ArchRProject <- addCellColData(
  ArchRProject,
  data = foo$group, 
  name = "group",
  cells = rownames(foo),
  force = TRUE
)

#######################
## Differential test ##
#######################

markerTest <- getMarkerFeatures(
  ArchRProject, 
  useMatrix = args$matrix,
  groupBy = "group",
  testMethod = args$statistical_test,
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
  .[,c("chr","start","end"):=NULL]

dt.2 <- data.table(
  # Log2FC = round(assay(markerTest,"Log2FC")[,1],3),
  MeanDiff = -assay(markerTest,"MeanDiff")[,1] %>% round(3),
  FDR = assay(markerTest,"FDR")[,1] %>% format(digits=3),
  AUC = assay(markerTest,"AUC")[,1] %>% round(3)
)

dt <- cbind(dt.1,dt.2) %>% setorder(-AUC)


##########
## Save ##
##########

if ("strand"%in%colnames(dt)) {
  dt[,strand:=NULL]
}

fwrite(dt, args$outfile, sep="\t")

##########
## TEST ##
##########

# R.utils::sourceDirectory("/Users/ricard/git/ArchR/R/", verbose=T, modifiedOnly=FALSE)

# ArchRProj = ArchRProject
# groupBy = "group"
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
