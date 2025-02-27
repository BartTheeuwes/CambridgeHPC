# TO-DO: USE OUTPUT OF SAVE ATAC MATRICES
here::i_am("atac/archR/feature_stats/archR_calculate_feature_stats.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

suppressPackageStartupMessages(library(ArchR))
suppressPackageStartupMessages(library(sparseMatrixStats))


######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--metadata',    type="character",    help='metadata file')
p$add_argument('--outfile',     type="character",    help='Output file')
p$add_argument('--group_by', type="character", help='Cell type label')
p$add_argument('--matrix',             type="character",     help='Matrix name')
p$add_argument('--pseudobulk_dir',             type="character",     help='Directory for the pseudobulk data')
p$add_argument('--ignore_small_celltypes',  action="store_true",  help='')
p$add_argument('--binarise',  action="store_true",  help='binarise single-cell ATAC counts?')
args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# args$metadata <- file.path(io$basedir,"results_new/atac/archR/celltype_assignment/sample_metadata_after_celltype_assignment.txt.gz")
# args$group_by <- "celltype.predicted"
# args$pseudobulk_dir <- file.path(io$basedir,"results_new/atac/archR/pseudobulk/celltype.mapped_mnn")
# args$ignore_small_celltypes <- TRUE
# args$matrix <- "PeakMatrix"
# args$binarise <- FALSE
# args$outfile <- file.path(io$basedir, sprintf("results_new/atac/archR/feature_stats/%s_feature_stats.txt.gz",args$matrix))
## END TEST ##

########################
## Load cell metadata ##
########################

if (grepl("genotype",args$group_by)) {
  sample_metadata <- fread(args$metadata) %>%
    .[,celltype_genotype:=sprintf("%s-%s",celltype.mapped,genotype)] %>%
    .[pass_atacQC==TRUE & doublet_call==FALSE]
} else {
  sample_metadata <- fread(args$metadata) %>%
    .[pass_atacQC==TRUE & doublet_call==FALSE & genotype=="WT"]
}

stopifnot(args$group_by%in%colnames(sample_metadata))
sample_metadata <- sample_metadata[!is.na(sample_metadata[[args$group_by]])]

# subset celltypes with sufficient number of cells
# if (args$ignore_small_celltypes) {
#   opts$celltypes <- names(which(table(sample_metadata[[args$group_by]])>30))
#   sample_metadata <- sample_metadata %>% .[eval(as.name(args$group_by))%in%opts$celltypes]
# }

print(table(sample_metadata[[args$group_by]]))

########################
## Load ArchR project ##
########################

source(here::here("atac/archR/load_archR_project.R"))

stopifnot(args$matrix%in%getAvailableMatrices(ArchRProject))

# Subset
ArchRProject.filt <- ArchRProject[sample_metadata$cell]

###################################
## Fetch single-cell ATAC Matrix ##
###################################

print(sprintf("Fetching single-cell ATAC %s matrix...",args$matrix))

atac.se <- getMatrixFromProject(ArchRProject.filt, useMatrix=args$matrix, binarize = args$binarise)
dim(atac.se)

# Define feature names
if (grepl("peak",tolower(args$matrix),ignore.case=T)) {
  row.ranges.dt <- rowRanges(atac.se) %>% as.data.table %>% 
    setnames("seqnames","chr") %>%
    .[,c("chr","start","end")] %>%
    .[,idx:=sprintf("%s:%s-%s",chr,start,end)]
  rownames(atac.se) <- row.ranges.dt$idx
} else if (grepl("gene",tolower(args$matrix),ignore.case=T)) {
  rownames(atac.se) <- rowData(atac.se)$name
} else {
  stop("Matrix not recognised")
}

########################################
## Fetch pseudobulk ATAC Peak Matrix ##
########################################

print(sprintf("Fetching pseudobulk ATAC %s matrix...",args$matrix))

atac_pseudobulk.se <- readRDS(file.path(args$pseudobulk_dir,sprintf("pseudobulk_%s_summarized_experiment.rds",args$matrix)))#[,opts$celltypes]

# Define feature names (already done in the pseudobulking script)
# if (grepl("peak",tolower(args$matrix),ignore.case=T)) {
#   rownames(atac_pseudobulk.se) <- rowData(atac_pseudobulk.se) %>% as.data.table %>% .[,idx:=sprintf("%s:%s-%s",seqnames,start,end)] %>% .$id
# } else if (grepl("gene",tolower(args$matrix),ignore.case=T)) {
#   rownames(atac_pseudobulk.se) <- rowData(atac_pseudobulk.se)$name
# }

###################
## Sanity checks ##
###################

stopifnot(rownames(atac.se)==rownames(atac_pseudobulk.se))

##########################
## Calculate peak stats ##
##########################

print("Calculating pseudobulk peak stats...")

featureStats.pseudobulk.dt <- data.table(
  feature = rownames(atac_pseudobulk.se),
  var = apply(assay(atac_pseudobulk.se),1,var) %>% round(5),
  mean = apply(assay(atac_pseudobulk.se),1,mean) %>% round(5)
  # min = apply(assay(atac_pseudobulk.se),1,min) %>% round(5),
  # max = apply(assay(atac_pseudobulk.se),1,max) %>% round(5)
)

print("Calculating singlecell peak stats...")

featureStats.singlecell.dt <- data.table(
  feature = rownames(atac.se),
  var = sparseMatrixStats::rowVars(assay(atac.se)) %>% round(5),
  mean = Matrix::rowMeans(assay(atac.se)) %>% round(5)
  # min = apply(assay(atac.se),1,min),
  # max = apply(assay(atac.se),1,max)
)

featureStats.dt <- merge(featureStats.singlecell.dt, featureStats.pseudobulk.dt, by="feature", suffixes=c("_singlecell","_pseudobulk"))

##########
## Save ##
##########

fwrite(featureStats.dt, args$outfile, sep="\t")
