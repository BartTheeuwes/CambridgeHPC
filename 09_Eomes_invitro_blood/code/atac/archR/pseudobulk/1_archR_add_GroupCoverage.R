# https://www.ArchRProject.com/bookdown/how-does-archr-make-pseudo-bulk-replicates.html
here::i_am("atac/archR/pseudobulk/1_archR_add_GroupCoverage.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

suppressPackageStartupMessages(library(ArchR))
suppressPackageStartupMessages(library(BSgenome.Mmusculus.UCSC.mm10))

#rhdf5::h5disableFileLocking()


######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--metadata',    type="character",    help='metadata file')
p$add_argument('--archr_directory',    type="character",    help='ArchR directory')
p$add_argument('--group_by',     type="character",    help='Metadata column to group by')
p$add_argument('--min_cells',     type="integer",    default=50,   help='Minimum number of cells')
p$add_argument('--max_cells',     type="integer",    default=1000,   help='Maximum number of cells')
p$add_argument('--threads',     type="integer",    default=1,    help='Number of threads')

args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# args$archr_directory <- file.path(io$basedir,"processed/atac/archR")
# args$metadata <- file.path(io$basedir,"results/atac/archR/qc/sample_metadata_after_qc.txt.gz")
# args$group_by <- "celltype_genotype"
# args$min_cells <- 50
# args$max_cells <- 9999
# args$threads <- 1
## END TEST ##

########################
## Load cell metadata ##
########################

cell_metadata.dt <- fread(args$metadata) %>%
  .[pass_atacQC==TRUE & doublet_call==FALSE | pass_atacQC==TRUE & is.na(doublet_call)]

# Filter groups by minimum number of cells
stopifnot(args$group_by%in%colnames(cell_metadata.dt))
cell_metadata.dt <- cell_metadata.dt[!is.na(cell_metadata.dt[[args$group_by]])]
cell_metadata.dt <- cell_metadata.dt[!grepl("NA",cell_metadata.dt[[args$group_by]])]
cell_metadata.dt <- cell_metadata.dt[,N:=.N,by=c(args$group_by)] %>% .[N>=args$min_cells] %>% .[,N:=NULL]

table(cell_metadata.dt[[args$group_by]])

########################
## Load ArchR project ##
########################

# source(here::here("atac/archR/load_archR_project.R"))

setwd(args$archr_directory)

addArchRGenome("mm10")
addArchRThreads(threads = args$threads)

ArchRProject <- loadArchRProject(args$archr_directory)

# Subset
ArchRProject.filt <- ArchRProject[cell_metadata.dt$cell]

###########################
## Update ArchR metadata ##
###########################

cell_metadata.dt.to.archr <- cell_metadata.dt %>% 
  .[cell%in%rownames(ArchRProject.filt@cellColData)] %>% setkey(cell) %>% .[rownames(ArchRProject.filt@cellColData)] %>%
  as.data.frame() %>% tibble::column_to_rownames("cell")

stopifnot(all(rownames(cell_metadata.dt.to.archr) == rownames(getCellColData(ArchRProject.filt))))
ArchRProject.filt <- addCellColData(
  ArchRProject.filt,
  data = cell_metadata.dt.to.archr[[args$group_by]],
  name = args$group_by,
  cells = rownames(cell_metadata.dt.to.archr),
  force = TRUE
)

# print cell numbers
table(getCellColData(ArchRProject.filt,"Sample")[[1]])
table(getCellColData(ArchRProject.filt,args$group_by)[[1]])

#########################
## Add Group Coverages ##
#########################
# Check if group Coverages already exist
# ArchRProject@projectMetadata$GroupCoverages

# This function will merge cells within each designated cell group for the generation of pseudo-bulk replicates 
# and then merge these replicates into a single insertion coverage file.
# Output: creates files in archR/GroupCoverages/celltype: [X]._.Rep[Y].insertions.coverage.h5
ArchRProject.filt <- addGroupCoverages(ArchRProject.filt, 
  groupBy = args$group_by,
  useLabels = TRUE,  # Need to use sample information, otherwise only 50 cells get used per condition. (bug in ArchR)
  minCells = args$min_cells,
  maxCells = args$max_cells,
  force = TRUE
)

##########
## Save ##
##########

saveRDS(ArchRProject.filt@projectMetadata, file.path(args$archr_directory,"projectMetadata.rds"))
#saveArchRProject(ArchRProject)