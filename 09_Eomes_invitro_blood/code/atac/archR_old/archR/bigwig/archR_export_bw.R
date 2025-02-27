here::i_am("atac/archR/bigwig/archR_export_bw.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

suppressPackageStartupMessages(library(ArchR))
#rhdf5::h5disableFileLocking()

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--metadata',    type="character",    help='metadata file')
p$add_argument('--group_by',     type="character",    help='Metadata column to group by')
p$add_argument('--norm_method',     type="character", default="ReadsInTSS",    help='Normalisation method')
p$add_argument('--tile_size',     type="integer", default=100,    help='Tile size')
p$add_argument('--min_cells',     type="integer", default=100,    help='Minimum number of cells per celltype')
p$add_argument('--threads',     type="integer",    default=1,    help='Number of threads')

args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# args$metadata <- file.path(io$basedir,"results/atac/archR/qc/sample_metadata_after_qc.txt.gz")
# args$group_by <- "celltype_genotype"
# args$norm_method <- c("ReadsInTSS")
# args$tile_size <- 100
# args$min_cells <- 50
# args$threads <- 1
## END TEST ##

########################
## Load cell metadata ##
########################

sample_metadata <- fread(args$metadata) %>%
  .[pass_atacQC==TRUE & sample%in%opts$samples]

# temporary
# sample_metadata[,day:=stringr::str_replace_all(sample,opts$sample2day)] %>%
#   .[,genotype:=stringr::str_replace_all(sample,opts$sample2genotype)] %>%
#   .[,day_genotype:=paste(day,genotype,sep="_")] %>%
#   .[,celltype_genotype:=paste(celltype.mapped_mnn,genotype,sep="_")]

# Filter celltypes by minimum number of cells
stopifnot(args$group_by%in%colnames(sample_metadata))
sample_metadata <- sample_metadata[!is.na(sample_metadata[[args$group_by]])]
sample_metadata <- sample_metadata[!grepl("NA",sample_metadata[[args$group_by]])]
sample_metadata <- sample_metadata[,N:=.N,by=c(args$group_by)] %>% .[N>=args$min_cells] %>% .[,N:=NULL]

table(sample_metadata[[args$group_by]])

########################
## Load ArchR project ##
########################

source(here::here("atac/archR/load_archR_project.R"))

addArchRThreads(threads = args$threads)

# Subset
ArchRProject.filt <- ArchRProject[sample_metadata$cell]

# sapply(getArrowFiles(ArchRProject),h5ls)

###########################
## Update ArchR metadata ##
###########################

sample_metadata.to.archr <- sample_metadata %>% 
  .[cell%in%rownames(ArchRProject.filt)] %>% setkey(cell) %>% .[rownames(ArchRProject.filt)] %>%
  as.data.frame() %>% tibble::column_to_rownames("cell")

stopifnot(all(rownames(sample_metadata.to.archr) == rownames(getCellColData(ArchRProject.filt))))
ArchRProject.filt <- addCellColData(
  ArchRProject.filt,
  data = sample_metadata.to.archr[[args$group_by]],
  name = args$group_by,
  cells = rownames(sample_metadata.to.archr),
  force = TRUE
)

# print cell numbers
table(getCellColData(ArchRProject.filt,args$group_by)[[1]])


###################
## Export bigwig ##
###################

# This function will group, summarize and export a bigwig for each group in an ArchRProject.
getGroupBW(
  ArchRProj = ArchRProject.filt,
  groupBy = args$group_by,
  # groupBy = "Sample",
  normMethod = args$norm_method,
  tileSize = args$tile_size,
  maxCells = 1000, # default
  ceiling = 4
)

# Create a completion token
file.create(file.path(io$archR.directory,sprintf("/GroupBigWigs/%s/completed.txt",args$group_by)))
