here::i_am("atac/archR/peak_calling/peak_calling_archR.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

suppressPackageStartupMessages(library(ArchR))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--archr_directory',    type="character",    help='ArchR directory')
p$add_argument('--metadata',    type="character",    help='metadata file')
#p$add_argument('--pathToMacs2',     type="character",    help='Path to MACS2 software')
p$add_argument('--group_by',     type="character",    help='Metadata column to group by')
p$add_argument('--pvalue_cutoff',     type="double",   help='MACS2 p-value cutoff')
p$add_argument('--extend_summits',     type="integer",   help='Number of bp to extend peak summits')
p$add_argument('--threads',     type="integer",    default=1,    help='Number of threads')

args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# args$archr_directory <- file.path(io$basedir,"processed/atac/archR")
# args$metadata <- "/bi/group/reik/ricard/data/eomes_10x_multiome/results/atac/archR/qc/sample_metadata_after_qc.txt.gz"
# args$pathToMacs2 <- "/bi/group/reik/ricard/software/miniconda3/envs/main/bin/macs2"
# args$group_by <- "celltype_genotype"
# args$pvalue_cutoff <- 1e-3
# args$extend_summits <- 300
# args$threads <- 1
## END TEST ##

args$min_cells <- 50

dir.create(file.path(args$archr_directory,"PeakCalls"), showWarnings=F)

########################
## Load cell metadata ##
########################

sample_metadata <- fread(args$metadata) %>%
  .[pass_atacQC==TRUE & doublet_call==FALSE | pass_atacQC==TRUE & is.na(doublet_call)]

# temporary
# sample_metadata[,day:=stringr::str_replace_all(sample,opts$sample2day)] %>%
#   .[,genotype:=stringr::str_replace_all(sample,opts$sample2genotype)] %>%
#   .[,day_genotype:=paste(day,genotype,sep="_")] %>%
#   .[,celltype_genotype:=paste(celltype.mapped_mnn,genotype,sep="_")]

# Filter groups by minimum number of cells
stopifnot(args$group_by%in%colnames(sample_metadata))
sample_metadata <- sample_metadata[!is.na(sample_metadata[[args$group_by]])]
sample_metadata <- sample_metadata[!grepl("NA",sample_metadata[[args$group_by]])]
sample_metadata <- sample_metadata[,N:=.N,by=c(args$group_by)] %>% .[N>=args$min_cells] %>% .[,N:=NULL]

table(sample_metadata[[args$group_by]])

########################
## Load ArchR project ##
########################

# source(here::here("atac/archR/load_archR_project.R"))

setwd(args$archr_directory)

addArchRGenome("mm10")
addArchRThreads(threads = args$threads)

ArchRProject <- loadArchRProject(args$archr_directory)
# Subset
ArchRProject.filt <- ArchRProject[sample_metadata$cell]

# Unnecesary, but just to make sure projectMetadata is updated
tmp <- file.path(args$archr_directory,"projectMetadata.rds")
if (file.exists(tmp)) ArchRProject.filt@projectMetadata <- readRDS(tmp)

###########################
## Update ArchR metadata ##
###########################

sample_metadata.to.archr <- sample_metadata %>% 
  .[cell%in%rownames(ArchRProject.filt@cellColData)] %>% setkey(cell) %>% .[rownames(ArchRProject.filt@cellColData)] %>%
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
# table(getCellColData(ArchRProject,"Sample")[[1]])
table(getCellColData(ArchRProject.filt,args$group_by)[[1]])

############################
## update projectMetadata ##
############################
# Contains path to groupcoverage files of ArchRProject.filt
ArchRProject.filt@projectMetadata = readRDS(file.path(args$archr_directory,"projectMetadata.rds"))

##################
## Peak calling ##
##################
args$pathToMacs2 = findMacs2()

ArchRProject.filt <- addReproduciblePeakSet(
  ArchRProj = ArchRProject.filt, 
  groupBy = args$group_by, 
  peakMethod = "Macs2",
  excludeChr = c("chrM", "chrY"),
  pathToMacs2 = args$pathToMacs2,
  cutOff = args$pvalue_cutoff,
  extendSummits = args$extend_summits,
  plot = FALSE,
  force = TRUE
)

################
## Save peaks ##
################

# Save PeakSet
# NOTE THAT THIS HAS TO GO TO ARCHR DIRECTORY
saveRDS(ArchRProject.filt@peakSet, file.path(args$archr_directory,"/PeakCalls/PeakSet.rds"))

# fetch peaks in data.table format
dt <- getPeakSet(ArchRProject.filt) %>% as.data.table() %>% setnames("seqnames","chr")

# Save peak metadata
fwrite(dt, file.path(args$archr_directory,"PeakCalls/peak_metadata.tsv.gz"), sep="\t")

# save peaks in bed format
fwrite(dt[,c("chr","start","end")], file.path(args$archr_directory,"PeakCalls/peaks_archR_macs2.bed.gz"), sep="\t", col.names = F)

#####################
## Add peak matrix ##
#####################

ArchRProject@peakSet <- ArchRProject.filt@peakSet
ArchRProject <- addPeakMatrix(ArchRProject, binarize = FALSE)
saveArchRProject(ArchRProject)