here::i_am("atac/archR/processing/0_create_arrow_files.R")

source(here::here("settings.R"))

suppressPackageStartupMessages(library(ArchR))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--samples',           type="character",  nargs='+',      help='Samples')
p$add_argument('--fragments_files',           type="character",  nargs='+',      help='ATAC Fragments files')
p$add_argument('--genome',           type="character", default="mm10",      help='Genome')
p$add_argument('--min_fragments',     type="integer",    default=1000,   help='Minimum number of ATAC fragments')
p$add_argument('--max_fragments',     type="integer",    default=1e7,    help='Maximum number of ATAC fragments')
p$add_argument('--min_tss_score',   type="double",     default=2.5,    help='Minimum TSS score threshold')
p$add_argument('--threads',     type="integer",    default=1,    help='Number of threads')
p$add_argument('--outdir',          type="character",                               help='Output directory')

args <- p$parse_args(commandArgs(TRUE))


## START TEST ##
# args$samples <- opts$samples
# # args$samples <- c("rv_eo_deg_day3_5_control")
# args$fragments_files <- file.path(io$basedir,sprintf("original/%s/out/atac_fragments.tsv.gz",args$samples))
# # args$fragments_files <- file.path(io$basedir,"original/rv_eo_deg_day3_5_dtag/atac_fragments.tsv.gz")
# args$genome <- "mm10"
# args$min_fragments <- 1000
# args$max_fragments <- 1000000
# args$min_tss_score <- 2.5
# args$threads <- 4
# args$outdir <- file.path(io$basedir,"processed/atac/archR")
## END TEST ##

dir.create(args$outdir, showWarnings=F, recursive=T)

#####################
## Define settings ##
#####################

setwd(args$outdir)

# ArchR options
addArchRThreads(threads=args$threads) 
addArchRGenome(args$genome)

#rhdf5::h5disableFileLocking()

# list of chromosomes has to be defined a priori
# stopifnot(!is.null(opts$chr))

########################
## create Arrow Files ##
########################
# test 
#args$fragments_files = paste0(io$basedir, '/original/', args$sample, '/outs/atac_fragments.tsv.gz')
args$fragments_files

ArrowFiles <- createArrowFiles(
  inputFiles = args$fragments_files,
  sampleNames = args$samples,
  outputNames = args$samples,
  addTileMat = TRUE,
  addGeneScoreMat = FALSE,
  excludeChr = c("chrM", "chrY"),

  #subThreading = FALSE, # parallel processing doesn't work well (https://github.com/GreenleafLab/ArchR/issues/248)
  force = TRUE,

  # QC metrics
  minFrags = args$min_fragments,  # The minimum number of fragments per cell
  maxFrags = args$max_fragments,  # The maximum number of fragments per cell
  minTSS = args$min_tss_score   # The minimum TSS enrichment score per cell
)

##################
## Sanity check ##
##################

# for (i in ArrowFiles) {
#   gr <- getFragmentsFromArrow(i)
#   stopifnot(sort(unique(seqnames(gr)))==sort(opts$chr))
# }


##########
## TEST ##
##########

# setwd("/bi/group/reik/ricard/data/eomes_10x_multiome/processed/atac/archR")

# library(rhdf5)
# ArrowFile <- "1A_Eo_DEG_G9_day3.arrow"
# ArrowFile <- "backup/rv_eo_deg_day3_5_control.arrow"
# fid <- H5Fopen(ArrowFile)
# h5ls(ArrowFile)
# h5read(fid, "Metadata")
# h5read(fid, "Metadata/CellNames") # 4470
# h5read(fid, "Metadata/Sample")
# h5write("E8.25_PijuanSala", file=fid, name="Metadata/Sample")
# h5closeAll()
