#################
## Description ##
#################

# makes Seurat object from accesibility data using merged fragment file
# data is quantified over bins so as to avoid using CellRanger peaks which don't
# match between samples

####################
## Load libraries ##
####################

library(Seurat)
library(Signac)
library(Matrix)
library(GenomicRanges)
library(EnsDb.Mmusculus.v79)
library(future)

#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
} else {
  stop("Computer not recognised")
}

# I/O
io$metadata              <- paste0(io$basedir, "/results/rna/mapping/sample_metadata_after_mapping.txt.gz")
io$merged_fragment_file  <- paste0(io$basedir, "/processed/atac/signac/merged_fragments.tsv.gz")
io$merged_metrics_file   <- paste0(io$basedir, "/processed/atac/signac/cell_metrics.tsv.gz")
io$matrix_out            <- paste0(io$basedir, "/processed/atac/signac/bins.mtx")
io$outfile               <- paste0(io$basedir, "/processed/atac/signac/signac.rds")

# Options
opts$cores               <- 1
opts$mem                 <- 10 # GB
opts$binsize             <- 5e3 # 5kb bins 
opts$block_size          <- 1000 # number of regions to keep in memory during processing

# Multiprocessing
plan("multiprocess", workers = opts$cores)
options(future.globals.maxSize = opts$mem * 1024 ^ 3)
plan()

########################
## Load cell metadata ##
########################

# load cell metadata, including RNA-based info
io$metadata <- paste0(io$basedir,"/results/rna/mapping/sample_metadata_after_mapping.txt.gz")
sample_metadata <- fread(io$metadata)

# Load CellRanger stats
cell_metrics <- fread(io$merged_metrics_file)
cols <- c("cell", colnames(cell_metrics)[!colnames(cell_metrics) %in% colnames(sample_metadata)])

stopifnot(sort(sample_metadata$cell)==sort(cell_metrics$cell))

# Merge sample metadaa with cellRanger stats
# sample_metadata <- merge(sample_metadata, cell_metrics[, .SD, .SDcol = cols], by = "cell", all = TRUE)
sample_metadata <- sample_metadata %>% merge(cell_metrics[,..cols], by = "cell")

############################
## load genome annotation ##
############################

annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotation) <- "UCSC"
genome(annotation) <- "mm10"

####################
## load fragments ##
####################

frags <- CreateFragmentObject(io$merged_fragment_file, cells = sample_metadata$cell)
frags

###########################
## load motif annotation ##
###########################

# library(JASPAR)
# library(TFBSTools)
# 
# pfm <- getMatrixSet(JASPAR, opts = list(species = "Homo sapiens"))
# 
# mapping <- fread("/Users/ricard/data/JASPAR/JASPAR_mapping.txt")
# stopifnot(names(pfm) %in% mapping$id)
# foo <- mapping$name; names(foo) <- mapping$id
# names(pfm) <- foo[names(pfm)] %>% paste0(.,"-motif")
# 
# motif.matrix <- CreateMotifMatrix(
#   features = peaks.granges,
#   pwm = pfm,
#   genome = 'mm10',
#   use.counts = FALSE
# ) %>% as.matrix

############################################
## Create accessibility matrix using bins ##
############################################

print("Creating matrix...")

mat <- GenomeBinMatrix(
  fragments = frags,
  genome = seqlengths(annotation),
  binsize = opts$binsize,
  process_n = opts$block_size
)

writeMM(mat, io$matrix_out)

##########################
## Create Signac object ##
##########################

# prepare cell metadata
# stopifnot(metadata$cell%in%colnames(chrom_assay))
# metadata.to.seurat <- metadata %>%
#   .[cell%in%colnames(chrom_assay)] %>%
#   setkey(cell) %>% .[colnames(chrom_assay)] %>%
#   as.data.frame %>% tibble::column_to_rownames("cell")
# stopifnot(rownames(metadata.to.seurat)==colnames(chrom_assay))
meta_df <- sample_metadata %>% setDF %>%
  tibble::column_to_rownames("cell") %>%
  .[colnames(mat),]

print("Creating Signac object...")

chrom_assay <- CreateChromatinAssay(
  counts = mat,
  fragments = frags,
  annotation = annotation
  # motifs = CreateMotifObject(motif.matrix, pfm)
)

seurat <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "bins",
  meta.data = meta_df
)


##########
## Save ##
##########

print("Saving...")
saveRDS(seurat, io$outfile)

