library(Matrix)
library(Seurat)
library(Signac)
library(EnsDb.Mmusculus.v79)

#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/public_datasets/Pijuan-Sala_2020/settings.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/public_datasets/Pijuan-Sala_2020/settings.R")
}
io$outfile <- paste0(io$basedir,"/data/processed/seurat.rds")

opts$test <- FALSE

###################
## Load metadata ##
###################

metadata <- fread(io$metadata)

###############
## Load data ##
###############

# Error in as(object = list_of_data[[j]], Class = "dgCMatrix"): no method or default for coercing "ngTMatrix" to "dgCMatrix"
# inputdata <- Read10X(paste0(io$basedir,"/data"), gene.column = 1)

# Load matrix, features and barcodes
# 305187 23838 108569464  ((305187*23838*8) / 1e9= 58GB in standard matrix format)
m <- Matrix::readMM(io$matrix)
barcodes <- fread(io$barcodes, header=F)[[1]]
features <- fread(io$features, header=F)[[1]]

if (opts$test) {
  n <- 1000
  m <- m[1:n,1:n]
  features <- features[1:n]
  barcodes <- barcodes[1:n]
}
rownames(m) <- features
colnames(m) <- barcodes

# Convert from ngTMatrix (logical) to dgCMatrix (numeric)
m <- m*1

# Define GRanges object using the features
granges <- StringToGRanges(features, sep = c("_", "_"))
granges <- granges[as.vector(seqnames(granges) %in% standardChromosomes(granges)),]

# Define Granges object with gene annotations from ENSEMBL
ensembl.annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(ensembl.annotations) <- 'UCSC'  # NCBI for 1,2,X,   UCSC for chr1,chr2,chrX
genome(ensembl.annotations) <- "mm10"


###################################
## Create Seurat Chromatin Assay ##
###################################

chrom_assay <- CreateChromatinAssay(
  counts = m,
  # sep = c("_", "_"),
  ranges = granges,
  genome = 'mm10',
  fragments = io$fragments,
  min.cells = 0,
  min.features = 0,
  annotation = ensembl.annotations
)


# Prepare metadata
stopifnot(metadata$cell%in%colnames(chrom_assay))
metadata.to.seurat <- metadata %>%
  .[cell%in%colnames(chrom_assay)] %>%
  setkey(cell) %>% .[colnames(chrom_assay)] %>%
  as.data.frame %>% tibble::column_to_rownames("cell")
stopifnot(rownames(metadata.to.seurat)==colnames(chrom_assay))

# Create Seurat object
seurat <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata.to.seurat
)

##########
## Test ##
##########

seurat[["peaks"]]

##########
## Save ##
##########

saveRDS(seurat, io$outfile)

# create a gene activity matrix from the peak matrix and GTF, using chromosomes 1:22, X, and Y.
# Peaks that fall within gene bodies, or 2kb upstream of a gene, are considered
# activity.matrix <- CreateGeneActivityMatrix(peak.matrix = peaks, annotation.file = "../data/Homo_sapiens.GRCh37.82.gtf", 
#                                             seq.levels = c(1:22, "X", "Y"), upstream = 2000, verbose = TRUE)
