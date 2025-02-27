library(SnapATAC)
library(purrr)
library(data.table)
library(GenomicRanges)

source(here::here("settings.R"))


######## needs to be run with R/3.5.2 ########

# script adapted from https://github.com/r3fang/SnapATAC/blob/master/examples/10X_brain_5k/README.md

### input/output ###

io$snap_file             <- file.path(io$basedir, "/processed/atac/snapatac/snapatac.rds")
io$gene_anno             <- "/bi/scratch/Stephen_Clark/annotations/Mmusculus_genes_BioMart.87.txt"

opts$cores               <- 8

snap <- readRDS(io$snap_file)

# load gene annotation
genes <- fread(io$gene_anno) %>%
  setnames("symbol", "name") %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)
genes

# add to snap
snap <- createGmatFromMat(snap, genes = genes, num.cores = opts$cores)
snap
snap <- makeBinary(snap, mat = "gmat")


saveRDS(snap, io$snap_file)