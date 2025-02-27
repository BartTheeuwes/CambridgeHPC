library(Seurat)
library(Signac)
library(purrr)
library(data.table)
library(Matrix)
library(GenomicRanges)
library(EnsDb.Mmusculus.v79)
library(future)

# makes Seurat object from accesibility data using merged fragment file
# data is quantified over bins so as to avoid using CellRanger peaks which don't
# match between samples



source(here::here("settings.R"))
source(here::here("atac/Signac/signac_settings.R"))


opts$cores               <- 8
opts$mem                 <- 5 # GB of global memory 
opts$binsize             <- 1e4 # 10kb bins 
opts$block_size          <- 1e4 # number of regions to keep in memory during processing




plan("multiprocess", workers = opts$cores)
options(future.globals.maxSize = opts$mem * 1024 ^ 3)
plan()

# load metadata
sample_metadata <- fread(io$metadata)
sample_metadata

# merge with CellRanger stats
cell_metrics <- fread(io$merged_metrics_file)
cols <- c("cell", colnames(cell_metrics)[!colnames(cell_metrics) %in% colnames(sample_metadata)])

merged_meta <- merge(sample_metadata, cell_metrics[, .SD, .SDcol = cols], by = "cell", all.y = TRUE)



# load genome annotation


annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotation) <- "UCSC"

# annotation <- fread(io$gene_metadata) %>%
#   setnames("symbol", "gene") %>%
#   makeGRangesFromDataFrame(keep.extra.columns = TRUE)


genome(annotation) <- "mm10"

# load fragments
frags <- CreateFragmentObject(
  io$merged_fragment_file,
  cells = merged_meta$cell
)
frags

# make bin x cell matrix
genome_lengths <- seqlengths(annotation)

mat <- GenomeBinMatrix(
  fragments = frags,
  genome = genome_lengths,
  binsize = opts$binsize,
  process_n = opts$block_size
)

mat[1:10, 1:10]

writeMM(mat, io$matrix_out)


meta_df <- copy(merged_meta) %>%
  setDF() %>%
  tibble::column_to_rownames("cell") %>%
  .[colnames(mat),]

meta_df[1:10, 1:10]



chrom_assay <- CreateChromatinAssay(
  counts = mat,
  fragments = frags,
  annotation = annotation
)

seurat <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "bins",
  meta.data = meta_df
)


seurat@meta.data$sample %>% unique()


saveRDS(seurat, io$signac_rds)

