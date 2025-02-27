library(Seurat)
library(Signac)
library(purrr)
library(data.table)
library(Matrix)
library(GenomicRanges)
library(EnsDb.Mmusculus.v79)
library(future)

# makes Seurat object from accesibility data using motif positions



source(here::here("settings.R"))


io$merged_fragment_file  <- file.path(io$rawdata, "/processed/atac/signac/merged_fragments.tsv.gz")
io$metrics_file          <- file.path(io$rawdata, "/processed/atac/signac/cell_metrics.tsv")
io$motif_bed             <- file.path(io$rawdata, "/processed/atac/signac/chromvar/motifbed.tsv.gz")

io$matrix_out            <- file.path(io$rawdata, "/processed/atac/signac/chromvar/motifmat.mtx")
io$outfile               <- file.path(io$rawdata, "/processed/atac/signac/chromvar/signac.rds")

opts$cores               <- 8
opts$mem                 <- 5 # GB of global memory 
opts$block_size          <- 1e4 # number of regions to keep in memory during processing

dir.create(dirname(io$outfile), recursive = TRUE)


plan("multiprocess", workers = opts$cores)
options(future.globals.maxSize = opts$mem * 1024 ^ 3)
plan()

# load metadata
sample_metadata <- fread(io$metadata)
sample_metadata

# merge with CellRanger stats
cell_metrics <- fread(io$metrics_file)
cols <- c("cell", colnames(cell_metrics)[!colnames(cell_metrics) %in% colnames(sample_metadata)])

merged_meta <- merge(sample_metadata, cell_metrics[, .SD, .SDcol = cols], by = "cell", all = TRUE)

cell_metrics[, passQC := FALSE] %>% 
  .[atac_raw_reads > 1e3, passQC := TRUE]

ggplot(cell_metrics, aes(sample, fill = passQC, colour = passQC)) +
  geom_bar()

# load genome annotation


annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotation) <- "UCSC"
genome(annotation) <- "mm10"

# load fragments
frags <- CreateFragmentObject(
  io$merged_fragment_file,
  cells = merged_meta$cell
)
frags

motifs <- fread(io$motifs_bed)

# remove loci on nonstandard chromosomes and in genomic blacklist regions
motifs <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
motifs <- subsetByOverlaps(x = peaks, ranges = blacklist_mm10, invert = TRUE)



meta <- fread(io$metadata) 
cells <- meta[, pass]


# quantify counts
macs2_counts <- FeatureMatrix(
  fragments = Fragments(signac),
  features = motifs,
  cells = cells
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


saveRDS(seurat, io$outfile)

