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

io$motif_bed             <- file.path(io$rawdata, "/processed/atac/signac/chromvar/motifbedchr19.tsv.gz")

io$cells_out             <- file.path(io$rawdata, "/processed/atac/signac/chromvar/cells.mtx")
io$motifs_out            <- file.path(io$rawdata, "/processed/atac/signac/chromvar/motifs.tsv.gz")
io$signac_out            <- file.path(io$rawdata, "/processed/atac/signac/chromvar/signacchr19.rds")

io$chromvar_out          <- file.path(io$rawdata, "/processed/atac/signac/chromvar/chromvarchr19.tsv.gz")


opts$cores               <- 8
opts$mem                 <- 5 # GB of global memory 
opts$block_size          <- 1e4 # number of regions to keep in memory during processing
opts$min_cells_per_locus <- 100


dir.create(dirname(io$cells_out), recursive = TRUE)


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

cells <- merged_meta[pass_atacQC == TRUE & pass_rnaQC == TRUE, cell]





# load fragments
frags <- CreateFragmentObject(
  io$merged_fragment_file,
  cells = cells
)
frags


motifs <- fread(io$motif_bed) %>% 
  .[, row := paste0(chr, "-", start, "-", end)]

# iterate over TF
tfs <- motifs[, unique(TF)]



tmpdir <- file.path(dirname(io$signac_out), "temp")
dir.create(tmpdir, recursive = TRUE)


tmpfiles <- map(tfs, ~{
  
  print(paste("creating matrix for", .x))
  
  gr <- makeGRangesFromDataFrame(motifs[TF == .x], keep.extra.columns = TRUE)
  
  # remove loci on nonstandard chromosomes and in genomic blacklist regions
  gr <- keepStandardChromosomes(gr, pruning.mode = "coarse")
  gr <- subsetByOverlaps(gr, ranges = blacklist_mm10, invert = TRUE)
  
  
  
  
  # quantify counts
  matrix <- FeatureMatrix(
    fragments = frags,
    features = gr,
    cells = cells,
    process_n = 1e4
  )
  
  #matrix[1:10, 1:10]
  
  # now filter to remove loci with zero reads
  print("filter out positions without low reads...")
  rows <- apply(matrix, 1, function(x) sum(x>0) > opts$min_cells_per_locus)
  table(rows)
  
  if (sum(rows) < 100) {
    print("not enough loci with good coverage. Skipping...")
    return(NULL)
  }
  
  print(paste("keeping", sum(rows), "loci"))
  matrix <- matrix[rows, ]
  
  # now construct motif x 'peak' matrix
  motifmat <- motifs[TF == .x] %>% 
    .[, row := paste0(chr, "-", start, "-", end)] %>% 
    .[row %in% rownames(matrix)] %>% 
    .[, .(row, TF)]
  
  # save temp files
  print("saving matrix..")
  tmpfile_mat <- paste0(tmpdir, "/", .x, ".mtx")
  writeMM(matrix, tmpfile_mat)
  
  print("saving motifs..")
  tmpfile_mot <- paste0(tmpdir, "/", .x, ".tsv.gz")
  fwrite(motifmat, tmpfile_mot, sep = "\t", na = "NA", quote = FALSE)
  
  tmpfile_cols <- paste0(tmpdir, "/", .x, "_cols.rds")
  saveRDS(colnames(matrix), tmpfile_cols)
  
  tmpfile_rows <- paste0(tmpdir, "/", .x, "_rows.rds")
  saveRDS(rownames(matrix), tmpfile_rows)
  
  c(tmpfile_mat, tmpfile_mot, tmpfile_cols, tmpfile_rows)
})

tmpfiles <- purrr::compact(tmpfiles)

cols <- map(tmpfiles, 3) %>% 
  map(readRDS)
rows <- map(tmpfiles, 4) %>% 
  map(readRDS)

print("checking colnames match with cells..")

map(cols, ~{
  all(.x == cells)
})

cellmat <- map(tmpfiles, 1) %>% 
  map(readMM) %>% 
  map2(cols, ~{
    colnames(.x) <- .y
    .x
  }) %>% 
  map2(rows, ~{
    rownames(.x) <- .y
    .x
  }) %>% 
  purrr::reduce(rbind)

cellmat[1:5,1:5]

motifmat <- map(tmpfiles, 2) %>% 
  map(fread) %>% 
  rbindlist() %>% 
  .[, v := 1] %>% 
  dcast(row ~ TF, value.var = "v", fill = 0) 






dim(cellmat)
dim(motifmat)


writeMM(cellmat, io$cells_out)
fwrite(motifmat, io$motifs_out, sep = "\t", na = "NA", quote = FALSE)
fwrite(as.data.table(rownames(cellmat)), gsub(".mtx", "_rownames.txt", io$cells_out), quote = FALSE)
fwrite(as.data.table(colnames(cellmat)), gsub(".mtx", "_colnames.txt", io$cells_out), quote = FALSE)



meta_df <- copy(merged_meta) %>%
  setDF() %>%
  tibble::column_to_rownames("cell") %>%
  .[colnames(cellmat),]

meta_df[1:10, 1:10]

positions <- motifs[row %in% rownames(cellmat)] %>% 
  setDF(rownames = .$row) %>% 
  .[rownames(cellmat), ] %>% 
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)


motifmat <- as.matrix(motifmat, rownames = "row") %>% 
  .[rownames(cellmat), ]
  




chrom_assay <- CreateChromatinAssay(
  counts = cellmat,
  motifs = CreateMotifObject(data = motifmat),
  ranges = positions
  #fragments = frags,
  # annotation = annotation
)

saveRDS(chrom_assay, io$signac_out)


chrom_assay <- RunChromVAR(object = chrom_assay, genome = BSgenome.Mmusculus.UCSC.mm10)

chromvar <- as.data.table(chrom_assay@data, keep.rownames = "tf")
fwrite(chromvar, io$chromvar_out, sep = "\t", na = "NA", quote = FALSE)


print("deleting temp files..")
file.remove(unlist(tmpfiles))
