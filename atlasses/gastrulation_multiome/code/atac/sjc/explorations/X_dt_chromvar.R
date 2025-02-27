library(Seurat)
library(Signac)
library(purrr)
library(data.table)
library(Matrix)
library(GenomicRanges)
library(EnsDb.Mmusculus.v79)
library(future)
library(BSgenome.Mmusculus.UCSC.mm10)

# makes Seurat object from accesibility data using motif positions



source(here::here("settings.R"))


io$merged_fragment_file  <- file.path(io$rawdata, "/processed/atac/signac/merged_fragments.tsv.gz")
#io$merged_fragment_file <- "/bi/scratch/Stephen_Clark/multiome/raw//processed/atac/signac/merged_fragments_100cells.tsv.gz"
io$metrics_file          <- file.path(io$rawdata, "/processed/atac/signac/cell_metrics.tsv")

io$motif_bed             <- file.path(io$rawdata, "/processed/atac/signac/chromvar/motifbedchr19.tsv.gz")

io$cells_out             <- file.path(io$rawdata, "/processed/atac/signac/chromvar/overlapchr19.tsv.gz")

io$motifs_out            <- file.path(io$rawdata, "/processed/atac/signac/chromvar/motifs.tsv.gz")
io$signac_out            <- file.path(io$rawdata, "/processed/atac/signac/chromvar/signacchr19.rds")

io$chromvar_out          <- file.path(io$rawdata, "/processed/atac/signac/chromvar/chromvarchr19.tsv.gz")


opts$cores               <- 8
opts$mem                 <- 5 # GB of global memory 
opts$block_size          <- 1e4 # number of regions to keep in memory during processing
opts$min_cells_per_locus <- 100
opts$binarise            <- TRUE


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
frags <- fread(io$merged_fragment_file, select = 1:4) %>% 
  setnames(c("chr", "start", "end", "cell")) %>% 
  .[cell %in% cells] 
frags


motifs <- fread(io$motif_bed)
  

motif_names <- unique(motifs[, .(ID, TF)])

# remove loci on nonstandard chromosomes and in genomic blacklist regions
motifs <- motifs[chr %in% paste0("chr", c(1:19, "X", "Y")), .(chr, start, end, TF)] %>% 
  setkey(chr, start , end)
blacklist <- as.data.table(blacklist_mm10) %>% 
  setnames("seqnames", "chr") %>% 
  setkey(chr, start, end)

motifs <- foverlaps(motifs, blacklist) %>% 
  .[is.na(start), .(chr, start = i.start, end = i.end, TF)] %>% 
  setkey(chr, start, end)

motifs


tmpdir <- file.path(dirname(io$cells_out), "tmp")
dir.create(tmpdir, recursive = TRUE)

# iterate over cells to save memory

overlap <- split(frags, by = "cell", keep.by = TRUE) %>% 
  map(~{
    foverlaps(.x, motifs, nomatch = 0L) %>% 
      .[, .N, .(chr, start, end, TF, cell)] %>% 
      # .[, ncells := .N, .(chr, start, end, TF)] %>% 
      # .[ncells >= opts$min_cells_per_locus] %>% 
      .[, .(id = paste0(chr, "_", start, "_", end, "_", TF),
            cell,
            N)]
  }) %>% rbindlist()



fwrite(overlap, io$cells_out, sep = "\t", na = "NA", quote = FALSE)



# convert to cell x locus matrix

if (opts$binarise) {
  overlap[N > 0, N := 1]
}


overlap <- dcast(overlap, id ~ cell, value.var = "N", fill = 0) %>% 
  as.matrix(rownames = "id") %>% 
  Matrix(sparse = TRUE)

motifmat <- data.table(id = rownames(overlap)) %>% 
  .[, tf := gsub("[^_]*_", "", id)] %>% 
  .[, val := 1] %>% 
  dcast(id ~ tf, value.var = "val", fill = 0) %>% 
  as.matrix(rownames = "id") %>% 
  .[rownames(overlap), ]


all(rownames(motifmat) == rownames(overlap))
nrow(motifmat)==nrow(overlap)

meta_df <- copy(merged_meta) %>%
  setDF() %>%
  tibble::column_to_rownames("cell") %>%
  .[colnames(overlap),]

meta_df[1:10, 1:10]

positions <- data.table(id = rownames(motifmat)) %>% 
  tidyr::separate("id", c("chr", "start", "end"), remove = FALSE, extra = "drop") %>% 
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)
 
all(rownames(motifmat) == positions$id)
  
  
  
# test <- CreateMotifObject(data = motifmat)
# all(rownames(test@data) == positions$id)
# 
# dim(test)
# dim(overlap)
# length(positions)
# 
# overlap[1:5,1:5]
# positions[1:5]
# test@data[1:5,1:2]



chrom_assay <- CreateChromatinAssay(
  counts = overlap,
  #motifs = test,
  ranges = positions
  #fragments = frags,
  # annotation = annotation
)
chrom_assay@motifs <- CreateMotifObject(data = motifmat)



saveRDS(chrom_assay, io$signac_out)


chrom_assay <- RunChromVAR(object = chrom_assay, genome = BSgenome.Mmusculus.UCSC.mm10)

chromvar <- as.data.table(chrom_assay@data, keep.rownames = "tf")
fwrite(chromvar, io$chromvar_out, sep = "\t", na = "NA", quote = FALSE)


print("deleting temp files..")
file.remove(unlist(tmpfiles))
