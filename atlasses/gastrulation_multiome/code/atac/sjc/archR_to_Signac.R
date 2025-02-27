library(data.table)
library(purrr)
library(ArchR)
library(Signac)
library(Seurat)
library(BSgenome.Mmusculus.UCSC.mm10)
library(EnsDb.Mmusculus.v79)

source(here::here("settings.R"))


io$archr_dir         <- file.path(io$rawdata, "/processed/atac/archR")
io$fragments_file    <- file.path(io$rawdata, "/processed/atac/signac/merged_fragments.tsv.gz")

io$outfile           <- file.path(io$rawdata, "/processed/atac/signac/archr_signac.rds")


# load archR project
ArchRProject <- loadArchRProject(io$archr_dir)




# extract counts matrix as summarised experiment (se)
se <- getMatrixFromProject(ArchRProject, useMatrix = "PeakMatrix")
# now extract matrix and metadata from se
mat <- assay(se)
gr <- se@rowRanges
meta <- se@colData


# re-format cell names to match fragments file
colnames(mat)  <- gsub("#", "_", colnames(mat))
rownames(meta) <- gsub("#", "_", rownames(meta))


# load fragments
frags <- CreateFragmentObject(
  io$fragments_file,
  cells = rownames(meta)
)
frags

# load genome annotation
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotation) <- "UCSC"
genome(annotation) <- "mm10"

chrom_assay <- CreateChromatinAssay(
  counts = mat,
  ranges = gr,
  genome = "mm10",
  fragments = frags,
  annotation = annotation
)


signac <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = as.data.frame(meta)
)

# normalise
signac <- RunTFIDF(signac)

# run svd
signac <- FindTopFeatures(signac, min.cutoff = 'q0')
signac <- RunSVD(signac)
# run umap
signac <- RunUMAP(object = signac, reduction = 'lsi', dims = 2:30)
signac <- FindNeighbors(object = signac, reduction = 'lsi', dims = 2:30)
signac <- FindClusters(object = signac, verbose = FALSE, algorithm = 3)


# save
saveRDS(signac, io$outfile)
