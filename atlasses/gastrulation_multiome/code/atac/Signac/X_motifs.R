library(Seurat)
library(Signac)
library(purrr)
library(data.table)
library(future)
library(JASPAR)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)
library(patchwork)



source(here::here("settings.R"))



io$signac        <- file.path(io$rawdata, "/processed/atac/signac/signac_peaks.rds")
io$outdir        <- file.path(io$rawdata, "/processed/atac/signac/motifs/")

opts$cores       <- 8
opts$mem         <- 50 # GB
opts$celltype1   <- "Gut"
opts$celltype2   <- "Def._endoderm"


io$outdir <- file.path(io$outdir, paste(opts$celltype1, opts$celltype2, sep = "-")) %>% 
  gsub("\\.", "", .)
dir.create(io$outdir, recursive = TRUE)



plan("multiprocess", workers = opts$cores)
options(future.globals.maxSize = opts$mem * 1024 ^ 3)
plan()

signac <- readRDS(io$signac)
signac


# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR,
  opts = list(species = 9606, all_versions = FALSE)
)

# add motif information
signac <- AddMotifs(
  object = signac,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = pfm
)


da_peaks <- FindMarkers(
  object = signac,
  ident.1 = opts$celltype1,
  ident.2 = opts$celltype2,
  only.pos = TRUE,
  test.use = 'LR',
  latent.vars = 'nCount_peaks'
)
head(da_peaks)
# get top differentially accessible peaks
top.da.peak <- rownames(da_peaks[da_peaks$p_val < 0.005, ])


# test enrichment
enriched.motifs <- FindMotifs(
  object = signac,
  features = top.da.peak
)
head(enriched.motifs)

MotifPlot(
  object = signac,
  motifs = head(rownames(enriched.motifs))
)


motifs_out <- file.path(io$outdir, "enriched_motifs.rds")

saveRDS(enriched.motifs, motifs_out)



