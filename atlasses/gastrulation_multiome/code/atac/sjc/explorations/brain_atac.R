library(Seurat)
library(Signac)
library(purrr)
library(data.table)
library(future)
library(GenomeInfoDb)
library(GenomicRanges)



source(here::here("settings.R"))

io$signac        <- file.path(io$rawdata, "/processed/atac/signac/signac_peaks.rds")

opts$celltypes <- c(
  "Surface_ectoderm",
  "Rostral_neurectoderm", 
  "Neural_crest", 
  "Forebrain_Midbrain_Hindbrain" , 
  "Spinal_cord", 
  "Epiblast",
  "Caudal_epiblast"
)

meta <- fread(io$metadata)
meta[, unique(celltype.predicted)]

meta <- meta[pass_rnaQC == TRUE & celltype.predicted %in% opts$celltypes]
cells <- meta[, cell]



signac <- readRDS(io$signac) %>% 
  .[, cells]
signac

da_peaks <- FindMarkers(
  object = signac,
  ident.1 = "Forebrain_Midbrain_Hindbrain",
  ident.2 = "Surface_ectoderm",
  min.pct = 0.2,
  test.use = 'LR',
  latent.vars = "nCount_peaks"
)


da_peaks <- setDT(da_peaks, keep.rownames = "locus")
da_peaks

brain <- da_peaks[p_val_adj < 0.05 & avg_logFC > 0, locus]
skin  <- da_peaks[p_val_adj < 0.05 & avg_logFC < 0, locus]

signac <- signac[c(brain, skin), ]
signac
