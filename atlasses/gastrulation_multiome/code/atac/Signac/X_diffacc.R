library(Seurat)
library(Signac)
library(purrr)
library(data.table)
library(future)




source(here::here("settings.R"))



io$signac        <- file.path(io$rawdata, "/processed/atac/signac/signac_peaks.rds")
io$outdir        <- file.path(io$rawdata, "/processed/atac/signac/diffacc")



opts$grouping    <- "celltype.predicted"
opts$group1      <- c("Erythroid1", "Erythroid2", "Erythroid3")
opts$group2      <- "Forebrain_Midbrain_Hindbrain"


dir.create(io$outdir, recursive = TRUE)

gr1 <- paste(opts$group1, collapse = "_")
gr2 <- paste(opts$group2, collapse = "_")


meta <- fread(io$metadata) %>% 
  .[, .SD, .SDcol = c("cell", opts$grouping)] %>% 
  setnames(c("cell", "group")) %>% 
  .[group %in% opts$group1, group := gr1] %>% 
  .[group %in% opts$group2, group := gr2] %>% 
  setDF(rownames = .$cell)


signac <- readRDS(io$signac)
signac

meta <- meta[colnames(signac),]
Idents(signac) <- meta$group

da_peaks <- FindMarkers(
  object = signac,
  ident.1 = gr1,
  ident.2 = gr2,
  min.pct = 0.2,
  test.use = 'LR',
  latent.vars = 'nCount_peaks' 
)


setDT(da_peaks, keep.rownames = "locus")

outfile <- paste0(io$outdir, "/", gr1, "_vs_", gr2, "_difacc.tsv.gz")


fwrite(da_peaks, outfile, sep = "\t", quote = FALSE, na = "NA")


fwrite_tsv <- function(path, ...){
  if (!dir.exists(dirname(path))){
    create.dir(dirname(path), recursive = TRUE)
  }
  fwrite(path, sep = "\t", na = "NA", quote = FALSE, ...)
}



