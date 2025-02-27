library(Seurat)
library(Signac)
library(purrr)
library(data.table)
library(furrr)




source(here::here("settings.R"))

io$signac           <- file.path(io$rawdata, "/processed/atac/signac/archr_signac.rds")
io$diffacc_out      <- file.path(io$rawdata, "/results/atac/signac/diffacc_each_vs_all.tsv.gz")
dir.create(dirname(io$diffacc_out), recursive = TRUE)







# load signac

signac <- readRDS(io$signac) 
signac
Idents(signac) <- signac$celltype.mapped

celltypes <- unique(signac$celltype.mapped)



options(future.globals.maxSize = object.size(signac) * 2)
cores <- availableCores()
cores
plan(multicore, workers = cores)

tmpdir <- paste0(dirname(io$diffacc_out), "/diffacctmpeach")
dir.create(tmpdir)

diffacc <- future_map(celltypes, ~{
  # differential testing
  
  celltype1 <- .x
  
  print(paste("testing", celltype1))
  
  types <- data.table(old = Idents(signac)) %>% 
    .[old == celltype1, new := celltype1] %>% 
    .[old != celltype1, new := "control"]
  
  Idents(signac) <- types$new
  
  
  
  da_peaks <- FindMarkers(
    object = signac,
    ident.1 = celltype1,
    ident.2 = "control",
    min.pct = 0.2,
    test.use = 'LR',
    latent.vars = "nCount_peaks"
  )
  
  tmpfile <- paste0(tmpdir, "/", celltype1, ".tsv.gz")
  
  da_peaks <- setDT(da_peaks, keep.rownames = "locus")
  da_peaks[, c("celltype1") := .(celltype1)]
  fwrite(da_peaks, tmpfile, sep = "\t", na = "NA", quote = FALSE)
  
  tmpfile
  
})  %>% 
  unlist()

rm(signac)


diffacc_dt <- map(diffacc, fread) %>% rbindlist()

fwrite(diffacc_dt, io$diffacc_out, sep = "\t", na = "NA", quote = FALSE)

map(diffacc, file.remove)


