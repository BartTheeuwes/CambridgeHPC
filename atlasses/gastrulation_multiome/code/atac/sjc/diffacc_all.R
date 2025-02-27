library(Seurat)
library(Signac)
library(purrr)
library(data.table)
library(furrr)




source(here::here("settings.R"))

io$signac           <- file.path(io$rawdata, "/processed/atac/signac/archr_signac.rds")
io$diffacc_out      <- file.path(io$rawdata, "/results/atac/signac/diffacc.tsv.gz")
dir.create(dirname(io$diffacc_out), recursive = TRUE)







# load signac

signac <- readRDS(io$signac) 
signac
Idents(signac) <- signac$celltype.mapped

celltypes <- unique(signac$celltype.mapped)

combs <- expand.grid(celltypes, celltypes) %>% 
  as.data.table() %>% 
  setnames(c("celltype1", "celltype2")) %>% 
  .[celltype1 != celltype2] %>% 
  .[, map(.SD, as.character)]

options(future.globals.maxSize = object.size(signac) * 2)
cores <- availableCores()
cores
plan(multicore, workers = cores)

tmpdir <- paste0(dirname(io$diffacc_out), "/diffacctmp")
dir.create(tmpdir)

diffacc <- future_map(1:nrow(combs), ~{
  # differential testing
  
  celltype1 <- combs[.x, celltype1]
  celltype2 <- combs[.x, celltype2]
  
  print(paste("testing:", celltype1, "vs", celltype2))
  
  da_peaks <- FindMarkers(
    object = signac,
    ident.1 = celltype1,
    ident.2 = celltype2,
    min.pct = 0.2,
    test.use = 'LR',
    latent.vars = "nCount_peaks"
  )
  
  tmpfile <- paste0(tmpdir, "/", celltype1, celltype2, ".tsv.gz")
  
  da_peaks <- setDT(da_peaks, keep.rownames = "locus")
  da_peaks[, c("celltype1", "celltype2") := .(celltype1, celltype2)]
  fwrite(da_peaks, tmpfile, sep = "\t", na = "NA", quote = FALSE)
  tmpfile
})  %>% unlist()

rm(signac)


diffacc_dt <- map(diffacc, fread) %>% rbindlist()

fwrite(diffacc_dt, io$diffacc_out, sep = "\t", na = "NA", quote = FALSE)

map(diffacc, file.remove)


