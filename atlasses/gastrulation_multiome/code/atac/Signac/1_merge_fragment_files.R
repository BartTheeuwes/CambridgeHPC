library(purrr)
library(data.table)


source(here::here("settings.R"))
source(here::here("atac/Signac/signac_settings.R"))

# script to merge multiple CellRanger fragment files
# barcodes are re-named to ensure they are unique between samples

# this is required since Signac has no way to handle multiple fragment files 
# with overlapping barcode names

# only cells which pass basic QC are kept to limit file size

# script requires 'htslib' to be installed/loaded in order to perform compression and indexing




opts$cores            <- 8
opts$min_fragments    <- 1e3



metrics <- map(seq_along(io$samples), ~{
  
  original_name <- io$samples[[.x]]
  new_name <- names(io$samples)[[.x]]
  metrics <- io$cell_metrics_files[[.x]]
    # filter to only include cells passing a lenient QC threshold
    # create cell-name to match Ricard's naming scheme
    
  fread(metrics) %>% 
     # .[is_cell == 1] %>% 
    .[atac_fragments >= opts$min_fragments] %>% 
    .[, c("sample", "filename") := .(new_name, original_name)] %>% 
    .[, cell := paste0(sample, "_", barcode)]
    
    
  }) %>% 
    rbindlist()

fwrite(metrics, io$merged_metrics_file, sep = "\t", na = "NA", quote = FALSE)
#fwrite(metrics, gsub(".tsv", "_all.tsv", io$merged_metrics_file), sep = "\t", na = "NA", quote = FALSE)


# .x=io$fragment_files[[1]]
# .y="multiome1"
fragments <- map(seq_along(io$samples), ~{
  
  original_name <- io$samples[[.x]]
  new_name <- names(io$samples)[[.x]]
  frags <- io$fragment_files[[.x]]
  
  
  
  cells <- metrics[sample == new_name, .(barcode, cell)]
  print(paste("number of cells to subset:", nrow(cells)))
  print(paste("reading fragment file:", original_name))
  
  dt <- fread(cmd = paste("zcat", frags))
  setnames(dt, c("chr", "start", "end", "barcode", "count"))
  
  setkey(dt, barcode)
  print(paste("number of cells in fragment file:", dt[, length(unique(barcode))]))
  setkey(cells, barcode)
  
  dt <- merge(dt, cells, by = "barcode")
  
  print(paste("number of cells in subsetted fragment file:", dt[, length(unique(cell))]))
  
  dt[, .(chr, start, end, cell, count)]
  
}) %>% 
  rbindlist() 

print(paste("number of cells in merged fragment file:", fragments[, length(unique(cell))]))

print("sorting by position...")
setorder(fragments, chr, start, end)

print("saving as tsv...")

tsv <- gsub(".gz", "", io$merged_fragment_file)

fwrite(fragments, tsv, col.names = FALSE, sep = "\t", quote = FALSE)

print("compressing...")
system(paste("bgzip -f -@", opts$cores, tsv))

print("indexing...")
system(paste0("tabix --preset=bed -f ", io$merged_fragment_file))
