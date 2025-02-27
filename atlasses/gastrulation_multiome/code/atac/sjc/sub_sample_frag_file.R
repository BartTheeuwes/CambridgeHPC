library(data.table)
library(purrr)

source(here::here("settings.R"))


fragments_file <- file.path(io$rawdata, "/processed/atac/signac/merged_fragments.tsv.gz")

opts$n_cells <- 100
#opts$rows_to_read <- opts$n_cells * 1e4

outfile <- gsub(".tsv.gz", paste0("_", opts$n_cells, "cells.tsv.gz"), fragments_file)


dt <- fread(fragments_file)#, nrows = opts$rows_to_read)
cells <- dt[, unique(V4)]
cells <- cells[1:opts$n_cells]
dt <- dt[V4 %in% cells]

print("saving")
print(outfile)
fwrite(dt, outfile, sep = "\t", na = "NA", quote = FALSE)
print("done")