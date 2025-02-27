library(data.table)
library(purrr)


source(here::here("settings.R"))

io$samples            <-        c(E8.5_rep1 = "multiome1",
                                  E8.5_rep2 = "multiome2",
                                  E7.5_rep1 = "rep1_L001_multiome",
                                  E7.5_rep2 = "rep2_L002_multiome") 


io$fragment_files     <- file.path(io$rawdata, 
                                   io$samples, 
                                   "/outs/atac_fragments.tsv.gz")

io$cell_metrics_files <- file.path(io$rawdata, 
                                   io$samples, 
                                   "/outs/per_barcode_metrics.csv")



io$merged_fragment_file  <- file.path(io$rawdata, "/processed/atac/signac/merged_fragments.tsv.gz")
io$merged_metrics_file          <- file.path(io$rawdata, "/processed/atac/signac/cell_metrics.tsv")

io$matrix_out            <- file.path(io$rawdata, "/processed/atac/signac/bins.mtx")
io$signac_rds               <- file.path(io$rawdata, "/processed/atac/signac/signac.rds")





dir.create(dirname(io$merged_fragment_file), recursive = TRUE)
dir.create(dirname(io$merged_metrics_file), recursive = TRUE)





