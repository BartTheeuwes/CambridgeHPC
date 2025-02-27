library(SnapATAC)
library(purrr)
library(data.table)



source(here::here("settings.R"))

io$snap_file            <- file.path(io$basedir, "/processed/atac/snapatac/snapatac.rds")

io$metadata             <- file.path(io$basedir, "/sample_metadata.csv")
io$barcode_metrics_file <- file.path(io$basedir, "/processed/atac/snapatac/barcode_metrics.tsv.gz")

io$outfile              <- file.path(io$basedir, "/sample_metadata_SnapATAC.csv")

opts$cols               <- c("pass_snapQC", 
                             "snap_cluster", 
                             "promoter_ratio", 
                             "atac_fragments", 
                             "atac_MT")



snap_meta <- copy(snap@metaData) %>%
  setDT() %>%
  .[, .(barcode, sample, cluster)] %>%
  setnames("cluster", "snap_cluster") %>%
  .[, pass_snapQC := TRUE] %>%
  .[, cell := paste0(sample, "_", gsub("-.*", "", barcode))] %>%
  .[, c("barcode", "sample") := NULL]


metrics <- fread(io$barcode_metrics_file) %>%
  .[, atac_MT := atac_mitochondrial_reads/atac_fragments] %>%
  .[, cell := paste0(sample, "_", gsub("-.*", "", barcode))] %>%
  .[, c("barcode", "sample") := NULL]

merged <- merge(snap_meta, metrics, by = "cell", all = TRUE) %>%
  .[is.na(pass_snapQC), pass_snapQC := FALSE] %>%
  .[, .SD, .SDcol = c("cell", opts$cols)]
  

sample_metadata <- fread(io$metadata) 

# first remove old data if it is there

cols <- colnames(sample_metadata)[colnames(sample_metadata) %in% opts$cols]
if (length(cols) > 0) sample_metadata[, cols := NULL]

# now add new metadata and save
new <- merge(sample_metadata, merged, all.x = TRUE, by = "cell")
new
fwrite(new, io$outfile, sep = "\t", na = "NA")


