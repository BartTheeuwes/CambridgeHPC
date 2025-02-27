# script containing i/o for snapatac scripts
source(here::here("settings.R"))
snapio <- list()

#snapio$frags_file <- file.path(io$basedir,  "public_datasets/Pijuan-Sala_/data/fragments.tsv.gz")
# for testing

# the name of the sample which is also the folder the raw data will be found in:
snapio$samples       <-        c("multiome1",
                                 "multiome2",
                                 "rep1_L001_multiome",
                                 "rep2_L002_multiome") 
#snapio$samples <- "multiome1"

snapio$frags_files   <- file.path(io$rawdata, 
                                  snapio$samples, 
                                  "/outs/atac_fragments.tsv.gz")

snapio$metrics_files <- file.path(io$rawdata, 
                                  snapio$samples, 
                                  "/outs/per_barcode_metrics.csv")
                                             
snapio$snap_files    <- file.path(io$rawdata, 
                                  snapio$samples, 
                                  "/parsed/snap.snap")



snapio$rds_file      <- file.path(io$basedir, "processed/atac/snapatac/snap.rds")
snapio$peaks         <- file.path(io$basedir, "processed/atac/snapatac/peaks.tsv.gz")
snapio$difacc        <- file.path(io$basedir, "processed/atac/snapatac/difacc.tsv.gz")
snapio$motifs        <- file.path(io$basedir, "processed/atac/snapatac/motifs.tsv.gz")

snapio$metadata      <- file.path(io$basedir, "processed/atac/snapatac/snap_metadata.tsv.gz")

snapio$plots_out     <- "/bi/home/clarks/plots/10X_multiome/snapatac/"


dir.create(snapio$plots_out, recursive = TRUE)

dir_create <- function(path){
  if (!dir.exists(path)) {
    print("creating directory")
    print(path)
    dir.create(path, recursive = TRUE)
  }
}
walk(unlist(snapio), ~dir_create(dirname(.x)))
  

