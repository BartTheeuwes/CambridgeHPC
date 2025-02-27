library(SnapATAC)
library(purrr)
library(data.table)
library(GenomicRanges)


source(here::here("settings.R"))
source(here::here("atac/SnapATAC/snapatac_settings.R"))


opts$cores    <- 8

snap <- readRDS(snapio$rds_file)

#snap@file <- gsub("/bi/scratch/Stephen_Clark/gastrulation_multiome_10x", io$basedir, snap@file)
snap

peaks <- fread(snapio$peaks) %>%
  makeGRangesFromDataFrame(seqnames.field = "V1", 
                           start.field = "V2", 
                           end.field = "V3")

(snaptools <- system("which snaptools", intern = TRUE))

runSnapAddPmat(snap, 
               peak = peaks,
               path.to.snaptools = snaptools,
               tmp.folder = tempdir(),
               num.cores = opts$cores)

snap <- addPmatToSnap(snap, do.par = TRUE, num.cores = opts$cores)
snap <- makeBinary(snap, mat = "pmat")
snap


saveRDS(snap, snapio$rds_file)
