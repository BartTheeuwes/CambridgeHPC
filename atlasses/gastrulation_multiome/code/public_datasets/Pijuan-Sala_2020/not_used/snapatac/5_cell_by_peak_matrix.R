library(SnapATAC)
library(purrr)
library(data.table)
library(GenomicRanges)

# requires snatptools, bedtools


source(here::here("settings.R"))
source(here::here("public_datasets/Pijuan-Sala_2020/snapatac/snapatac_settings.R"))

opts$cores    <- 8

print("loading rds file...")

snap <- readRDS(snapio$rds_file)

#snap@file <- gsub("/bi/scratch/Stephen_Clark/gastrulation_multiome_10x", io$basedir, snap@file)
snap

print("loading bed file...")

peaks <- fread(snapio$peaks) %>%
  makeGRangesFromDataFrame(seqnames.field = "V1", 
                           start.field = "V2", 
                           end.field = "V3")

peaks

(snaptools <- system("which snaptools", intern = TRUE))

print("adding peaks...")

runSnapAddPmat(snap,
               peak = peaks,
               path.to.snaptools = snaptools,
               tmp.folder = tempdir(),
               num.cores = opts$cores)

snap <- addPmatToSnap(snap, do.par = TRUE, num.cores = opts$cores)

print("binarising....")

snap <- makeBinary(snap, mat = "pmat")
snap

print("saving...")
saveRDS(snap, snapio$rds_file)
