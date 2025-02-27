library(Seurat)
library(Signac)
library(purrr)
library(data.table)
library(future)
library(GenomeInfoDb)
library(GenomicRanges)



source(here::here("settings.R"))


io$signac        <- file.path(io$rawdata, "/processed/atac/signac/signac.rds")
io$signac_out    <- file.path(io$rawdata, "/processed/atac/signac/signac_peaks_celltype.rds")
io$peaks_out     <- file.path(io$rawdata, "/processed/atac/signac/signac_peaks_per_celltype_predicted.rds")

opts$cores       <- 8
opts$mem         <- 5 # GB
opts$group_by    <- "celltype.predicted"#"seurat_clusters"#"celltype.predicted"#"seurat_clusters"#"celltype.mapped" # split data by this before peak calling (or NULL)

opts$run_macs2     <- TRUE # if FALSE, load the peaks from the peaks file
opts$add_to_signac <- FALSE 
opts$update_signac <- FALSE



plan("multiprocess", workers = opts$cores)
options(future.globals.maxSize = opts$mem * 1024 ^ 3)
plan()

signac <- readRDS(io$signac)
signac
signac@meta.data$sample %>% unique()

if (opts$run_macs2) {
  macs2 <- system("which macs2", intern = TRUE)
  macs2
  
  peaks <- CallPeaks(
    object = signac,
    group.by = opts$group_by,
    macs2.path = macs2
  )
  
  # remove peaks on nonstandard chromosomes and in genomic blacklist regions
  peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
  peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_mm10, invert = TRUE)
  
  saveRDS(peaks, io$peaks_out)
} else {
  peaks <- readRDS(io$peaks_out)
}



if (opts$add_to_signac){
  # quantify counts in each peak
  macs2_counts <- FeatureMatrix(
    fragments = Fragments(signac),
    features = peaks,
    cells = colnames(signac)
  )
  
  # create a new assay using the MACS2 peak set and add it to the Seurat object
  signac[["peaks"]] <- CreateChromatinAssay(
    counts = macs2_counts,
    fragments = Fragments(signac),
    annotation = Annotation(signac)
  )
  
  if (opts$update_signac){
    print("saving signac...")
    saveRDS(signac, io$signac)
  }
  
  
  print("saving peaks-only signac...")
  DefaultAssay(signac) <- "peaks"
  signac[["bins"]] <- NULL
  signac@meta.data$sample %>% unique()
  
  saveRDS(signac, io$signac_out)
}






