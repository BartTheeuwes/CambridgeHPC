library(Seurat)
library(Signac)
library(purrr)
library(data.table)
library(future)
library(GenomeInfoDb)
library(GenomicRanges)



source(here::here("settings.R"))


io$signac        <- file.path(io$rawdata, "/processed/atac/signac/signac.rds")
io$anno_dir      <- "/bi/scratch/Stephen_Clark/annotations/gastrulation"
io$anno_files    <- file.path(io$anno_dir, c("H3K27ac_distal_E7.5_Ect_intersect12.bed",
                                             "H3K27ac_distal_E7.5_End_intersect12.bed",
                                             "H3K27ac_distal_E7.5_Mes_intersect12.bed"))
  
io$anno_out      <- file.path(io$rawdata, "/processed/atac/signac/k3k27ac/h3k27ac_anno.rds")
io$signac_out    <- file.path(io$rawdata, "/processed/atac/signac/k3k27ac/signac_h3k27ac.rds")

dir.create(dirname(io$anno_out), recursive = TRUE)

opts$cores       <- 8
opts$mem         <- 5 # GB

plan("multiprocess", workers = opts$cores)
options(future.globals.maxSize = opts$mem * 1024 ^ 3)
plan()

signac <- readRDS(io$signac)
signac
signac@meta.data$sample %>% unique()

h3k27ac <- map(io$anno_files, fread) %>% 
  rbindlist() %>% 
  setnames(c("chr", "start", "end","strand", "id", "anno")) %>% 
  .[, chr := paste0("chr", chr)] %>% 
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

h3k27ac

saveRDS(h3k27ac, io$anno_out)

# quantify counts in each locus
counts <- FeatureMatrix(
  fragments = Fragments(signac),
  features = h3k27ac,
  cells = colnames(signac)
)

# create a new assay using the annottion and add it to the Seurat object
signac[["h3k27ac"]] <- CreateChromatinAssay(
  counts = counts,
  fragments = Fragments(signac),
  annotation = Annotation(signac)
)
DefaultAssay(signac) <- "h3k27ac"
signac[["bins"]] <- NULL
signac[["peaks"]] <- NULL
signac

print("selecting features and normalising..")
signac <- FindTopFeatures(signac, min.cutoff = 5)
signac <- RunTFIDF(signac)
signac <- RunSVD(signac)

print("saving signac...")

saveRDS(signac, io$signac_out)






