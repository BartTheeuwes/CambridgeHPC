library(Seurat)
library(Signac)
library(purrr)
library(data.table)
library(future)
library(GenomeInfoDb)
library(GenomicRanges)



source(here::here("settings.R"))


io$signac        <- file.path(io$rawdata, "/processed/atac/signac/signac.rds")
io$bed_files     <- dir("/bi/scratch/Stephen_Clark/annotations/public_datasets/Metzis_2018_NMPs/",
                        f = TRUE,
                        pattern = ".bed")
io$signac_out    <- file.path(io$rawdata, "/processed/atac/signac/signac_over_nmp_anno.rds")


opts$group_by    <- "celltype.predicted"#"seurat_clusters"#"celltype.predicted"#"seurat_clusters"#"celltype.mapped" # split data by this before peak calling (or NULL)





plan("multiprocess", workers = availableCores())
options(future.globals.maxSize = 1 * 1024 ^ 3)
plan()

signac <- readRDS(io$signac)
signac
signac@meta.data$sample %>% unique()

anno <- map(io$bed_files, fread) %>% 
  map2(basename(io$bed_files), ~.x[, c("anno", "id") := .(.y, paste0(.y, "_", .I))]) %>% 
  rbindlist() %>% 
  setnames(paste0("V", 1:3), c("chr", "start", "end")) %>% 
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)




# quantify counts in each peak
counts <- FeatureMatrix(
    fragments = Fragments(signac),
    features = anno,
    cells = colnames(signac)
  )
  
# create a new assay using the MACS2 peak set and add it to the Seurat object
signac[["anno"]] <- CreateChromatinAssay(
    counts = counts,
    fragments = Fragments(signac),
    annotation = Annotation(signac)
  )
  

    

  
  
print("saving peaks-only signac...")
DefaultAssay(signac) <- "anno"

if (!is_null(signac[["bins"]])) signac[["bins"]] <- NULL
if (!is_null(signac[["peaks"]])) signac[["peaks"]] <- NULL





signac@meta.data$sample %>% unique()
  
saveRDS(signac, io$signac_out)







