library(Seurat)
library(Signac)
library(purrr)
library(data.table)
library(future)
library(chromVAR)

library(motifmatchr)
library(JASPAR)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)



source(here::here("settings.R"))



io$signac        <- file.path(io$rawdata, "/processed/atac/signac/signac_peaks.rds")
io$signac_out    <- file.path(io$rawdata, "/processed/atac/signac/signacMotifs_human.rds")
io$metadata_out  <- file.path(io$rawdata, "/processed/atac/signac/signac_chromVar.tsv.gz")
io$plots_out     <- "/bi/home/clarks/plots/10X_multiome/signac/footprints"


opts$cores       <- 8
opts$mem         <- 10 # GB

opts$species     <- "human"

if (toupper(opts$species) == "HUMAN") opts$species <- 9606
if (toupper(opts$species) == "MOUSE") opts$species <- 10090


plan("multiprocess", workers = opts$cores)
options(future.globals.maxSize = opts$mem * 1024 ^ 3)
plan()



signac <- readRDS(io$signac)
signac





# extract position frequency matrices for the motifs
pfm <- getMatrixSet(
  x = JASPAR,
  opts = list(species = opts$species, all_versions = FALSE)
)


length(pfm)

pfm_names <- as.data.table(name(pfm), keep.rownames = TRUE) %>%
  setnames(c("motif", "tf")) %>%
  .[, tf := gsub("::|\\(|\\)|\\.", "_", tf) %>% gsub("_$", "", .)] %>%
  .[, rows := motif] %>%
  setDF() %>%
  tibble::column_to_rownames("rows") 
  



signac <- AddMotifs(
  object = signac,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = pfm
)
signac

signac <- RunChromVAR(
  object = signac,
  genome = BSgenome.Mmusculus.UCSC.mm10
)
signac

DefaultAssay(signac) <- "chromvar"
pfm_names <- pfm_names[rownames(signac),]


signac[["chromvar"]][["tf"]] <- pfm_names$tf


# save output 
print("saving...")
saveRDS(signac, io$signac_out)
fwrite(pfm_names, paste0(io$signac_out, "_motif_names.tsv"), sep = "\t", na = "NA")

sig_meta <- signac@meta.data %>% 
  cbind(signac[["umap"]]@cell.embeddings) %>% 
  cbind(as.data.frame(t(signac[["chromvar"]][]))) %>% 
  setDT(keep.rownames = "cell")

fwrite(sig_meta, io$metadata_out, sep = "\t", na="NA")






