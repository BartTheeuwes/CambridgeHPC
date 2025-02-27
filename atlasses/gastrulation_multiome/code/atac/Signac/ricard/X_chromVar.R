
#################
## Description ##
#################

# loads a metadata file and incorporates it into the Seurat object

####################
## Load libraries ##
####################

library(Seurat)
library(Signac)
library(future)
library(chromVAR)
library(motifmatchr)
library(JASPAR)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)


#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
} else {
  stop("Computer not recognised")
}

# I/O
io$metadata <- paste0(io$basedir, "/sample_metadata.txt.gz")
io$outdir <- paste0(io$basedir, "/results/atac/signac/chromvar")
io$signac_out    <- file.path(io$basedir, "/processed/atac/signac/signac_chromvar.rds")
# io$plots_out     <- "/bi/home/clarks/plots/10X_multiome/signac/footprints"

# Options
opts$cores       <- 1
opts$mem         <- 5 # GB
opts$assay       <- "peaks"

# Multiprocessing
# plan("multiprocess", workers = opts$cores)
# options(future.globals.maxSize = opts$mem * 1024 ^ 3)
# plan()

#################
## Load Signac ##
#################

signac <- readRDS(io$signac)

# switch assay
DefaultAssay(signac) <- opts$assay
signac

############################################
## Load Position-specific Weight Matrices ##
############################################

# extract position frequency matrices for the motifs
pfm <- getMatrixSet(JASPAR, opts = list(species = "Homo sapiens"))
length(pfm)

# convert to data.table
# pfm_names <- as.data.table(name(pfm), keep.rownames = TRUE) %>%
#   setnames(c("motif", "tf")) %>%
#   .[, tf := gsub("::|\\(|\\)|\\.", "_", tf) %>% gsub("_$", "", .)] %>%
#   .[, rows := motif] %>%
#   setDF() %>%
#   tibble::column_to_rownames("rows") 
  
mapping <- fread("ftp://ftp.ebi.ac.uk/pub/databases/mofa/10x_rna_atac_vignette/JASPAR_mapping.txt")
# mapping <- fread("/Users/ricard/data/JASPAR/JASPAR_mapping.txt")
stopifnot(names(pfm) %in% mapping$id)
foo <- mapping$name; names(foo) <- mapping$id
names(pfm) <- foo[names(pfm)] %>% paste0(.,"-motif")

###################################
# Add motifs to the Signac oject ##
###################################

signac <- AddMotifs(
  object = signac,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = pfm
)
signac

##################
## Run ChromVAR ##
##################

signac <- RunChromVAR(
  object = signac,
  genome = BSgenome.Mmusculus.UCSC.mm10
)
signac

saveRDS(signac, io$signac_out)
stop()

# ???
DefaultAssay(signac) <- "chromvar"
pfm_names <- pfm_names[rownames(signac),]
signac[["chromvar"]][["tf"]] <- pfm_names$tf


#################
## Save output ##
#################

print("saving...")
saveRDS(signac, io$signac_out)
fwrite(pfm_names, paste0(io$signac_out, "_motif_names.tsv"), sep = "\t", na = "NA")

sig_meta <- signac@meta.data %>% 
  cbind(signac[["umap"]]@cell.embeddings) %>% 
  cbind(as.data.frame(t(signac[["chromvar"]][]))) %>% 
  setDT(keep.rownames = "cell")

fwrite(sig_meta, io$metadata_out, sep = "\t", na="NA")

