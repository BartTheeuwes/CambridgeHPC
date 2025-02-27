library(Seurat)
library(Signac)
library(purrr)
library(data.table)
library(future)
library(cowplot)



source(here::here("settings.R"))



io$signac        <- file.path(io$basedir, "/processed/atac/signac/signac_peaks.rds")
        

io$plots_out     <- "/bi/home/clarks/plots/10X_multiome/Signac"

opts$cores       <- 8
opts$mem         <- 50 # GB




dir.create(io$plots_out, recursive = TRUE)


plan("multiprocess", workers = opts$cores)
options(future.globals.maxSize = opts$mem * 1024 ^ 3)
plan()

seurat <- readRDS(io$seurat)
seurat

signac <- readRDS(io$signac)
signac

# remove bins for memoery efficiency
signac[["bins"]] <- NULL
signac

# match up cells 
cells <- colnames(seurat)[colnames(seurat) %in% colnames(signac)]

signac <- signac[, cells]
seurat <- seurat[, cells]

signac[["RNA"]] <- seurat[["RNA"]]

