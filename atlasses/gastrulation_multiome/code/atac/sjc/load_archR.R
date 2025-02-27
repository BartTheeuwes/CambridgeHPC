library(data.table)
library(purrr)
library(ArchR)

source(here::here("settings.R"))


io$archr_dir <- file.path(io$rawdata, "/processed/atac/archR")

dir(io$archr_dir, recursive = T)

setwd(io$archr_dir)

bulk <- readRDS("pseudobulk/pseudobulk_PeakMatrix_summarized_experiment.rds" ) %>% 
  assays() %>% 
  .$PeakMatrix %>% 
  as.data.table(keep.rownames = "locus") %>% 
  melt(id.vars = "locus", value.name = "acc", variable.name = "celltype")


addArchRThreads(threads = 2) 
addArchRGenome("mm10")

ArchRProject <- loadArchRProject(io$archr_dir)
head(ArchRProject@cellColData)

ArchRProject[1:100]

mat <- getMatrixFromProject(ArchRProject[1:100], useMatrix = "PeakMatrix")
mat[1:15,1:5]
