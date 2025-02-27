suppressPackageStartupMessages({
    library(ArchR) 
    library(data.table)
    library(purrr)
    library(parallel)
    library(dplyr)
    library(ggpubr)
    library(BSgenome)
    library(biomaRt)
    library(argparse)
    library(gridExtra)
    library(viridis)
})

# I/O
io = list()
io$basedir='/rds/project/rds-SDzz0CATGms/users/bt392/04_Rabbit_ATAC_final'
io$archR.directory = file.path(io$basedir, 'ArchR/Project')
io$output.directory <- file.path(io$basedir,"ArchR")
setwd(io$output.directory)

addArchRThreads(threads = 1) 
addArchRVerbose(verbose = FALSE)

options(repr.plot.width=15, repr.plot.height=5)


# options
opts = list()
opts$samples = c('rabbit_BGRGP1', 
                 'rabbit_BGRGP2',
                 'rabbit_BGRGP3', 
                 'rabbit_BGRGP4', 
                 'rabbit_BGRGP5', 
                 'rabbit_BGRGP6', 
                 'rabbit_BGRGP7', 
                 'rabbit_BGRGP8') 

opts$stages = c('GD7',
                'GD8', 
                'GD9')

opts$stage.colors = c("GD7" = "#D56958",
                  "GD8" = "#6EC280",
                  "GD9" = "#6494D8",
                  "GD9_ExE" = "#BEB44D")