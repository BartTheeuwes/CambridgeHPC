# Batch correction for Mixl1 KO/WT chimaera + atlas
suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(Matrix)
    library(scran)
    library(Rtsne)
    library(BiocParallel)
    require(irlba)
    library(batchelor)
    library(stringr)
})

ncores = 4
mcparam = SnowParam(workers = ncores)
register(mcparam)

# out folder
out_folder = "/rds/project/bg200/rds-bg200-hphi-gottgens/users/bt392/mouse/Mixl1_KO/data/"

print('loading data')
# load data
big_sce <- readRDS(paste0(out_folder, "big_sce_norm.rds"))
big_meta <- read.table(paste0(out_folder,"big_meta.tab"), header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "$")
hvgs <- read.csv(paste0(out_folder,"big_hvgs_3000.csv"))$x

print('getting order')
# get order: oldest to youngest; most cells to least cells
order_df = big_meta[!duplicated(big_meta$sample), c("stage", "sample")]
order_df$ncells = sapply(order_df$sample, function(x) sum(big_meta$sample == x))
order_df$stage = factor(order_df$stage, 
                        levels = rev(c("E8.5", 
                                   "E8.25", 
                                   "E8.0", 
                                   "E7.75", 
                                   "E7.5", 
                                   "E7.25", 
                                   "mixed_gastrulation", 
                                   "E7.0", 
                                   "E6.75", 
                                   "E6.5")))
order_df = order_df[order(order_df$stage, order_df$ncells, decreasing = TRUE),]
order_df$stage = as.character(order_df$stage)

print('performing batch correction')
# perform batch correction
source("/rds/project/bg200/rds-bg200-hphi-gottgens/users/bt392/mouse/Mixl1_KO/core_own.R")
all_correct = doBatchCorrect(counts = logcounts(big_sce)[rownames(big_sce) %in% hvgs,], 
                             timepoints = big_meta$stage, 
                             samples = big_meta$sample, 
                             timepoint_order = order_df$stage, 
                             sample_order = order_df$sample, 
                             npc = 50,
                             BPPARAM = mcparam)


save(all_correct, file = paste0(out_folder,"all_correction.rds"))