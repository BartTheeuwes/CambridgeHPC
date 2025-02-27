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

# Load in data
out_folder = "/rds/project/bg200/rds-bg200-hphi-gottgens/users/bt392/mouse/Mixl1_KO/data/"
big_meta = read.table(paste0(out_folder,"big_meta.tab"), header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "$")
pca_before = as.data.frame(readRDS(paste0(out_folder, 'uncorrected_PCs.rds')))
pca_before$cell <- rownames(pca_before)
pca_after = as.data.frame(read.csv(paste0(out_folder, 'corrected_pc_complete.csv'), header=T, row.names=1,sep=","))
pca_after$cell <- rownames(pca_after)

# add metadata to PCs
pca_before <- merge(pca_before, big_meta, by='cell')
rownames(pca_before) <- pca_before$cell
pca_after <- merge(pca_after, big_meta, by='cell')
rownames(pca_after) <- pca_after$cell

# calculate tsne's on whole object and chimeara seperately
### chimaera only
# before correction
pca_before$origin <- sapply(strsplit(pca_before$cell,"_"), `[`, 1)
chim_before <- pca_before %>% filter(origin=='chim') %>% select(2:31) # use first 30 PCs
chim_before_tsne = Rtsne(chim_before, pca = FALSE)$Y
chim_before_tsne <- as.data.frame(chim_before_tsne)
colnames(chim_before_tsne) <- c('tsne1','tsne2')
chim_before_tsne$cell <- rownames(chim_before)
chim_before_tsne <- merge(chim_before_tsne, big_meta, by='cell')
# write
write.csv(chim_before_tsne, paste0(out_folder, 'tsne_before_correction_chim.csv')) 
print('finished first tsne')


# after correction
pca_after$origin <- sapply(strsplit(pca_after$cell,"_"), `[`, 1)
chim_after <- pca_after %>% filter(origin=='chim') %>% select(2:31) # use first 30 PCs
chim_after_tsne = Rtsne(chim_after, pca = FALSE)$Y
chim_after_tsne<- as.data.frame(chim_after_tsne)
colnames(chim_after_tsne) <- c('tsne1','tsne2')
chim_after_tsne$cell <- rownames(chim_after)
chim_after_tsne <- merge(chim_after_tsne, big_meta, by='cell')
# write
write.csv(chim_after_tsne, paste0(out_folder, 'tsne_after_correction_chim.csv')) 
print('finished second tsne')

### whole object
# before correction
pca_before$origin <- sapply(strsplit(pca_before$cell,"_"), `[`, 1)
pca_before <- pca_before %>% select(2:31) # use first 30 PCs
pca_before_tsne = Rtsne(pca_before, pca = FALSE)$Y
pca_before_tsne <- as.data.frame(pca_before_tsne)
colnames(pca_before_tsne) <- c('tsne1','tsne2')
pca_before_tsne$cell <- rownames(pca_before)
pca_before_tsne <- merge(pca_before_tsne, big_meta, by='cell')
# write
write.csv(pca_before_tsne, paste0(out_folder, 'tsne_before_correction.csv')) 
print('finished third tsne')

# after correction
pca_after$origin <- sapply(strsplit(pca_after$cell,"_"), `[`, 1)
pca_after <- pca_after %>% select(2:31) # use first 30 PCs
pca_after_tsne = Rtsne(pca_after, pca = FALSE)$Y
pca_after_tsne <- as.data.frame(pca_after_tsne)
colnames(pca_after_tsne) <- c('tsne1','tsne2')
pca_after_tsne$cell <- rownames(pca_after)
pca_after_tsne <- merge(pca_after_tsne, big_meta, by='cell')
# write
write.csv(pca_after_tsne, paste0(out_folder, 'tsne_after_correction.csv')) 
print('finished fourth tsne')