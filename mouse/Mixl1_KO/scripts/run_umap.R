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
    library(M3C)
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

print('start running UMAPs')
# calculate umap's on whole object and chimeara seperately
### chimaera only
# before correction
pca_before$origin <- sapply(strsplit(pca_before$cell,"_"), `[`, 1)
chim_before <- pca_before %>% filter(origin=='chim') %>% select(2:51)
chim_before_umap = umap(t(chim_before))$data
chim_before_umap <- as.data.frame(chim_before_umap)
colnames(chim_before_umap) <- c('umap1','umap2')
chim_before_umap$cell <- rownames(chim_before)
chim_before_umap <- merge(chim_before_umap, big_meta, by='cell')
# write
write.csv(chim_before_umap, paste0(out_folder, 'umap_before_correction_chim.csv')) 
print('finished first umap')


# after correction
pca_after$origin <- sapply(strsplit(pca_after$cell,"_"), `[`, 1)
chim_after <- pca_after %>% filter(origin=='chim') %>% select(2:51) 
chim_after_umap = umap(t(chim_after))$data
chim_after_umap<- as.data.frame(chim_after_umap)
colnames(chim_after_umap) <- c('umap1','umap2')
chim_after_umap$cell <- rownames(chim_after)
chim_after_umap <- merge(chim_after_umap, big_meta, by='cell')
# write
write.csv(chim_after_umap, paste0(out_folder, 'umap_after_correction_chim.csv')) 
print('finished second umap')

### whole object
# before correction
pca_before$origin <- sapply(strsplit(pca_before$cell,"_"), `[`, 1)
pca_before <- pca_before %>% select(2:51)
pca_before_umap = umap(t(pca_before))$data
pca_before_umap <- as.data.frame(pca_before_umap)
colnames(pca_before_umap) <- c('umap1','umap2')
pca_before_umap$cell <- rownames(pca_before)
pca_before_umap <- merge(pca_before_umap, big_meta, by='cell')
# write
write.csv(pca_before_umap, paste0(out_folder, 'umap_before_correction.csv')) 
print('finished third umap')

# after correction
pca_after$origin <- sapply(strsplit(pca_after$cell,"_"), `[`, 1)
pca_after <- pca_after %>% select(2:51)
pca_after_umap = umap(t(pca_after))$data
pca_after_umap <- as.data.frame(pca_after_umap)
colnames(pca_after_umap) <- c('umap1','umap2')
pca_after_umap$cell <- rownames(pca_after)
pca_after_umap <- merge(pca_after_umap, big_meta, by='cell')
# write
write.csv(pca_after_umap, paste0(out_folder, 'umap_after_correction.csv')) 
print('finished fourth umap')