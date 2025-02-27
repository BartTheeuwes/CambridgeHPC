# preparing extended atlas

suppressPackageStartupMessages({
    library(data.table)
    library(Matrix)
    library(dplyr)
})

main = '/rds/project/rds-SDzz0CATGms/users/bt392/atlasses/extended/'

print('reading meta')
meta = readRDS(paste0(main, 'integrated_meta_celltype_clus.rds'))
meta = meta[, c('cell', 'sample', 'stage', 'celltype.descendant', 'celltype.clustering')]
colnames(meta) = c('cell', 'sample', 'stage', 'celltype_original', 'celltype')

# +
#Downsample to 10k cells per sample (or maximum number of cells in sample)
set.seed(42)
keep = lapply(unique(meta$stage), function(x){
  if(x == "mixed_gastrulation"){
    return(c())
  } else if(sum(meta$stage == x) < 10000) {
    return(which(meta$stage == x))
  } else {
    hits = which(meta$stage == x)
    return(sample(hits, 10000))
  }
})
keep = do.call(c, keep)

meta = meta[keep,]
# -

print('writing meta')
write.table(meta, paste0(main, 'meta.tab'), sep='\t',row.names = FALSE)

print('reading counts')
counts = fread(paste0(main, 'integrated_counts_datatable.tab')) %>% .[rowSums(.[,-1]) > 0, ] %>% select(1, meta$cell)
genes = counts[,1]
colnames(genes) = 'V1'

print('writing genes')
write.table(genes, paste0(main, 'genes.tsv'), sep='\t',row.names = FALSE)

print('writing counts')
counts = as(as.matrix(counts[,-1]), "dgTMatrix")
writeMM(counts, paste0(main, 'raw_counts.mtx'))

# +
# Calculate size factors
# load packages
suppressPackageStartupMessages({
    library(Matrix)
    library(scran)
    library(scater)
    library(igraph)
    library(BiocParallel)
})
ncores = 16
mcparam = SnowParam(workers = ncores)
register(mcparam)

atlas_in = '/rds/project/rds-SDzz0CATGms/users/bt392/atlasses/extended/'
counts = readMM(paste0(atlas_in, 'raw_counts.mtx')) # readMM for .mtx input 


sce = SingleCellExperiment(assays = list("counts" = counts))

lib.sizes = Matrix::colSums(counts(sce))
sce = sce[calculateAverage(sce)>0.1,]

clusts = as.numeric(quickCluster(sce, method = "igraph", min.size = 100, BPPARAM = mcparam))

#now run the normalisation
#number of cells in each cluster should be at least twice that of the largest 'sizes'
min.clust = min(table(clusts))/2
new_sizes = c(floor(min.clust/3), floor(min.clust/2), floor(min.clust))
sce = computeSumFactors(sce, clusters = clusts, sizes = new_sizes, max.cluster.size = 3000)

write.table(sizeFactors(sce), quote = F, col.names = F, row.names = F, file = paste0(main, "sizefactors.tab"))
# -

ggplot(data = data.frame(X = lib.sizes, Y = sizeFactors(sce)),
              mapping = aes(x = X, y = Y)) +
  geom_point() +
  scale_x_log10(breaks = c(5000, 10000, 50000, 100000), labels = c("5,000", "10,000", "50,000", "100,000") ) +
  scale_y_log10(breaks = c(0.2, 1, 5)) +
  labs(x = "Number of UMIs", y = "Size Factor") +
  theme_bw() + 
  theme(text = element_text(size=20)) 


