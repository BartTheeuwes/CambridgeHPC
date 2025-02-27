# preparing extended atlas

suppressPackageStartupMessages({
    library(data.table)
    library(Matrix)
    library(dplyr) 
    library(scran)
    library(scater)
    library(igraph)
    library(BiocParallel)
    library(purrr)
    library(stringr)
})

ncores = 16
mcparam = SnowParam(workers = ncores)
register(mcparam)

main = '/rds/project/rds-SDzz0CATGms/users/bt392/atlasses/extended/'
out = file.path(main, 'blood/')

ct_keep = c('Haematoendothelial progenitors',
  'Blood progenitors',
  'Embryo proper endothelium',
  'Venous endothelium',
  'EMP',
  'Erythroid',
  'MEP',
  'Megakaryocyte progenitors',
  'Chorioallantoic-derived erythroid progenitors',
  'YS mesothelium-derived endothelial progenitors')

print('reading meta')
meta = readRDS(paste0(main, 'integrated_meta_celltype_clus.rds')) %>% as.data.table(keep.rownames=T) %>% setnames('rn', 'cell')
meta = meta[, c('cell', 'sample', 'stage', 'celltype.descendant', 'celltype.clustering')]
colnames(meta) = c('cell', 'sample', 'stage', 'celltype_original', 'celltype')

meta = meta[celltype %in% ct_keep] %>%
    .[,dataset:=str_split(as.character(cell),'_') %>%  map_chr(1)] %>%
    .[,dataset:=ifelse(dataset=='ext', 'ext', 'blanca')]

print('writing meta')
fwrite(meta, paste0(out, 'meta_blood.txt.gz'), sep='\t',row.names = FALSE)

print('reading counts')
counts = readMM(file.path(main, 'integrated_counts_sparse.mtx'))
counts = counts[rowSums(counts) > 100, meta$cell]

print('writing genes')
fwrite(rownames(counts), paste0(out, 'genes.txt.gz'), sep='\t',row.names = FALSE)

print('writing counts')
counts = as(as.matrix(counts[,-1]), "dgTMatrix")
writeMM(counts, paste0(out, 'raw_counts.mtx'))

sce = SingleCellExperiment(assays = list("counts" = counts))

lib.sizes = Matrix::colSums(counts(sce))
sce = sce[calculateAverage(sce)>0.1,]

coldata(sce) = meta

clusts = as.numeric(quickCluster(sce, method = "igraph", min.size = 100, BPPARAM = mcparam))

#now run the normalisation
#number of cells in each cluster should be at least twice that of the largest 'sizes'
min.clust = min(table(clusts))/2
new_sizes = c(floor(min.clust/3), floor(min.clust/2), floor(min.clust))
sce = computeSumFactors(sce, clusters = clusts, sizes = new_sizes, max.cluster.size = 3000)

fwrite(sizeFactors(sce), quote = F, col.names = F, row.names = F, file = paste0(out, "sizefactors.txt.gz"))

saveRDS(sce, file.path(out, 'sce_blood.rds'))


