## Subset data from Cao et al. to only haematopoietic cells
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(R.utils)
  library(Matrix)
})

main = '/rds/project/rds-SDzz0CATGms/users/bt392/atlasses/schendure/'

# Read metadata and select cells with haem related identity
meta = read.csv2(paste0(main, 'cell_annotate.csv'), header=TRUE, sep=',') %>% select(sample, development_stage, Size_Factor,detected_doublet, Main_cell_type)
clusters_keep = c('Primitive erythroid lineage',
                  'Definitive erythroid lineage', 
                  'Megakaryocytes', 
                  'Endothelial cells')
meta = meta %>% mutate(cell_id = row_number()) %>% filter(Main_cell_type %in% clusters_keep & detected_doublet == FALSE)
nrow(meta)

# 230k cells are way too much, so randomly subset it!
keep = lapply(unique(meta$development_stage), function(x){
    if(sum(meta$development_stage == x) < 5000) {
    return(which(meta$development_stage == x))
  } else {
    hits = which(meta$development_stage == x)
    return(sample(hits, 5000))
  }
})
keep = do.call(c, keep)

meta = meta[keep,]

# Select genes
genes = read.csv2(paste0(main, 'GSE119945_gene_annotate.csv'), header=TRUE, sep=',') 
genes = genes %>% mutate(gene_id = row_number()) %>%
    filter(gene_type %in% c('protein_coding'))
genes = genes %>%  filter(!gene_id %in% genes[duplicated(genes$gene_short_name),]$gene_id)
nrow(genes)

# +
# open count file
count = fread(paste0(main, 'GSE119945_gene_count.txt'))[-1,] %>% select(1,2,3) 
colnames(count) = c('gene_id', 'cell_id', 'UMI')

# subset count file for cells & genes
count = count[cell_id %in% meta$cell_id & gene_id %in% genes$gene_id,]
count = dcast(count, gene_id ~ cell_id, value.var = 'UMI')

# NA -> 0
count[is.na(count)] <- 0

# Make sure genes and cells are in correct order
genes = genes[genes$gene_id %in% count$gene_id,]
genes = genes[match(genes$gene_id, count$gene_id),]
write.csv(genes, paste0(main, 'genes.csv'), row.names=FALSE)

meta = meta[meta$cell_id %in% colnames(count),]
meta = meta[match(meta$cell_id, meta$cell_id),]
write.csv(meta, paste0(main, 'meta.csv'), row.names=FALSE)

# Write new count file
print('writing counts')
count = as(as.matrix(count[,-1]), "dgTMatrix")
rownames(count) = genes$gene_short_name
colnames(count) = meta$cell_id
writeMM(count, paste0(main, 'raw_counts.mtx'))
