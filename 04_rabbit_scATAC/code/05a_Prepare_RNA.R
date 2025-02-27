here::i_am("01_create_arrow.R")
source(here::here("settings.R"))

geneAnnotation = readRDS(file.path(io$basedir, 'geneAnnotation.rds'))

geneAnnotation

rna_in = '/rds/project/rds-SDzz0CATGms/users/mlnt2/PhD_MT06/'
rna_out = file.path(io$basedir, 'RNA')
counts_in = file.path(rna_in, '3_cellqc/raw_counts.mtx')
genes_in = file.path(rna_in, '2_cellcalling/genes.tsv')
meta_in1 = file.path(rna_in, '5_doublet/meta.tab')
meta_in2 = file.path('/rds/project/rds-SDzz0CATGms/users/bt392/04_Rabbit_ATAC_2/RNA/metadata.csv')
sizefactors_in = file.path(rna_in, 'sizefactors.tab')

dir.create(file.path(rna_out), showWarnings = FALSE)


# Load data
counts = readMM(counts_in)
genes = fread(genes_in, header=FALSE)
meta = fread(meta_in1)

# Change col & row names
colnames(counts) = meta$cell
rownames(counts) = genes$V1

# Keep only cells that are not doublet or stripped nuclei
keep = meta[meta$doublet==FALSE & meta$stripped==FALSE,]$cell
counts = counts[,keep]

# load new metadata
celltypes = fread(file.path(rna_out, 'rna_meta.tsv')) %>%
     setnames(c('V1', 'celltype'), c('cell', 'celltype')) %>%
    .[order(match(cell, colnames(counts))),]
stages = fread(file.path(rna_out, 'metadata.csv')) %>%
    setnames('index', 'cell') %>%
     .[,c('cell', 'stage')] %>%
     .[order(match(cell, colnames(counts))),] 
meta = cbind(celltypes, stages[,2])

# Extract Granges of genes
rowranges = geneAnnotation$genes
genes_keep = genes[genes$V1 %in% rowranges$gene_id,]$V1 # keep only genes in ATAC data
counts = counts[genes_keep, ] 
rowranges = rowranges[rowranges$gene_id %in% rownames(counts),]
rowranges = rowranges[order(match(rowranges$gene_id, rownames(counts))),] # Match order between RNA & ATAC (not needed?)

# Create Ranged SummarizedExperiment object of scRNA-seq data
sceRNA <- SummarizedExperiment(assays=SimpleList(counts=counts),rowRanges=rowranges, colData=meta)
rownames(sceRNA) = sceRNA@rowRanges$symbol

# Save object
saveRDS(sceRNA, file.path(rna_out, 'RangedSummarizedExperiment.rds'))
