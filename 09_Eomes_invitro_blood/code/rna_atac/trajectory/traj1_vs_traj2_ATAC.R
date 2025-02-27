## Load packages
suppressPackageStartupMessages(library(destiny))
suppressPackageStartupMessages(library(tradeSeq))
suppressPackageStartupMessages(library(slingshot))
suppressPackageStartupMessages(library(BiocParallel))

here::i_am("rna_atac/trajectory/traj1_vs_traj2.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

# Options
BPPARAM <- BiocParallel::bpparam()
BPPARAM$workers = 32
set.seed(6)

## I/O
args = list()
# Metadata
args$metadata = file.path(io$basedir, 'results/rna_atac/clustering/metadata_celltype_annotated.txt.gz')
# atac sce
args$peakmtx = file.path(io$basedir, 'processed/atac/archR/Matrices/PeakMatrix_summarized_experiment.rds')
# outdir
args$outdir = file.path(io$basedir, 'results/rna_atac/trajectory/')
# Slingshot
args$slingshot = file.path(args$outdir, 'slingshot.rds')
# peak gene linkage
args$peak2gene = file.path(args$outdir, 'acc_gene_correlation.txt.gz')

## Load data
print('loading data')
# Load meta
meta = fread(args$metadata)[day%in%c('D3.5', 'D4', 'D4.5', 'D5')] # & genotype=='WT']
# load sce
atac.sce = readRDS(args$peakmtx)[,meta$cell]
colData(atac.sce) = meta %>% as.data.frame() %>% tibble::column_to_rownames('cell') %>% DataFrame()

# load slingshot
slingshot = readRDS(args$slingshot)
pseudotime = slingPseudotime(slingshot, na=F)

## Subset Cells
print('Subsetting cells')
# Keep only WT cells
sce_filt = atac.sce[,colData(atac.sce)$genotype == 'WT']
slingshot = slingshot[colnames(sce_filt),]
pseudotime = pseudotime[colnames(sce_filt),]

# Only keep YS vs Allantois trajectory
slingshot = slingshot[,c(2,3)]
pseudotime = pseudotime[,c(2,3)]

# Subset to cells within one of the trajectories
cells_keep = ifelse(pathStats(slingshot)$weights[,1] == 0 & pathStats(slingshot)$weights[,2] == 0, FALSE, TRUE)
slingshot = slingshot[cells_keep,]
sce_filt = sce_filt[, cells_keep]
pseudotime = pseudotime[cells_keep,]

# Adjust pseudotime to match between the two lineages
pseudodifference = as.data.table(slingPseudotime(slingshot)) %>% .[,diff:=Lineage3-Lineage2]
pseudotime[,2] = pseudotime[,2] - mean(pseudodifference$diff, na.rm=T)

# Subset by pseudotime
min_pseudo = 5
max_pseudo = 10
within_pseudotime_cells = as.data.table(pseudotime, keep.rownames=T) %>% 
    .[(Lineage2 >= min_pseudo & Lineage3 >= min_pseudo) | (Lineage2 >= min_pseudo & is.na(Lineage3)) | (is.na(Lineage2) & Lineage3 >= min_pseudo)] %>% 
    .[(Lineage2 <= max_pseudo & Lineage3 <= max_pseudo) | (Lineage2 <= max_pseudo & is.na(Lineage3)) | (is.na(Lineage2) & Lineage3 <= max_pseudo)]

slingshot = slingshot[within_pseudotime_cells$rn,]
sce_filt = sce_filt[, within_pseudotime_cells$rn]
pseudotime = pseudotime[within_pseudotime_cells$rn,]

# Minmax normalise both pseudotimes so they range from 0 to 1
minmax = function(x){(x-min(x))/(max(x)-min(x))}
pseudotime[,1] = minmax(pseudotime[,1])
pseudotime[,2] = minmax(pseudotime[,2])

sce_filt$slingshot = slingshot

## Subset peaks to those significantly linking with genes
print('Subsetting genes')
# load peak2gene linkage
p2g = fread(args$peak2gene)
# Filter to only significant peaks
sce_filt = sce_filt[unique(p2g[pval<0.05, peak]), ]

# Test if same cells object
summary(rownames(slingshot) == colnames(sce_filt))


# Run trade-seq
print('Run Trade-seq')
U <- model.matrix(~sce_filt$replicate)

sceGAM <- fitGAM(counts = assay(sce_filt),
                 pseudotime = pseudotime,
                 cellWeights =pathStats(slingshot)$weights,
                 U = U,
                 nknots=8, 
                 parallel = TRUE,
                 BPPARAM=BPPARAM,
                 verbose=T)

print('Saving Trade-seq results')
saveRDS(sceGAM, paste0(args$outdir,"/sceGAM_WT_Traj1_vs_Traj2_ATAC.rds"))