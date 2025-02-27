here::i_am("rna_atac/clustering/05_get_atlas_markers.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(edgeR))

BPPARAM <- BiocParallel::bpparam()
BPPARAM$workers = 21

# Multi core using future - built in to seurat
plan("multicore", workers = 16)
options(future.globals.maxSize = 50 * 1024 ^ 3) # for 50 Gb RAM

# outdir
args = list()
args$outdir = file.path(io$basedir, 'results/rna_atac/celltype_markers/')
dir.create(args$outdir, recursive=TRUE, showWarnings =FALSE)

# Load cell metadata
# Update cell-types --> only change ExE mesoderm to LPM for this comparison
meta_atlas = fread(paste0(io$basedir, '/../atlasses/gastrulation_multiome/results/atac/archR/qc/sample_metadata_after_qc.txt.gz')) %>% 
    .[pass_atacQC == T & doublet_call == F & genotype == 'WT'] %>% 
    .[,celltype_updated := celltype] %>% 
    .[,celltype_updated := ifelse(celltype_updated == 'ExE_mesoderm', 'Lateral_plate_mesoderm', celltype_updated)]

celltypes_keep = meta_atlas$celltype_updated
meta_atlas = meta_atlas[celltype_updated %in% celltypes_keep]

# Load Atlas

vivo_peaks = readRDS(paste0(io$basedir, '/../atlasses/gastrulation_multiome/data/processed/atac/archR/Matrices/PeakMatrix_summarized_experiment.rds'))[,meta_atlas$cell]

# Marker gene identification using edgeR (https://www.bioconductor.org/packages/devel/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf)
library(BiocParallel)
BPPARAM = MulticoreParam(22)
vivo_peaks = as(vivo_peaks, 'SingleCellExperiment')
atlas.pb = aggregateAcrossCells(vivo_peaks, id = meta_atlas[,c("celltype_updated", "sample")], BPPARAM = BPPARAM)

# discard samples with too few cells
discarded <- atlas.pb$ncells < 15
atlas.pb <- atlas.pb[,!discarded]
atlas.pb$sample = factor(atlas.pb$sample)

# Convert to edgeR object
atlas.edgeR = scran::convertTo(atlas.pb, type="edgeR")

keep.genes = filterByExpr(atlas.edgeR, group=atlas.edgeR$samples$celltype_updated, min.count=10, min.total.count=50)
atlas.edgeR <- atlas.edgeR[keep.genes, , keep=FALSE]
table(keep.genes)
atlas.edgeR = calcNormFactors(atlas.edgeR)
design = model.matrix(~ factor(atlas.edgeR$samples$celltype_updated) + factor(atlas.edgeR$samples$sample))

# Estimate dispersion
atlas.edgeR <- estimateDisp(atlas.edgeR, design, robust=TRUE)

fit <- glmQLFit(atlas.edgeR, design, robust=TRUE)
saveRDS(fit, file.path(args$outdir, 'Atlas_ATAC_edgeR_fit.rds'))
#fit = readRDS(file.path(args$outdir, 'Atlas_edgeR_fit.rds'))

cluster = factor(atlas.edgeR$samples$celltype_updated)
ncls = nlevels(cluster)
contr = rbind( matrix(1/(1-ncls), ncls, ncls),
matrix(0, ncol(design)-ncls, ncls) )
diag(contr) = 1
contr[1,] = 0
rownames(contr) = colnames(design)
colnames(contr) = paste0("cluster", levels(cluster))

# Get markers for every cell-type
markers = mclapply(1:ncls, function(i){
    a = Sys.time()
    tmp = glmQLFTest(fit, contrast=contr[,i])
    tmp = as.data.table(topTags(tmp, n=1e10)$table, keep.rownames=T) %>% 
        setnames('rn', 'gene') %>% 
        .[,celltype := levels(cluster)[i]]
    b = Sys.time()
    print(paste0('Cluster ', i, ' done, duration: ', b-a))
    return(tmp)
}, mc.cores=24)

saveRDS(markers, file.path(args$outdir, 'Atlas_ATAC_markers_edgeR.rds'))
