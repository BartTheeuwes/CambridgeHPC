here::i_am("rna_atac/differential/RNA.ipynb")

source(here::here("settings.R"))
source(here::here("utils.R"))

suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(destiny))
suppressPackageStartupMessages(library(destiny))
suppressPackageStartupMessages(library(tradeSeq))
suppressPackageStartupMessages(library(slingshot))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(BiocParallel))


BPPARAM <- BiocParallel::bpparam()
BPPARAM$workers = 70

# outdir
args = list()
args$groupA = 'WT'
args$groupB = 'KO'
opts$groups = c(args$groupA, args$groupB)

args$outdir = file.path(io$basedir, 'results/rna_atac/differential/Accessibility/')
dir.create(args$outdir, recursive=TRUE, showWarnings =FALSE)


sce_filt = readRDS(sprintf('%s/sce_al.rds', args$outdir))

# Set all cells to belong to one trajectory
cellWeights <- rep(1,ncol(sce_filt))

# create a model matrix -> Comparing across replicates
U <- model.matrix(~sce_filt$replicate)

sceGAM <- fitGAM(counts = assay(sce_filt, 'counts'), 
                 conditions = factor(sce_filt$genotype, levels=c('KO', 'WT')),
                 U = U,
                 pseudotime= round(sce_filt$pseudotime,2),
                 cellWeights=cellWeights,
                 nknots=3, 
                 parallel = TRUE,
                 BPPARAM=BPPARAM,
                 verbose=F)
sceGAM
print(sprintf('%s/sceGAM_v1.rds', args$outdir))
saveRDS(sceGAM, sprintf('%s/sceGAM_v1.rds', args$outdir))