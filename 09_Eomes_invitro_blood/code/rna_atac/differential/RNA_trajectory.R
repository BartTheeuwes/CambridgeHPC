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
BPPARAM$workers = 24

# outdir
args = list()
args$groupA = 'WT'
args$groupB = 'KO'
opts$groups = c(args$groupA, args$groupB)

args$outdir = file.path(io$basedir, 'results/rna_atac/differential/expression/')
dir.create(args$outdir, recursive=TRUE, showWarnings =FALSE)

sce_al = readRDS(sprintf('%s/sce_al.rds', args$outdir))

# Min values set very low to keep genes that are very specific
min.cdr = 0.05
min.expr = 0.05

expr.dt <- data.table(
    gene = rownames(sce_al),
    mean_groupA = rowMeans(logcounts(sce_al[,sce_al$genotype == args$groupA])) %>% round(2),
    mean_groupB = rowMeans(logcounts(sce_al[,sce_al$genotype == args$groupB])) %>% round(2),
    cdr_groupA = rowMeans(logcounts(sce_al[,sce_al$genotype == args$groupA])>0) %>% round(2),
    cdr_groupB = rowMeans(logcounts(sce_al[,sce_al$genotype == args$groupB])>0) %>% round(2)
)

#######################
## Feature selection ##
#######################
# filter genes
genes.to.use <- expr.dt[mean_groupA>=min.expr | mean_groupB>=min.expr,][cdr_groupA>=min.cdr | cdr_groupB>=min.cdr,gene]
genes.to.use = genes.to.use[!genes.to.use %in% genes.to.use[grep("*Rik|^Gm|^Rps|^Rpl|^Olfr", genes.to.use)]]  
length(genes.to.use)

# Subset genes
sce_filt = sce_al[genes.to.use,]

# Set all cells to belong to one trajectory
cellWeights <- rep(1,ncol(sce_filt))

# create a model matrix -> Comparing across replicates
U <- model.matrix(~sce_filt$replicate)

sceGAM <- fitGAM(counts = counts(sce_filt),
                 conditions = factor(sce_filt$genotype, levels=c('KO', 'WT')),
                 U = U,
                 pseudotime= sce_filt$pseudotime,
                 cellWeights=cellWeights,
                 nknots=6, 
                 parallel = TRUE,
                 BPPARAM=BPPARAM,
                 verbose=F)
sceGAM
print(sprintf('%s/sceGAM_v1.rds', args$outdir))
saveRDS(sceGAM, sprintf('%s/sceGAM_v1.rds', args$outdir))