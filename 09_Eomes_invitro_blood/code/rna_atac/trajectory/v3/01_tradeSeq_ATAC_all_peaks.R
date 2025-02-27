here::i_am("rna_atac/trajectory/01_trajectory_inference_v2.ipynb")

source(here::here("settings.R"))
source(here::here("utils.R"))

suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(destiny))
suppressPackageStartupMessages(library(tradeSeq))
suppressPackageStartupMessages(library(slingshot))
suppressPackageStartupMessages(library(BiocParallel))

BPPARAM <- BiocParallel::bpparam()
BPPARAM$workers = 12

# Multi core using future - built in to seurat
plan("multicore", workers = 16)
options(future.globals.maxSize = 50 * 1024 ^ 3) # for 50 Gb RAM
set.seed(6)

# I/O
args = list()

# atac sce
args$peakmtx = file.path(io$basedir, 'processed/atac/archR/Matrices/PeakMatrix_summarized_experiment.rds')

# outdir
args$outdir = file.path(io$basedir, 'results/rna_atac/trajectory/v3')
dir.create(args$outdir, recursive=TRUE, showWarnings =FALSE)

print('Loading objects')
R_script_list = readRDS(file.path(args$outdir, 'ATAC_tradeseq_objects.rds'))

print('loading peak mtx')
peaks.sce = readRDS(args$peakmtx)[,R_script_list$cells]

peaks.mtx = peaks.sce@assays@data$PeakMatrix
rownames(peaks.mtx) = R_script_list$peaks
colnames(peaks.mtx) = R_script_list$cells

# Get model matrix
U <- model.matrix(~R_script_list$replicate)

print('Running fitGAM')
sceGAM <- fitGAM(counts = peaks.mtx,
              #   sds = rna.sce$slingshot,
                 pseudotime = R_script_list$pseudotime,
                 cellWeights = R_script_list$weight,
                 U = U,
                 nknots=8, 
            #     genes = unique(R_script_list$peak_gene_correlation[padj<0.05][,1][[1]]), # Keep only peaks that are linked to the differential genes
                 genes = R_script_list$peak_stats[FDR<0.05, rn],
                 parallel = TRUE,
                 BPPARAM=BPPARAM,
                 verbose=T) 

saveRDS(sceGAM, file.path(args$outdir, 'sceGAM_ATAC_all_peaks.rds'))