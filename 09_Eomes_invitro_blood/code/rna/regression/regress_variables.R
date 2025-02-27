here::i_am("rna/regression/regress_variables.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(Seurat))

# Multicore
BPPARAM <- BiocParallel::bpparam()
BPPARAM$workers = 24

# Multi core using future - built in to seurat
plan("multicore", workers = 24)
options(future.globals.maxSize = 50 * 1024 ^ 3) # for 50 Gb RAM

# Set args
args = list()
args$rna_metadata = file.path(io$basedir, 'results/rna/doublet_detection/sample_metadata_after_doublets.txt.gz')
args$rna_sce = file.path(io$basedir, 'processed/rna/SingleCellExperiment.rds')
args$rna_nfeatures = 4000
args$vars_to_regress = c("nFeature_RNA", "nCount_RNA", "mitochondrial_percent_RNA", "ribosomal_percent_RNA")
args$regression_out = file.path(io$basedir, 'results/rna/regression/')
dir.create(args$regression_out, recursive=TRUE, showWarnings =FALSE)

# Load metadata
metadata_rna <- fread(args$rna_metadata) %>%
  .[pass_rnaQC==TRUE & doublet_call==FALSE] %>% 
  .[,exp := str_replace_all(sample, opts$sample2exp)]

# Load sce
rna.sce <- load_SingleCellExperiment(args$rna_sce, normalise = TRUE, cells = metadata_rna$cell)

# Add sample metadata to the colData of the SingleCellExperiment
colData(rna.sce) <- metadata_rna %>% as.data.frame %>% tibble::column_to_rownames("cell") %>%
  .[colnames(rna.sce),] %>% DataFrame()

# Filter features manually
rna.sce <- rna.sce[grep("*Rik|^Gm|^Mt-|^Rps|^Rpl|^Olfr",rownames(rna.sce), invert=T),]

# Remove genes with very low variance
gene_vars = rowVars(logcounts(rna.sce))
names(gene_vars) = rownames(rna.sce)
keep_genes = names(gene_vars[gene_vars>0.1])
rna.sce = rna.sce[keep_genes, ]

# Regress out variables 
logcounts_regressed.mtx <- RegressOutMatrix(
    mtx = logcounts(rna.sce),
    covariates = metadata_rna[,args$vars_to_regress,with=F]
  )


fwrite(as.data.table(logcounts_regressed.mtx, keep.rownames=T), file.path(args$regression_out,"logcounts_regressed_mtx.txt.gz"))