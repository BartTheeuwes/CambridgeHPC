# +
# Load packages
suppressPackageStartupMessages({
    library(Seurat)
    library(scran)
    library(scater)
    library(batchelor)
    library(SingleCellExperiment)
    library(scater)
    library(Matrix)
    library(reshape2)
    library(data.table)
    library(BiocParallel)
    library(dplyr)
    library(gridExtra)
})

    ncores = 8
    mcparam = MulticoreParam(workers = ncores)
    register(mcparam)
    BPPARAM = SerialParam()

options(repr.plot.width=15, repr.plot.height=8)


source("/rds/project/rds-SDzz0CATGms/users/bt392/atlasses/Erythroid/mapping_functions.R")

out_dir = '/rds/project/rds-SDzz0CATGms/users/bt392/atlasses/Erythroid/'

# load cao
cao_in = '/rds/project/rds-SDzz0CATGms/users/bt392/atlasses/schendure/'

counts = readMM(paste0(cao_in, 'raw_counts.mtx')) # readMM for .mtx input 
genes = read.csv(paste0(cao_in, 'genes.csv'), stringsAsFactors = F)
cao_meta = read.table(paste0(cao_in, 'meta.csv'), header = TRUE, sep = ",", stringsAsFactors = FALSE, comment.char = "$")

rownames(counts) = genes[,3] #ensembl
cao_meta$cell = paste0('cao_', cao_meta$cell)
colnames(counts) = cao_meta$cell

cao_sce = SingleCellExperiment(assays = list("counts" = counts))
cao_sce = cao_sce[Matrix::rowSums(counts(cao_sce)) > 0,]

genes = read.csv(paste0(cao_in, 'genes.csv'), stringsAsFactors = F)

cao_meta$sample = cao_meta$development_stage
cao_meta$stage = paste0('E', cao_meta$development_stage)
cao_meta$development_stage = NULL
cao_meta$sample = paste0('cao_', cao_meta$sample)
cao_meta$exp = 'cao'


sizeFactors(cao_sce) = cao_meta$Size_Factor

cao_meta = cao_meta %>% select(cell, sample, stage, Main_cell_type, exp)
colnames(cao_meta) = c('cell', 'sample', 'stage', 'celltype', 'exp')

# load pijuan-sala atlas
atlas_in = '/rds/project/rds-SDzz0CATGms/users/bt392/atlasses/extended/'

counts = readMM(paste0(atlas_in, 'raw_counts.mtx')) # readMM for .mtx input 
genes = read.csv(paste0(atlas_in, 'genes_name.tsv'), stringsAsFactors = F, header=TRUE)

pijuan_meta = read.table(paste0(atlas_in, 'meta.tab'), header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "$")
pijuan_meta = pijuan_meta[, c('cell', 'sample', 'stage', 'celltype')]
pijuan_meta$cell = paste0('pijuan_', pijuan_meta$cell)
pijuan_meta$sample = paste0('pijuan_', pijuan_meta$sample)
pijuan_meta$exp = 'pijuan'


rownames(counts) = genes[,2] #ensembl
colnames(counts) = pijuan_meta$cell

pijuan_sce = SingleCellExperiment(assays = list("counts" = counts))
pijuan_sce = pijuan_sce[Matrix::rowSums(counts(pijuan_sce)) > 0,]
rm(counts)

sfs = read.table(paste0(atlas_in, 'sizefactors.tab'), stringsAsFactors = F)[,1]
sizeFactors(pijuan_sce) = sfs

keep = c('Blood progenitors',
            'Haematoendothelial progenitors', 
            'Embryo proper endothelium',
            'Intermediate mesoderm',
            'Nascent mesoderm',
            'Somitic mesoderm',
            'Presomitic mesoderm',
            'NMPs/Mesoderm-biased',
            'Allantois endothelium',
            'Erythroid',
            'Somitic mesoderm',
            'YS endothelium',
            'Venous endothelium',
            'MEP',
            'Megakaryocyte progenitors', 
            'Chorioallantoic-derived erythroid progenitors',
            'YS mesothelium-derived endothelial progenitors'
        )

pijuan_meta = pijuan_meta[pijuan_meta$celltype %in% keep,]
pijuan_sce = pijuan_sce[, pijuan_meta$cell]

# genes shared across datasets
shared_genes = merge(data.frame('V1'=names(cao_sce)), data.frame('V1'=names(pijuan_sce)), by='V1')

# filter cells & genes chimera dataset
cao_sce = cao_sce[shared_genes$V1,]

# filter cells & genes atlas dataset
pijuan_sce = pijuan_sce[shared_genes$V1,]

atlas_meta = rbind(pijuan_meta, cao_meta)



### START BATCH CORRECTION ###
npcs = 50

message("Normalizing joint dataset...")

#easier to avoid directly binding sce objects as it is a lot more likely to have issues
big_sce <- multiBatchNorm(cbind(pijuan_sce, cao_sce),
                        batch=c(atlas_meta$exp))

#big_sce <- scater::normalize(sce_all)
#big_sce <- scater::logNormCounts(sce_all) # edited 09.02. because normalize deprecated in favour of logNormCounts
big_sce <- multiBatchNorm(big_sce, batch=atlas_meta$sample) # edited 17.02. because now multibatchnorm exists
message("Done\n")

hvgs <- getHVGs(big_sce, block=atlas_meta$sample)
message("Done\n")

message("Performing PCA...")
atlas_pca <- multiBatchPCA(big_sce,
                       batch=atlas_meta$sample,
                       subset.row = hvgs,
                       d = npcs,
                       preserve.single = TRUE,
                       assay.type = "logcounts")[[1]]
rownames(atlas_pca) <- colnames(big_sce) 
message("Done\n")

message("Batch effect correction for the atlas...")  
order_df        <- atlas_meta[!duplicated(atlas_meta$sample), c("stage", "sample")]
order_df$ncells <- sapply(order_df$sample, function(x) sum(atlas_meta$sample == x))

order_df$stage  <- factor(order_df$stage, 
                    levels = rev(c("E13.5",
                                   "E12.5",
                                   "E11.5",
                                   "E10.5",
                                   "E9.5",
                                   "E9.25",
                                   "E9.0",
                                   "E8.75",
                                   "E8.5",
                                   "E8.25",
                                   "E8.0",
                                   "E7.75",
                                   "E7.5",
                                   "E7.25",
                                   "mixed_gastrulation",
                                   "E7.0",
                                   "E6.75",
                                   "E6.5")))

order_df       <- order_df[order(order_df$stage, order_df$ncells, decreasing = TRUE),]
order_df$stage <- as.character(order_df$stage)

set.seed(42)
atlas_corrected <- doBatchCorrect(counts         = logcounts(big_sce[hvgs,]), 
                                timepoints      = atlas_meta$stage, 
                                samples         = atlas_meta$sample, 
                                timepoint_order = order_df$stage, 
                                sample_order    = order_df$sample, 
                                pc_override     = atlas_pca,
                                npc             = npcs)
message("Done\n")

message("Correct between datasets...")                        
correct <- reducedMNN(atlas_corrected,
                  batch=atlas_meta$exp)$corrected
# -

write.csv(atlas_corrected, file =  paste0(out_dir,"uncorrected_PCs.csv"),row.names=TRUE)
write.csv(correct, file =  paste0(out_dir,"corrected_PCs.csv"),row.names=TRUE)
saveRDS(big_sce, file =  paste0(out_dir,"erythroid_atlas.RDS"))
write.csv(atlas_meta, file =  paste0(out_dir,"erythroid_atlas_meta.csv"),row.names=FALSE)
