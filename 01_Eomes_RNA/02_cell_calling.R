# Load packages
suppressPackageStartupMessages({
    library(DropletUtils)
    library(ggplot2)
    library(Matrix)
    library(BiocParallel)
    library(knitr)
    library(reshape2)
})

ncores = 3
mcparam = MulticoreParam(workers = ncores)
register(mcparam)


##===Read in files===##
main = '/rds/project/rds-SDzz0CATGms/users/bt392/01_Eomes_RNA/'
subfolders = paste0(main, '01_index_correction/SLX-20795_SITT', c('A11', 'A12', 'B11', 'F11', 'G10', 'G11', 'H10', 'H11'), '_HKTG2DRXY/')

out_dir <- paste0(main, '02_calledCells/')
dir.create(out_dir, showWarnings = FALSE)

mtx_loc <- paste0(subfolders, "matrix_unswapped.mtx")
bc_loc <- paste0(subfolders, "barcodes_unswapped.tsv")
gene_loc <- paste0(subfolders, "genes_unswapped.tsv")

# read matrices and tables into a single variable
matrices = bplapply(mtx_loc, readMM)
bcs = bplapply(bc_loc, function(x) read.table(x, header = FALSE, stringsAsFactors = FALSE)[,1])


#correct barcode sample number
for(i in 1:length(bcs)){
  bcs[[i]] = paste0(bcs[[i]], "-", i)
}


##==Open dataframe with the sample metadata/parameters==##
exp_design = read.csv(paste0(main, 'SLX-20795_fastq/SLX-20795.HKTG2DRXY.s_1.contents.csv'))
exp_design$stage = 'E8.5'
exp_design$batch = 1:nrow(exp_design)
exp_design$batch_name = paste0(exp_design$Pool, '_', exp_design$Barcode, '_HKTG2DRXY')
exp_design$tdTom = rep(c('positive','negative'), 4)

exp_design


# ### ##==Do cell calling==##

set.seed(42)
#do call
outs = lapply(matrices, emptyDrops, niters = 20000, ignore = 4999, BPPARAM = mcparam, lower = 100, retain = Inf)


#identify cells
sigs = lapply(outs, function(x) x$FDR <= 0.01 & !is.na(x$FDR))

#subset the cells
cells = lapply(1:length(matrices), function(i) matrices[[i]][, sigs[[i]]])
barcodes = lapply(1:length(bcs), function(i) bcs[[i]][sigs[[i]]])

#append
counts = do.call(cbind, cells)
barcodes = do.call(c, barcodes)


##==save==##
writeMM(counts, file = paste0(out_dir,"raw_counts.mtx"))
write.table(barcodes, file = paste0(out_dir, "barcodes.tsv"), col.names = FALSE, row.names = FALSE, quote = FALSE)
file.copy(from = paste0(subfolders[1], 'genes_unswapped.tsv'),
          to = paste0(out_dir, "genes.tsv"), overwrite = TRUE)
write.csv2(exp_design, paste0(main, 'experimental_design_meta_eomes.csv'), row.names = FALSE)

sessionInfo()
