# +
# Swapped molecule removal

suppressPackageStartupMessages({
    library(DropletUtils)
    library(ggplot2)
    library(cowplot)
    library(Matrix)
})

### Change only this:
main = '/rds/project/rds-SDzz0CATGms/users/bt392/03_Stat3_RNA/'
subfolders = paste0(main, 'SLX-21143_fastq/SLX-21143_SITT', c('A2', 'A4', 'B2', 'B4', 'C4', 'D4', 'D5', 'E3','E5', 'F3', 'F5', 'G1', 'G3', 'G5', 'H1', 'H3'), '_HTJH3DSX2/')

out_dir <- paste0(main, '01_index_correction/')
out_subfolders = paste0(out_dir, 'SLX-21143_SITT', c('A2', 'A4', 'B2', 'B4', 'C4', 'D4', 'D5', 'E3','E5', 'F3', 'F5', 'G1', 'G3', 'G5', 'H1', 'H3'), '_HTJH3DSX2/')

nr_samples = 8

###

dir.create(out_dir, showWarnings = FALSE)
for(i in 1:length(out_subfolders)){
    dir.create(out_subfolders[i], showWarnings = FALSE)
}


mol_loc = paste0(subfolders, "outs/molecule_info.h5")

out_loc = paste0(out_subfolders, "matrix_unswapped.mtx")
bc_loc = paste0(out_subfolders, "barcodes_unswapped.tsv")
gene_loc = paste0(out_subfolders, "genes_unswapped.tsv")

unswapped = swappedDrops(mol_loc, get.swapped = TRUE)

ratios = sapply(1:length(unswapped$cleaned), function(i){
  sum(unswapped$swapped[[i]])/(sum(unswapped$cleaned[[i]]) + sum(unswapped$swapped[[i]]))
})

for(i in 1:length(mol_loc)){
  null_holder = writeMM(unswapped$cleaned[[i]], file = out_loc[i])
  write.table(colnames(unswapped$cleaned[[i]]), file = bc_loc[i], col.names = FALSE, row.names = FALSE, quote = FALSE)
  write.table(rownames(unswapped$cleaned[[i]]), file = gene_loc[i], col.names = FALSE, row.names = FALSE, quote = FALSE)
}





