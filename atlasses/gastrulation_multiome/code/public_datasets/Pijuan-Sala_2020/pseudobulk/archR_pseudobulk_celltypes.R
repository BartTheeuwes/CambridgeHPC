
########################
## Load ArchR project ##
########################

source("/Users/ricard/gastrulation_multiome_10x/public_datasets/Pijuan-Sala_2020/archR/load_archR_project.R")

#####################
## Define settings ##
#####################

io$outdir <- paste0(io$basedir,"/data/processed/archR/pseudobulk_matrices")

################
## Pseudobulk ##
################

# https://www.ArchRProject.com/bookdown/how-does-archr-make-pseudo-bulk-replicates.html

# This function will merge cells within each designated cell group for the generation of pseudo-bulk replicates 
# and then merge these replicates into a single insertion coverage file.
# Output: creates files in archR/GroupCoverages/celltype: [X]._.Rep[Y].insertions.coverage.h5
ArchRProject <- addGroupCoverages(ArchRProject, groupBy = "celltype")

###################
## Sanity checks ##
###################

# getGroupBW(ArchRProject, groupBy = "celltype", tileSize = 500)

###################################################
## Pseudobulk into a SummarizedExperiment object ##
###################################################

se_list <- list()
matrices <- getAvailableMatrices(ArchRProject)
for (i in matrices) {
  
  # summarise
  se_list[[i]] <- getGroupSE(ArchRProject, groupBy = "celltype", useMatrix = i)
  
  # save
  outfile <- sprintf("%s/pseudobulk_%s_summarized_experiment.rds",io$outdir,i)
  saveRDS(se_list[[i]], outfile)
}


##########
## Plot ##
##########

for (i in matrices) {
  m <- as.matrix(assay(se_list[[i]],i))
  
  # PCA
  pca <- prcomp(t(m))$x[,1:10]
  
  # correlation
  r <- cor(t(pca),t(pca))
  
  # plot heatmap
  pdf(sprintf("%s/pseudobulk_%s_similarity.pdf",io$outdir,i), width=8, height=6)
  pheatmap::pheatmap(r)
  dev.off()
  
}

