library(SnapATAC)
library(purrr)
library(data.table)
library(viridisLite)
library(ggplot2)
library(GenomicRanges)

source(here::here("settings.R"))
source(here::here("atac/SnapATAC/snapatac_settings.R"))


######## needs to be run with R/3.5.2 ########

# script adapted from https://github.com/r3fang/SnapATAC/blob/master/examples/10X_brain_5k/README.md

### input/output ###


opts$cores               <- 8


bin_size <- 5000 # prefix to label plots



snap <- readRDS(snapio$rds_file)
snap

########################################
### Step 5. Dimensionality reduction ###
########################################

# We compute diffusion maps for dimentionality reduction.

snap <- runDiffusionMaps(obj=snap,
                         input.mat="bmat", 
                         num.eigs=50)

################################################
### Step 6. Determine significant components ###
################################################

# We next determine the number of reduced dimensions to include for downstream 
# analysis. We use an ad hoc method by simply looking at a pairwise plot and 
# select the number of dimensions in which the scatter plot starts looking like 
# a blob. In the below example, we choose the first 20 dimensions.

plotDimReductPW(obj=snap, 
                eigs.dims=1:50,
                point.size=0.3,
                point.color="grey",
                point.shape=19,
                point.alpha=0.6,
                down.sample=5000,
                pdf.file.name=paste0(snapio$plots_out, "/", bin_size ,"_dim_red_components.pdf"), 
                pdf.height=7, 
                pdf.width=7)

### Step 7. Graph-based clustering
# Using the selected significant dimensions, we next construct a K Nearest 
# Neighbor (KNN) Graph (k=15). Each cell is a node and the k-nearest neighbors 
# of each cell are identified according to the Euclidian distance and edges are 
# draw between neighbors in the graph.

snap <- runKNN(obj=snap,
               eigs.dims=1:18, 
               k=15)

snap <- runCluster(obj=snap,
                   tmp.folder=tempdir(),
                   louvain.lib="R-igraph",
                   seed.use=10)

snap@metaData$cluster <- snap@cluster

#############################
### Step 8. Visualization ###
#############################




snap <- runViz(obj=snap, 
               tmp.folder=tempdir(),
               dims=2,
               eigs.dims=1:18, 
               method="umap",
               num.cores = opts$cores,
               seed.use=10)

# save as rds

saveRDS(snap,  snapio$rds_file)

#par(mfrow = c(2, 2))

print("plotting clusters")

plotViz(obj=snap,
        method="umap", 
        main="Clusters",
        point.color=snap@cluster, 
        point.size=0.2, 
        point.shape=19, 
        point.alpha=0.8, 
        text.add=TRUE,
        text.size=1.5,
        text.color="black",
        text.halo.add=TRUE,
        text.halo.color="white",
        text.halo.width=0.2,
        down.sample=10000,
        legend.add=FALSE,
        pdf.file.name = paste0(snapio$plots_out, "/", bin_size , "_umap_clusters.pdf"),
        pdf.width = 7,
        pdf.height = 7)


print("plotting fragments")


plotFeatureSingle(obj=snap,
                  feature.value=log(snap@metaData[,"atac_fragments"]+1,10),
                  method="umap", 
                  main="Read Depth",
                  point.size=0.2, 
                  point.shape=19, 
                  down.sample=10000,
                  quantiles=c(0.01, 0.99),
                  pdf.file.name = paste0(snapio$plots_out, "/", bin_size , "_umap_fragments.pdf"),
                  pdf.width = 7,
                  pdf.height = 7)

print("plotting frip")


plotFeatureSingle(obj=snap,
                  feature.value=snap@metaData$atac_peak_region_fragments / snap@metaData$atac_fragments,
                  method="umap", 
                  main="FRiP",
                  point.size=0.2, 
                  point.shape=19, 
                  down.sample=10000,
                  pdf.file.name = paste0(snapio$plots_out, "/", bin_size , "_umap_frip.pdf"),
                  pdf.width = 7,
                  pdf.height = 7,
                  quantiles=c(0.01, 0.99)) # remove outliers)

print("plotting duplication")



plotFeatureSingle(obj=snap,
                  feature.value=snap@metaData$atac_dup_reads / snap@metaData$atac_raw_reads,
                  method="umap", 
                  main="Duplicates",
                  point.size=0.2, 
                  point.shape=19, 
                  down.sample=10000,
                  pdf.file.name = paste0(snapio$plots_out, "/", bin_size , "_umap_duplicates.pdf"),
                  pdf.width = 7,
                  pdf.height = 7,
                  quantiles=c(0.01, 0.99)) # remove outliers

print("plotting sample")



plotViz(obj=snap, 
        method="umap", 
        main = "Sample",
        point.color=snap@sample, 
        text.add=FALSE,
        point.size=0.2, 
        point.shape=19, 
        point.alpha=0.7,
        legend.add=TRUE,
        pdf.file.name = paste0(snapio$plots_out, "/", bin_size , "_umap_sample.pdf"))

print("plotting cell type")


plotViz(obj=snap,
        method="umap", 
        main="Cell type mapped from RNA-seq",
        point.color=snap@metaData$cell_type, 
        point.size=0.2, 
        point.shape=19, 
        point.alpha=0.7, 
        text.add=FALSE,
        legend.add = TRUE,
        legend.pos = "top",
        legend.text.size = 0.3,
        #down.sample = 1e4,
        pdf.file.name = paste0(snapio$plots_out, "/", bin_size , "_umap_celltype.pdf"))

print("plotting RNA QC")


plotViz(obj=snap,
        method="umap", 
        main="Cells passing QC in RNA-seq",
        point.color=snap@metaData$pass_rnaQC, 
        point.size=0.2, 
        point.shape=19, 
        point.alpha=0.7, 
        text.add=FALSE,
        legend.add = TRUE,
        legend.pos = "top",
        legend.text.size = 0.3,
        pdf.file.name = paste0(snapio$plots_out, "/", bin_size , "_umap_passRNAQC.pdf"))






