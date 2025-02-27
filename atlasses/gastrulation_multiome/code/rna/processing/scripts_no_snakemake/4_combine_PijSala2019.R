library(scran)
library(scuttle)

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
} else {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
}

io$outfile.sce <- paste0(io$basedir,"/processed/rna/SingleCellExperiment_combined.rds")
io$outfile.metadata <- paste0(io$basedir,"/processed/rna/sample_metadata_combined.txt.gz")

###################
## Load Multiome ##
###################

cat("Loading multiome...\n")

sce.multiome <- readRDS(io$sce)
dim(sce.multiome)

############################
## Load PijSala2019 atlas ##
############################

cat("Loading atlas...\n")

sce.atlas <- readRDS(io$rna.atlas.sce)
dim(sce.atlas)

# Rename ensemble IDs to gene names
gene_metadata <- fread(io$gene_metadata) %>% .[,c("ens_id","symbol")] %>%
  .[symbol!="" & ens_id%in%rownames(sce.atlas)] %>%
  .[!duplicated(symbol)]

sce.atlas <- sce.atlas[rownames(sce.atlas)%in%gene_metadata$ens_id,]
foo <- gene_metadata$symbol; names(foo) <- gene_metadata$ens_id
rownames(sce.atlas) <- foo[rownames(sce.atlas)]

# Sanity cehcks
stopifnot(sum(is.na(rownames(sce.atlas)))==0)
stopifnot(sum(duplicated(rownames(sce.atlas)))==0)

#################
## Match genes ##
#################

genes <- intersect(rownames(sce.multiome), rownames(sce.atlas))

#######################
## Create merged SCE ##
#######################

cat("Merging...\n")

sce <- SingleCellExperiment(assays = list(counts = cbind(counts(sce.multiome[genes,]), counts(sce.atlas[genes,]))))
dim(sce)

###########################
## Merge sample metadata ##
###########################

metadata.multiome <- fread(io$metadata) %>%
    .[,c("cell", "sample", "stage", "nFeature_RNA", "nCount_RNA", "pass_rnaQC", "celltype.mapped", "celltype.score", "closest.cell", 
         "cxds_score", "cxds_call",  "bcds_score", "bcds_call", "hybrid_score", "hybrid_call", "doublet_call", 
        "TSSEnrichment_atac", "PromoterRatio_atac", "nFrags_atac", "pass_atacQC"
    )] %>%
  .[,dataset:="Multiome"] %>%
  setnames("celltype.mapped","celltype")

metadata.atlas <- fread(io$rna.atlas.metadata) %>%
  .[,c("cell","sample","stage","celltype","umapX","umapY","nFeature_RNA","nCount_RNA")] %>%
  .[,pass_rnaQC:=TRUE] %>%
  .[,dataset:="PijuanSala2019"]# %>%
  # .[,sample:=sprintf("%s_%s",stage,sample)]

metadata <-  plyr::rbind.fill(metadata.multiome, metadata.atlas)

# Add metadata to the ColData slot
colData(sce) <- metadata %>% as.data.frame %>% tibble::remove_rownames() %>%
  tibble::column_to_rownames("cell") %>% .[colnames(sce),] %>% DataFrame


##########################
## Compute size factors ##
##########################

cat("Computing size factors...\n")

clusts <- as.numeric(quickCluster(sce, method = "igraph", min.size = 100, BPPARAM = mcparam))
min.clust <- min(table(clusts))/2
new_sizes <- c(floor(min.clust/3), floor(min.clust/2), floor(min.clust))
sce <- computeSumFactors(sce, clusters = clusts, sizes = new_sizes, max.cluster.size = 3000)

###################
## Log Normalise ##
###################

if (args$normalise) {
  cat("Normalising...\n")
  sce <- logNormCounts(sce)
}


##########
## Save ##
##########

fwrite(metadata, io$outfile.metadata, sep="\t", quote=F, na="NA")
saveRDS(sce, io$outfile.sce)
