###################
## Sanity checks ##
###################


###############################
## Load SingleCellExperiment ##
###############################

sce <- load_SingleCellExperiment(
  file = io$rna.sce, 
  cells = sample_metadata$cell, 
  normalise = TRUE, 
  remove_non_expressed_genes = FALSE
)

rownames(sce) <- toupper(rownames(sce))

##########################
## Load chromVAR scores ##
##########################

# Load SummarizedExperiment
atac.chromvar.se <- readRDS(sprintf("%s/deviations_summarized_experiment_%s.rds",io$archr.chromvar.dir,opts$motif_annotation))

# Subset cells
# atac.chromvar.se <- atac.chromvar.se[,colnames(atac.chromvar.se)%in%sample_metadata$cell]
atac.chromvar.se <- atac.chromvar.se[,sample_metadata$cell]

# Rename motifs
rownames(atac.chromvar.se) <- rowData(atac.chromvar.se)$name %>% toupper %>% stringr::str_split(.,"_") %>% map_chr(1)
rownames(atac.chromvar.se) <- gsub("TCFAP","TFAP",rownames(atac.chromvar.se))
rownames(atac.chromvar.se) <- gsub("NKX2","NKX2-",rownames(atac.chromvar.se))
rownames(atac.chromvar.se) <- gsub("NKX3","NKX3-",rownames(atac.chromvar.se))
rownames(atac.chromvar.se) <- gsub("NKX6","NKX6-",rownames(atac.chromvar.se))

# Remove duplicated motifs
atac.chromvar.se <- atac.chromvar.se[!duplicated(rownames(atac.chromvar.se)),]

################################
## Load motif2gene annotation ##
################################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/atac/archR/load_motif_annotation.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/atac/archR/load_motif_annotation.R")
} else {
  stop("Computer not recognised")
}

motif2gene.dt <- motif2gene.dt %>%
  .[gene%in%rownames(atac.chromvar.se) & gene%in%rownames(sce)] %>%
  .[,N:=length(unique(motif)),by="gene"] %>% .[N==1] %>% .[,N:=NULL]

sce <- sce[motif2gene.dt$gene,]
atac.chromvar.se <- atac.chromvar.se[motif2gene.dt$gene,]

#################
## Smooth data ##
#################

if (opts$denoise & opts$knn>1) {

  # trajectory.mtx <- trajectory.dt[,c("cell","V1")] %>% tibble::column_to_rownames("cell") %>% as.matrix
  
  # RNA
  pca.rna <- fread(io$pca.rna) %>% matrix.please %>% .[sample_metadata$cell,]
  rna.mtx <- smoother_aggregate_nearest_nb(mat=as.matrix(logcounts(sce)), D=pdist(pca.rna), k=opts$knn)
  
  # ATAC chromVAR
  pca.atac <- fread(io$pca.atac) %>% matrix.please %>% .[sample_metadata$cell,]
  atac.chromvar.mtx <- smoother_aggregate_nearest_nb(mat=as.matrix(assay(atac.chromvar.se,"z")), D=pdist(pca.atac), k=opts$knn)
  
} else {
  rna.mtx <- as.matrix(logcounts(sce))
  atac.chromvar.mtx <- as.matrix(assay(atac.chromvar.se,"z"))
}

colnames(rna.mtx) <- colnames(sce)
colnames(atac.chromvar.mtx) <- colnames(atac.chromvar.se)

################
## parse data ##
################

# chromvar_dt <- atac.chromvar.mtx %>% t %>% 
#   as.data.table(keep.rownames = T) %>%
#   setnames("rn","cell") %>% 
#   melt(id.vars=c("cell"), variable.name="gene", value.name="chromvar_zscore")

# rna_dt <- rna.mtx %>%
#   as.data.table(keep.rownames = T) %>%
#   setnames("rn","gene") %>%
#   melt(id.vars="gene", variable.name="cell", value.name="expr")
