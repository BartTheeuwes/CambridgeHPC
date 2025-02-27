
################
## Load atlas ##
################

# Load cell metadata
meta_atlas <- fread(args$atlas_metadata) %>%
  .[stripped==F & doublet==F & stage%in%args$atlas_stages]

# Filter
if (isTRUE(args$test)) meta_atlas <- head(meta_atlas,n=5000)

# Load SingleCellExperiment
# sce_atlas  <- readRDS(args$atlas_sce)[,meta_atlas$cell]
# sce_atlas <- logNormCounts(sce_atlas)
sce_atlas <- load_SingleCellExperiment(args$atlas_sce, normalise = TRUE, cells = meta_atlas$cell)

# Add cell metadata
colData(sce_atlas) <- meta_atlas %>% as.data.frame %>% tibble::column_to_rownames("cell") %>%
  .[colnames(sce_atlas),] %>% DataFrame()

if (!"cell"%in%colnames(colData(sce_atlas))) {
  sce_atlas$cell <- colnames(sce_atlas)
}


print("Atlas cell type diversity:")
print(table(meta_atlas$celltype))

################
## Load query ##
################

# Load cell metadata
meta_query <- fread(args$query_metadata) %>% .[pass_QC==T & sample%in%args$query_samples]
# if (isTRUE(args$test)) meta_query <- head(meta_query,n=1000)

# Load SingleCellExperiment
sce_query <- load_SingleCellExperiment(args$query_sce, cells = meta_query$cell)

if (!"cell"%in%colnames(colData(sce_query))) {
  sce_query$cell <- colnames(sce_query)
}

###################
## Sanity checks ##
###################



#############
## Prepare ## 
#############

# Filter out non-expressed genes
sce_query <- sce_query[rowSums(counts(sce_query))>10,]
sce_atlas <- sce_atlas[rowSums(counts(sce_atlas))>10,]

# Rename ensemble IDs to gene names
gene_metadata <- fread(io$gene_metadata) %>% .[,c("ens_id","symbol")] %>%
  .[symbol!="" & ens_id%in%rownames(sce_atlas)] %>%
  .[!duplicated(symbol)]

sce_atlas <- sce_atlas[rownames(sce_atlas)%in%gene_metadata$ens_id,]
foo <- gene_metadata$symbol; names(foo) <- gene_metadata$ens_id
rownames(sce_atlas) <- foo[rownames(sce_atlas)]

# Sanity cehcks
stopifnot(sum(is.na(rownames(sce_atlas)))==0)
stopifnot(sum(duplicated(rownames(sce_atlas)))==0)

# Intersect genes
genes.intersect <- intersect(rownames(sce_query), rownames(sce_atlas))

# Remove mitochondrial genes
genes.intersect <- genes.intersect[grep("mt-",genes.intersect,invert = T)]
genes.intersect <- genes.intersect[grep("Rik",genes.intersect,invert = T)]

# Subset SingleCellExperiment objects
sce_query  <- sce_query[genes.intersect,]
sce_atlas <- sce_atlas[genes.intersect,]

#################
## Subset HVGs ##
#################

# Select HVGs
# hvg <- getHVGs(sce.atlas, block=as.factor(sce.atlas$sample), p.value = 0.10)
# sce.query <- sce.query[hvg,]
# sce.atlas <- sce.atlas[hvg,]

#########################
## Subset marker genes ##
#########################

# # Load gene markers to be used as HVGs
# marker_genes.dt <- fread(io$marker_genes)
# marker_genes.dt <- marker_genes.dt[,head(.SD,n=50),by="celltype"]
# marker_genes <- unique(marker_genes.dt$ens_id)
# 
# stopifnot(all(marker_genes%in%rownames(sce.atlas)))
# stopifnot(all(marker_genes%in%rownames(sce.query)))
# 
# # Update SingleCellExperiment objects
# sce.query <- sce.query[marker_genes,]
# sce.atlas <- sce.atlas[marker_genes,]

