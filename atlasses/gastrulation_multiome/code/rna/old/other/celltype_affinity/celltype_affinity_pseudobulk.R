
##############
## Settings ##
##############

source("/Users/ricard/gastrulation_multiome_10x/settings.R")
io$outdir <- paste0(io$basedir,"/results/rna/celltype_affinity")

opts$stages <- c(
  "E8.5"
)

#####################
## Update metadata ##
#####################

# sample_metadata <- sample_metadata %>% 
#   .[pass_rnaQC==TRUE & hybrid_call_RNA==FALSE] %>%
#   .[celltype.mapped%in%opts$celltypes & stage%in%opts$stages]

# foo <- table(sample_metadata$celltype)<100
# if (any(foo)) {
#   warning("There are cell types with very small amount of cells, removing them:")
#   warning(paste(names(which(foo)),collapse=",  "))
#   sample_metadata <- sample_metadata[!celltype%in%names(which(foo))]
# }

###############
## Load data ##
###############

# Load gene markers
marker_genes.dt <- fread(io$atlas.marker_genes) %>%
  .[celltype%in%opts$celltypes]
length(unique(marker_genes.dt$celltype))

# Load average expression per celltype and gene
pseudobulk_expr.dt <- fread(io$average_expression_per_celltype) %>%
  .[ens_id%in%unique(marker_genes.dt$ens_id) & gene!=""] %>%
  .[celltype%in%opts$celltypes]

##################
## Computations ##
##################

# Calculate correlation coefficient between each pair of cell types, across genes 
m <- pseudobulk_expr.dt %>% 
  .[,id:=paste(ens_id,gene,sep="_")] %>%
  dcast(id~celltype, value.var="mean_expr") %>%
  matrix.please

r <- cor(m)
# diag(r) <- NA

##########
## Plot ##
##########

pdf(sprintf("%s/correlation_celltypes_heatmap.pdf",io$outdir), width = 11, height = 8)
pheatmap::pheatmap(r)
dev.off()

