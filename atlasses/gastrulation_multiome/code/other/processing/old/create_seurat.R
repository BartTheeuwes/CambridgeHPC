
#####################
## Define settings ##
#####################

source("/Users/ricard/gastrulation_multiome_10x/settings.R")

# io$outfile <- io$seurat
io$outfile <- "/Users/ricard/data/gastrulation_multiome_10x/multiome2/processed/seurat.rds"

###############
## Load data ##
###############

counts <- Read10X(c(
  # "multiome1" = paste0(io$basedir,"/multiome1/original/filtered_feature_bc_matrix"),
  "multiome2" = paste0(io$basedir,"/multiome2/original/filtered_feature_bc_matrix")
), strip.suffix = TRUE)
lapply(counts,dim)

###################
## Create Seurat ##
###################

seurat <- CreateSeuratObject(
  counts = counts["Gene Expression"][[1]],
  project = "Gastrulation Multiome 10x",
  min.cells = 1
)

# Add ATAC modality
seurat[["ATAC"]] <- CreateAssayObject(counts = counts["Peaks"][[1]])

seurat

##################
## Add metadata ##
##################

# metadata <- fread("/Users/ricard/data/10x_rna_atac/sample_metadata.csv") %>%
#   .[,barcode:=gsub("-1","",barcode)]

# dt <- data.table(barcode=colnames(seurat)) %>%
  # merge(metadata,by="barcode", all.x=TRUE) %>%
  # .[,c("pass_rnaQC","pass_accQC"):=FALSE] %>%
  # .[!is.na(celltype),c("pass_rnaQC","pass_accQC"):=TRUE] %>%
  # tibble::column_to_rownames("barcode")

# seurat <- AddMetaData(seurat, dt)

head(seurat@meta.data)

# metadata <- seurat@meta.data %>%
#   tibble::rownames_to_column("barcode") %>%
#   as.data.table %>%
#   .[,orig.ident:=NULL]
# fwrite(metadata, io$metadata, sep="\t", quote=F)

##########
## Save ##
##########

saveRDS(seurat, io$outfile, compress = FALSE)
