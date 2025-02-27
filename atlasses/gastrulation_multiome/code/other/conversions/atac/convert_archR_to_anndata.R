suppressPackageStartupMessages({
  library("reticulate")
})


#####################################
## Reticulate connection to scanpy ##
#####################################

use_python("/bi/group/reik/ricard/software/miniconda3/envs/main/bin/python")
sc <- import("scanpy")

#####################
## Define settings ##
#####################

source(here::here("settings.R"))
source(here::here("utils.R"))

# I/O
io$outfile <- file.path(io$basedir,"processed_new/atac/scanpy/anndata_atac_peak_matrix.h5ad")

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(io$metadata) %>% 
  .[pass_rnaQC==TRUE & pass_atacQC==TRUE & doublet_call==FALSE & !is.na(celltype.mapped_mnn)] %>%
  .[,c("cell", "sample", "stage", "celltype.mapped_mnn", "TSSEnrichment_atac", "nFrags_atac")] %>%
  setnames("celltype.mapped_mnn","celltype")

########################
## Load ArchR project ##
########################

source(here::here("atac/archR/load_archR_project.R"))

##################
## Subset ArchR ##
##################

ArchRProject.filt <- ArchRProject[sample_metadata$cell,]
ArchRProject.filt

###############
## Load data ##
###############

# Load PeakMatrix as a SummarizedExperiment object
atac.peak.se <- getMatrixFromProject(ArchRProject.filt, useMatrix="PeakMatrix", binarize = FALSE)
dim(atac.peak.se)

# Update peak metadata (rowData)
peak.df <- getPeakSet(ArchRProject.filt) %>% as.data.table %>% 
  .[,c("seqnames","start","end","score")] %>%
  setnames(c("chr","start","end","score")) %>%
  .[,idx:=sprintf("%s:%s-%s",chr,start,end)] %>%
  DataFrame
rownames(peak.df) <- sprintf("%s:%s-%s",peak.df$chr,peak.df$start,peak.df$end)

rownames(atac.peak.se) <- rownames(peak.df)
rowData(atac.peak.se) <- peak.df

# Update cell metadata (colData)
# colnames(atac.peak.se) <-  sample_metadata %>%
#   .[cell%in%colnames(atac.peak.se)] %>% setkey(cell) %>%
#   .[colnames(atac.peak.se)] %>% .$cell
# colData(atac.peak.se) <- sample_metadata %>% 
#   .[cell%in%colnames(atac.peak.se)] %>% setkey(cell) %>%
#   .[colnames(atac.peak.se)] %>% .[,cell:=NULL] %>%
#   tibble::column_to_rownames("cell") %>% DataFrame

#############################################
## Convert SingleCellExperiment to AnnData ##
#############################################

adata <- sc$AnnData(
    X   = t(assay(atac.peak.se)),
    obs = as.data.frame(colData(atac.peak.se)),
    var = as.data.frame(rowData(atac.peak.se))
)

head(adata$obs)
head(adata$var)

# Add cell type colors
# celltype.colors <- c('#532C8A', '#c19f70', '#f9decf', '#c9a997', '#B51D8D', '#3F84AA', '#9e6762', '#354E23', '#F397C0', '#ff891c', '#635547', '#C72228', '#f79083', '#EF4E22', '#989898', '#7F6874', '#8870ad', '#647a4f', '#EF5A9D', '#FBBE92', '#139992', '#cc7818', '#DFCDE4', '#8EC792', '#C594BF', '#C3C388', '#0F4A9C', '#FACB12', '#8DB5CE', '#1A1A1A', '#C9EBFB', '#DABE99', '#65A83E', '#005579', '#CDE088', '#f7f79e', '#F6BFCB')
# adata$uns$update(celltype.mapped_colors = celltype.colors)


##########
## Save ##
##########

adata$write_h5ad(io$outfile)
