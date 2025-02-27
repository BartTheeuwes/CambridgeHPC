library(muscat)
# library(DESeq2)
library(scuttle)

#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
} else {
  stop("Computer not recognised")
}

# I/O
# io$cell2metacell <- "/Users/ricard/data/gastrulation_multiome_10x/results/rna/metacells/cell2metacell_1000metacells.txt.gz"
io$cell2metacell <- paste0(io$basedir,"/results/rna_atac/metacells/cell2metacell_2500metacells.txt.gz")
io$outfile <- paste0(io$basedir,"/results/rna_atac/metacells/SingleCellExperiment_2500metacells.rds")

########################
## Load cell2metacell ##
########################

cell2metacell <- fread(io$cell2metacell) 

# cell2metacell[,.N,by="Metacell"]

###############
## Load data ##
###############

# Load cell metadata
sample_metadata <- fread(io$metadata) %>%
  .[cell%in%cell2metacell$cell] %>%
  merge(cell2metacell,"cell")

# Load SingleCellExperiment
sce <- load_SingleCellExperiment(file=io$rna.sce, cells=sample_metadata$cell)
colData(sce) <- sample_metadata %>% tibble::column_to_rownames("cell") %>% DataFrame

###################################
## Aggregate counts per celltype ##
###################################

# assays(sce)$cpm <- edgeR::cpm(assay(sce), normalized.lib.sizes = FALSE, log = FALSE)

sce_pseudobulk <- aggregateData(
  sce,
  assay = "counts",
  by = c("Metacell"),
  fun = c("sum"),
  scale = FALSE # Should pseudo-bulks be scaled with the effective library size & multiplied by 1M?
)
assayNames(sce_pseudobulk) <- "counts"

colData(sce_pseudobulk) <- sample_metadata[cell==Metacell] %>% 
  .[,c("cell","Metacell","stage","sample","barcode","celltype.mapped","celltype.score","closest.cell")] %>%
  merge(cell2metacell[,.(ncells=.N),by="Metacell"],by="Metacell") %>%
  tibble::column_to_rownames("Metacell") %>% DataFrame

###############
## Normalise ##
###############

# sce_pseudobulk <- logNormCounts(sce_pseudobulk)

###################
## Sanity checks ##
###################

# cor(
#   colMeans(logcounts(sce_pseudobulk)),
#   metadata(sce_pseudobulk)$n_cells
# )

##########
## Save ##
##########

saveRDS(sce_pseudobulk, io$outfile)
