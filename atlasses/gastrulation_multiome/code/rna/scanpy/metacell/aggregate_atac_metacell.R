#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
  source("/Users/ricard/gastrulation_multiome_10x/utils.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
  source("/homes/ricard/gastrulation_multiome_10x/utils.R")
} else {
  stop("Computer not recognised")
}

# I/O
# io$cell2metacell <- "/Users/ricard/data/gastrulation_multiome_10x/results/rna/metacells/cell2metacell_1000metacells.txt.gz"
io$cell2metacell <- paste0(io$basedir,"/results/rna_atac/metacells/cell2metacell_2500metacells.txt.gz")
io$outfile <- paste0(io$basedir,"/results/rna_atac/metacells/atac_SingleCellExperiment_2500metacells.rds")

########################
## Load cell metadata ##
########################

sample_metadata <- fread(io$metadata) %>%
  .[pass_atacQC==TRUE & pass_rnaQC==TRUE & doublet_call==FALSE]

########################
## Load cell2metacell ##
########################

cell2metacell <- fread(io$cell2metacell) 

# cell2metacell[,.N,by="Metacell"] %>% View

########################
## Load ArchR project ##
########################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/atac/archR/load_archR_project.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/atac/archR/load_archR_project.R")
} else {
  stop("Computer not recognised")
}

ArchRProject.filt <- ArchRProject[sample_metadata$cell]


stopifnot(cell2metacell$cell %in% rownames(ArchRProject.filt))
# Add metacell information to the archR's colData

tmp <- cell2metacell %>% setkey(cell) %>% .[rownames(ArchRProject.filt)] %>%
  as.data.frame() %>% tibble::column_to_rownames("cell")
mean(is.na(tmp$Metacell))

ArchRProject.filt <- addCellColData(
  ArchRProject.filt,
  data = tmp[,1], 
  name = "Metacell",
  cells = rownames(tmp)
)

###################################
## Aggregate counts per metacell ##
###################################

atac_metacell.se <- getGroupSE(ArchRProject.filt, groupBy = "Metacell", useMatrix = "PeakMatrix", divideN = FALSE)

# Define peak names
peak_names <- rowData(atac_metacell.se) %>% as.data.table %>% .[,idx:=sprintf("%s:%s-%s",seqnames,start,end)] %>% .$id
rownames(atac_metacell.se) <- peak_names

# Add colData
colData(atac_metacell.se) <- sample_metadata[cell==Metacell] %>%
  .[,c("cell","Metacell","stage","sample","barcode","celltype.mapped","celltype.score","closest.cell")] %>%
  merge(cell2metacell[,.(ncells=.N),by="Metacell"],by="Metacell") %>%
  tibble::column_to_rownames("Metacell") %>% DataFrame


###############
## Normalise ##
###############

# SingleCellExperiment::SingleCellExperiment()

# clusts <- as.numeric(quickCluster(sce, method = "igraph", min.size = 100, BPPARAM = mcparam))
# # clusts <- as.numeric(quickCluster(sce))
# min.clust <- min(table(clusts))/2
# new_sizes <- c(floor(min.clust/3), floor(min.clust/2), floor(min.clust))
# sce <- computeSumFactors(sce, clusters = clusts, sizes = new_sizes, max.cluster.size = 3000)

# tmp <- colSums(assay(atac_metacell.se))

#   sce <- logNormCounts(sce)


##########
## Save ##
##########

saveRDS(atac_metacell.sce, io$outfile)
