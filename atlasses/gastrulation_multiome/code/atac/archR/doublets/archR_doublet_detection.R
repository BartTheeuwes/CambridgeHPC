
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

names(ArchRProject@embeddings)

#####################
## Define settings ##
#####################

io$metadata <- paste0(io$basedir,"/processed/atac/archR/sample_metadata_after_archR.txt.gz")
io$outdir <- paste0(io$basedir,"/results/atac/archR/doublets")

opts$samples <- c(
  "E7.5_rep1",
  "E7.5_rep2",
  "E8.5_rep1",
  "E8.5_rep2"
)

########################
## Load cell metadata ##
########################

sample_metadata <- fread(io$metadata) %>%
  .[pass_atacQC==TRUE] %>%
  .[,archR_cell:=sprintf("%s#%s",sample,barcode)] %>%
  .[sample%in%opts$samples]

stopifnot(sample_metadata$archR_cell %in% rownames(ArchRProject))

##################
## Subset ArchR ##
##################

ArchRProject.filt <- ArchRProject[sample_metadata$archR_cell]
table(getCellColData(ArchRProject.filt,"Sample")[[1]])

##############################
## Calculate Doublet Scores ##
##############################

# E8.5_rep1 (1 of 1) : UMAP Projection R^2 = 0.9947 (should be >0.9)
ArchRProject.filt <- addDoubletScores(
    input = ArchRProject.filt,
    useMatrix = "TileMatrix",
    k = 10,             # Refers to how many cells near a "pseudo-doublet" to count.
    knnMethod = "UMAP", # Refers to the embedding to use for nearest neighbor search with doublet projection.
    LSIMethod = 1,
    UMAPParams = list(n_neighbors = 40, min_dist = 0.4, metric = "euclidean", verbose =FALSE),
)

# Adding doublet scores will create plots in the "QualityControl" directory. There are 3 plots associated with each of your samples in this folder:
# - Doublet Enrichments - These represent the enrichment of simulated doublets nearby each single cell compared to the expected if we assume a uniform distribution.
# - Doublet Scores - These represent the significance (-log10(binomial adjusted p-value)) of simulated doublets nearby each single cell compared to the expected if we assume a uniform distribution. We have found this value to be less consistent than the doublet enrichments and therefore use doublet enrichments for doublet identification.
# - Doublet Density - This represents the density of the simulated doublet projections. This allows you to visualize where the synthetic doublets were located after projection into your 2-dimensional embedding.

##########
## Save ##
##########

dt <- getCellColData(ArchRProject.filt, c("DoubletScore","DoubletEnrichment")) %>% 
  as.data.table(keep.rownames = T) %>% setnames("rn","cell") %>%
  .[,cell:=stringr::str_replace_all(cell,"#","_")]

fwrite(dt, sprintf("%s/doublet_scores.txt.gz",io$outdir), quote=F, na="NA", sep="\t")
