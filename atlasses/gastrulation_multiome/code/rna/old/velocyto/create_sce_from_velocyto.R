library(Matrix)
library(SingleCellExperiment)

#####################
## define settings ##
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
io$loom.dir <- "/hps/nobackup2/research/stegle/users/ricard/gastrulation_multiome_10x/processed/rna/velocyto/loom"
io$outfile <- "/hps/nobackup2/research/stegle/users/ricard/gastrulation_multiome_10x/processed/rna/velocyto/SingleCellExperiment_velocyto.rds"

# Options
opts$samples <- c(
  "E7.5_rep1",
  "E7.5_rep2",
  "E8.0_rep1",
  "E8.0_rep2",
  "E8.5_rep1",
  "E8.5_rep2"
)

opts$rename.samples <- c(
  "^multiome1:" = "E8.5_rep1#",
  "^multiome2:" = "E8.5_rep2#",
  "^E8_0_rep1_multiome:" = "E8.0_rep1#",
  "^E8_0_rep2_multiome:" = "E8.0_rep2#",
  "^rep1_L001_multiome:" = "E7.5_rep1#",
  "^rep2_L002_multiome:" = "E7.5_rep2#",
  "x$" = "-1"
)

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(io$metadata) %>%
  .[pass_rnaQC==TRUE & doublet_call==FALSE & sample%in%opts$samples]
cells <- sample_metadata$cell

#####################
## Load loom files ##
#####################

loom.mtx <- opts$samples %>% map(~ {
  tmp <- velocyto.R::read.loom.matrices(sprintf("%s/%s.loom",io$loom.dir,.))
  tmp[["ambiguous"]] <- NULL
  # doesnt work because cells have a different naming
  # subset.cells <- sample_metadata[sample==.,cell]
  # tmp[["spliced"]] <- tmp[["spliced"]][,subset.cells]
  # tmp[["unspliced"]] <- tmp[["unspliced"]][,subset.cells]
  return(tmp)
}); names(loom.mtx) <- opts$samples

###########################
## Concatenate and parse ##
###########################

concat_spliced.mtx <- do.call("cbind",opts$samples %>% map(~ loom.mtx[[.]][["spliced"]]))
concat_unspliced.mtx <- do.call("cbind",opts$samples %>% map(~ loom.mtx[[.]][["unspliced"]]))

# Rename cells
colnames(concat_spliced.mtx) <- str_replace_all(colnames(concat_spliced.mtx),opts$rename.samples) 
colnames(concat_unspliced.mtx) <- str_replace_all(colnames(concat_unspliced.mtx),opts$rename.samples) 
  
concat_spliced.mtx <- concat_spliced.mtx[,cells]
concat_unspliced.mtx <- concat_unspliced.mtx[,cells]

#################################
## create SingleCellExperiment ##
#################################

sce <- SingleCellExperiment(
  assays = list(spliced = concat_spliced.mtx, unspliced = concat_unspliced.mtx)
)

colData(sce) <- sample_metadata %>% as.data.frame %>% tibble::column_to_rownames("cell") %>%
  .[colnames(sce),] %>% DataFrame()

############################
## Calculate size factors ##
############################

# Save
saveRDS(sce, io$outfile)
