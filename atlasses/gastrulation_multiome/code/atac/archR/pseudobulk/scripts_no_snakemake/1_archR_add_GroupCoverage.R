# https://www.ArchRProject.com/bookdown/how-does-archr-make-pseudo-bulk-replicates.html

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

#####################
## Define settings ##
#####################

# I/O
# io$metadata <- paste0(io$basedir,"/results/atac/archR/celltype_assignment/sample_metadata_after_archR.txt.gz")
# io$metadata <- paste0(io$basedir,"/sample_metadata.txt.gz")

# Options
opts$samples <- c(
  "E7.5_rep1",
  "E7.5_rep2"
  # "E8.0_rep1",
  # "E8.0_rep2",
  # "E8.5_rep1",
  # "E8.5_rep2"
)

opts$to.merge <- c(
  "Erythroid3" = "Erythroid",
  "Erythroid2" = "Erythroid",
  "Erythroid1" = "Erythroid",
  "Blood_progenitors_1" = "Blood_progenitors",
  "Blood_progenitors_2" = "Blood_progenitors"
  # "Intermediate_mesoderm" = "Mixed_mesoderm",
  # "Paraxial_mesoderm" = "Mixed_mesoderm",
  # "Nascent_mesoderm" = "Mixed_mesoderm",
  # "Pharyngeal_mesoderm" = "Mixed_mesoderm"
  # "Visceral_endoderm" = "ExE_endoderm"
)

opts$remove.small.lineages <- TRUE

########################
## Load cell metadata ##
########################

sample_metadata <- fread(io$metadata) %>%
  .[pass_atacQC==TRUE & !is.na(celltype.predicted) & sample%in%opts$samples] %>%
  .[,celltype.predicted:=stringr::str_replace_all(celltype.predicted,opts$to.merge)]

if (opts$remove.small.lineages) {
  opts$min.cells <- 50
  sample_metadata <- sample_metadata %>%
    .[,N:=.N,by=c("celltype.predicted")] %>% .[N>opts$min.cells] %>% .[,N:=NULL]
}

##################
## Subset ArchR ##
##################

# Subset
ArchRProject.filt <- ArchRProject[sample_metadata$cell]

# Update archR metadata
sample_metadata.to.archr <- sample_metadata %>% 
  .[cell%in%rownames(ArchRProject.filt)] %>% setkey(cell) %>% .[rownames(ArchRProject.filt)] %>%
  as.data.frame() %>% tibble::column_to_rownames("cell")

stopifnot(all(sample_metadata.to.archr$TSSEnrichment_atac == getCellColData(ArchRProject.filt,"TSSEnrichment")[[1]]))

ArchRProject.filt <- addCellColData(
  ArchRProject.filt,
  data = sample_metadata.to.archr[["celltype.predicted"]], 
  name = "celltype.predicted",
  cells = rownames(sample_metadata.to.archr),
  force = TRUE
)

# print cell numbers
table(getCellColData(ArchRProject.filt,"Sample")[[1]])
table(getCellColData(ArchRProject.filt,"celltype.predicted")[[1]])

#########################
## Add Group Coverages ##
#########################

# Check if group Coverages already exist
ArchRProject.filt@projectMetadata$GroupCoverages

# This function will merge cells within each designated cell group for the generation of pseudo-bulk replicates 
# and then merge these replicates into a single insertion coverage file.
# Output: creates files in archR/GroupCoverages/celltype: [X]._.Rep[Y].insertions.coverage.h5
ArchRProject.filt <- addGroupCoverages(ArchRProject.filt, groupBy = "celltype.predicted", force = TRUE)

##########
## Save ##
##########

io$archR.projectMetadata <- paste0(io$archR.directory,"/projectMetadata.rds")
saveRDS(ArchRProject.filt@projectMetadata, io$archR.projectMetadata)

# saveArchRProject(ArchRProject.filt)