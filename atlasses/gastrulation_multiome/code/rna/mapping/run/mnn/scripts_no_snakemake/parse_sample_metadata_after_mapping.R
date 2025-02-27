###################
## Load settings ##
###################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
} else {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
}

#########
## I/O ##
#########

io$metadata <- paste0(io$basedir,"/results/rna/qc/sample_metadata_after_qc.txt.gz")
io$output.metadata <- paste0(io$basedir,"/results/rna/mapping/sample_metadata_after_mapping.txt.gz")
io$mapping.dir <- paste0(io$basedir,"/results/rna/mapping")

#############
## Options ##
#############

# opts$samples <- c("E8.5_rep1-E8.5_rep2")
opts$samples <- c(
	"E7.5_rep1",
	"E7.5_rep2", 
	"E8.0_rep1",
	"E8.0_rep2",
	"E8.5_rep1",
	"E8.5_rep2"
)

###############
## Load data ##
###############

# Load mapping results
mapping.dt <- opts$samples %>% map(function(x) 
  readRDS(sprintf("%s/mapping_mnn_%s.rds",io$mapping.dir,x))$mapping %>% .[,c("cell","celltype.mapped","celltype.score","closest.cell")] %>% as.data.table
) %>% rbindlist

###########
## Merge ##
###########

sample_metadata <- fread(io$metadata) %>% 
  # .[,c("celltype.mapped.y","celltype.score.y","closest.cell.y"):=NULL] %>%
  merge(mapping.dt,by="cell",all.x=TRUE)

head(sample_metadata)

#################
## Save output ##
#################

fwrite(sample_metadata, io$output.metadata, sep="\t", na="NA", quote=F)


