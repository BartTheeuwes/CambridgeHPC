#########
## I/O ##
#########

source("/Users/ricard/gastrulation_multiome_10x/settings.R")

io$mapping.dir <- paste0(io$basedir,"/results/rna/mapping")


##############################
## Load data from mapping 1 ##
##############################

opts$batches <- c("E8.5_rep1","E8.5_rep2")

mapping.dt.1 <- opts$batches %>% map(function(x) 
  readRDS(sprintf("%s/mapping_mnn_%s.rds",io$mapping.dir,x))$mapping %>% .[,c("cell","celltype.mapped","celltype.score","closest.cell")] %>% as.data.table
) %>% rbindlist

##############################
## Load data from mapping 2 ##
##############################

opts$batches <- c("E8.5_rep1-E8.5_rep2")

mapping.dt.2 <- opts$batches %>% map(function(x) 
  readRDS(sprintf("%s/mapping_mnn_%s.rds",io$mapping.dir,x))$mapping %>% .[,c("cell","celltype.mapped","celltype.score","closest.cell")] %>% as.data.table
) %>% rbindlist

###########
## Merge ##
###########

mapping.dt <- mapping.dt.1 %>% 
  merge(mapping.dt.2,by="cell",all.x=TRUE)

mapping.dt[celltype.mapped.x==celltype.mapped.y]
mapping.dt[celltype.mapped.x!=celltype.mapped.y] %>% View

#################
## Save output ##
#################

foo <- sample_metadata %>% merge(mapping.dt,by="cell",all.x=TRUE)
# io$metadata <- "/Users/ricard/data/gastrulation_multiome_10x/multiome2/sample_metadata.csv"
fwrite(foo, io$metadata, sep="\t", na="NA", quote=F)


