#########
## I/O ##
#########

# Load default settings
if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
}

io$mapping.dir <- paste0(io$basedir,"/results/rna/mapping")
io$mapping.dir.soupX <- paste0(io$basedir,"/results/rna_soupX/mapping")


opts$samples <- c(
  # "E7.5_rep1",
  # "E7.5_rep2",
  # "E8.0_rep1",
  # "E8.0_rep2",
  "E8.5_rep1",
  "E8.5_rep2"
)

####################################
## Load mapping results, no soupX ##
####################################

mapping.dt.standard <- opts$samples %>% map(function(x) 
  readRDS(sprintf("%s/mapping_mnn_%s.rds",io$mapping.dir,x))$mapping %>% .[,c("cell","celltype.mapped","celltype.score","closest.cell")] %>% as.data.table
) %>% rbindlist

######################################
## Load mapping results after soupX ##
######################################

mapping.dt.soupX <- opts$samples %>% map(function(x) 
  readRDS(sprintf("%s/mapping_mnn_%s.rds",io$mapping.dir.soupX,x))$mapping %>% .[,c("cell","celltype.mapped","celltype.score","closest.cell")] %>% as.data.table
) %>% rbindlist

###########
## Merge ##
###########

mapping.dt <- mapping.dt.standard %>% 
  merge(mapping.dt.soupX,by="cell",all.x=TRUE)

mean(mapping.dt$celltype.mapped.x == mapping.dt$celltype.mapped.y,na.rm=T)

mapping.dt[celltype.mapped.x != celltype.mapped.y,c("celltype.mapped.x","celltype.score.x","celltype.mapped.y","celltype.score.y")]
# mapping.dt[celltype.mapped.x==celltype.mapped.y]
# mapping.dt[celltype.mapped.x!=celltype.mapped.y] %>% View
