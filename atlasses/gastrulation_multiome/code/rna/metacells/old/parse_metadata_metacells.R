here::i_am("rna/scanpy/metacell/run/parse_metadata_metacells.R")

source(here::here("settings.R"))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
# p$add_argument('--metadata',    type="character",  help='Metadata file to use as input')
p$add_argument('--cell2metacell',    type="character",  nargs="+", help='Metacell results')
p$add_argument('--outfile',          type="character",               help='Output file')
args <- p$parse_args(commandArgs(TRUE))

###################
## Load settings ##
###################

## START TEST ##
args$metadata <- file.path(io$basedir,"results/rna/mapping/sample_metadata_after_mapping.txt.gz")
args$cell2metacell <- file.path(io$basedir,sprintf("results/rna/metacells/%s/cell2metacell_assignment.txt.gz",opts$samples))
args$outfile <- file.path(io$basedir,"results/rna/metacells/metacells_metadata.txt.gz")
## END TEST ##

stopifnot(file.exists(args$cell2metacell))

###########################
## Load metacell results ##
###########################

cell2metacell.dt <- args$cell2metacell %>% map(~ fread(.)) %>% rbindlist
# stopifnot(mapping_mnn.dt$cell%in%sample_metadata$cell)

table(cell2metacell.dt$sample)

###################
## Load metadata ##
###################

sample_metadata <- fread(args$metadata) %>%
  .[cell%in%cell2metacell.dt$metacell] %>%
  .[,c("cell","sample","stage","genotype","celltype.mapped")] %>%
  setnames("celltype.mapped","celltype")

###########
## Merge ##
###########

# to.save <- sample_metadata %>% 
#   merge(mapping_mnn.dt, by=c("cell","sample","class"))

#################
## Save output ##
#################

fwrite(cell2metacell.dt, args$outfile, sep="\t", na="NA", quote=F)
