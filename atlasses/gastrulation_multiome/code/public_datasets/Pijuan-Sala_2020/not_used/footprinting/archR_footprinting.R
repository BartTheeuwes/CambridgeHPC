
########################
## Load ArchR project ##
########################

source("/Users/ricard/gastrulation_multiome_10x/public_datasets/Pijuan-Sala_2020/archR/load_archR_project.R")

#####################
## Define settings ##
#####################

io$outdir <- paste0(io$basedir,"/results/chromVAR/archR")

#####################
## Update metadata ##
#####################
