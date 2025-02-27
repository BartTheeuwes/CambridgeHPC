
####################
## Load libraries ##
####################

suppressPackageStartupMessages(library(Signac))

#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/public_datasets/Pijuan-Sala_2020/settings.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/public_datasets/Pijuan-Sala_2020/settings.R")
} else {
  stop("Computer not recognised")
}

# Define I/O
io$outdir <- paste0(io$basedir,"/data/processed")
io$outfile <- paste0(io$basedir,"/data/processed/fragments.tsv.gz")

###############
## Load data ##
###############

m <- Matrix::readMM(io$matrix)
barcodes <- fread(io$barcodes, header=F)[[1]]
features <- fread(io$features, header=F)[[1]]

#################################
## Create fragments data.table ##
#################################

# alternative: SAVE .TXT FILE PER BARCODE AND THEN CONCATENATE

dt <- barcodes %>% head(n=10) %>% map(function(i) {
  rownames(m[as.logical(m[,i]==1),i]) %>% 
    stringr::str_split(.,pattern="-") %>%
    map(as.data.frame.list) %>%
    rbindlist(use.names=F) %>%
    setnames(c("chr","start","end")) %>%
    .[,barcode:=i]
}) %>% rbindlist

##########
## Save ##
##########

fwrite(dt, io$outfile, sep="\t", col.names = F)
