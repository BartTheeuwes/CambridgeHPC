# TO-DO: TEST if PARALLEL PROCESSING OWKRS
library(furrr)

#################
## Description ##
#################

# script to merge multiple CellRanger fragment files
# barcodes are re-named to ensure they are unique between samples
# this is required since Signac has no way to handle multiple fragment files 
# with overlapping barcode names
# only cells which pass cellranger QC are kept
# script requires 'htslib' to be installed/loaded in order to perform compression and indexing

#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
} else {
  stop("Computer not recognised")
}

#########
## I/O ##
#########

io$fragment_files <- opts$samples %>%
  map_chr(~ sprintf("%s/original/%s/atac_fragments.tsv.gz",io$basedir,.)) %>%
  set_names(opts$samples)
# io$fragment_files <- c("E7.5_rep1" = sprintf("%s/original/E7.5_rep1/atac_fragments.tsv.gz",io$basedir))
  
io$cell_metrics_files <- opts$samples %>%
  map_chr(~ sprintf("%s/original/%s/per_barcode_metrics.csv",io$basedir,.)) %>%
  set_names(opts$samples)
# io$cell_metrics_files <- c("E7.5_rep1" = sprintf("%s/original/E7.5_rep1/per_barcode_metrics.csv",io$basedir))

io$fragments_outfile  <- file.path(io$signac.dir,"merged_fragments_test.tsv")
io$metrics_outfile    <- file.path(io$signac.dir,"cell_metrics_test.tsv.gz")

#############
## Options ##
#############

opts$cores <- 2
opts$chr <- paste0("chr",c(1:19,"X","Y"))

if (opts$cores > 1){
  plan(multisession, workers = opts$cores)
} else {
  plan(sequential)
}


###############
## Load data ##
###############

# Load and concatenate barcode QC metrics
metrics <- future_map2(io$cell_metrics_files, names(io$cell_metrics_files), ~{
  fread(.x) %>%
    .[is_cell == 1] %>% # filter to only include cells passing 10X QC (which is lenient)
    # .[, cell := paste0(.y, "_", gsub("-.*", "", barcode))] %>%
    .[, cell := sprintf("%s_%s",.y,barcode)] %>%
    .[, sample := .y]
}) %>% rbindlist

# Load and concatenate fragment files
fragments <- future_map2(io$fragment_files, names(io$fragment_files), ~{
  dt <- fread(.x, sep="\t", header=F, verbose=F, showProgress=F, 
              colClasses=c("V1"="factor", "V2"="integer", "V3"="integer", "V4"="factor" ,"V5"="integer")) %>% 
    setnames(c("chr", "start", "end", "barcode", "count")) %>%
    .[barcode%in%unique(metrics[sample==.y,barcode])] %>%
    .[chr%in%opts$chr] %>%
    .[, cell := sprintf("%s_%s",.y,barcode)] %>%
    .[, .(chr,start,end,cell,count)]
  print(paste("number of cells in fragment file:", dt[, length(unique(cell))]))
  return(dt)
}) %>% rbindlist %>% setorder(chr,start,end)

print(paste("number of cells in merged fragment file:", fragments[, length(unique(cell))]))

##########
## Save ##
##########

fwrite(metrics, io$metrics_outfile, col.names = TRUE, sep = "\t")
fwrite(fragments, io$fragments_outfile, col.names = FALSE, sep = "\t")

print("compressing...")
system(paste("bgzip -f -@", opts$cores, io$fragments_outfile))

print("indexing...")
system(sprintf("tabix --preset=bed -f %s.gz", io$fragments_outfile))
