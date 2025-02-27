################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--distance',  type="integer",            default=1e5,      help='Maximum distance for a linkage between a peak and a gene')
p$add_argument('--outfile',       type="character",                help='Output file')
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

# load default setings
source(here::here("settings.R"))
source(here::here("utils.R"))

## START TEST
args$distance <- 5e4
args$outfile <- file.path(io$basedir,"results_new/rna_atac/gene_regulatory_networks/pseudobulk/TF2gene_after_virtual_chip.txt.gz")
## END TEST

###########################
## Load virtual ChIP-seq ##
###########################

io$virtual_chip.mtx <- file.path(io$basedir,"results_new/rna_atac/virtual_chipseq/CISBP/virtual_chip.mtx")
virtual_chip.mtx <- readRDS(io$virtual_chip.mtx)

#########################################################
## Load peak2gene linkages using only genomic distance ##
#########################################################

peak2gene.dt <- fread(io$archR.peak2gene.all) %>%
  # peak2gene.dt <- fread(io$archR.peak2gene.nearest) %>%
  .[dist<=args$distance] %>%
  .[,peak:=sprintf("chr%s:%s-%s",chr,peak.start,peak.end)]

# Sanity checks
stopifnot(length(intersect(rownames(virtual_chip.mtx),unique(peak2gene.dt$peak)))>1e5)

#########################################################
## Link TFs to target genes using the virtual ChIP-seq ##
#########################################################

opts$min.chip.threshold <- 0.25

tf2gene_chip.dt <- colnames(virtual_chip.mtx) %>% map(function(i) {
  print(i)
  # Select target peaks (note that we only take positive correlations into account)
  target_peaks_i <- names(which(virtual_chip.mtx[,i]>=opts$min.chip.threshold))
  
  if (length(target_peaks_i)>=1) {
    
    tmp <- data.table(
      tf = i,
      peak = target_peaks_i,
      chip_score = virtual_chip.mtx[target_peaks_i,i]
    ) %>% merge(peak2gene.dt[peak%in%target_peaks_i,c("peak","gene","dist")], by="peak")
    
    # print(sprintf("%s: %s target peaks & %s target genes",i,length(target_peaks_i),length(unique(tmp$gene))))
    return(tmp)
  }
}) %>% rbindlist


# Save
fwrite(tf2gene_chip.dt, args$outfile, sep="\t")
