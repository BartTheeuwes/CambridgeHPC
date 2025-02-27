suppressPackageStartupMessages(library(SoupX))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(furrr))
suppressPackageStartupMessages(library(argparse))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--inputdir',       type="character",                    help='Output directory')
p$add_argument('--outputdir',       type="character",                    help='Output directory')
p$add_argument('--samples',         type="character",       nargs="+",   help='Samples')
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

# Load default settings
if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/10x_TET_TKO_EBs/settings.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/10x_TET_TKO_EBs/settings.R")
}

## START TEST ##
# args <- list()
# args$inputdir <- paste0(io$basedir,"/original")
# args$outputdir <- paste0(io$basedir,"/processed/test/SoupX")
# args$samples <- opts$batches[1]
## END TEST ##

# Define cores
# opts$cores <- snakemake@threads[[1]]
# if (opts$cores > 1){
#   plan(multisession, workers = opts$cores)
# }

# Sanity checks
stopifnot(args$samples%in%opts$batches)

###############
## Run SoupX ##
###############

future_map(args$samples, ~{
  
  # Define I/O
  io$sampledir <- sprintf("%s/%s_results",args$inputdir,.x)
  io$features_file <- sprintf("%s/%s_features.tsv.gz",io$sampledir,.x,.x) 
  io$clusters_file <- sprintf("%s/clustering/graphclust/clusters.csv",io$sampledir,.x,.x) 
  io$mat.outfile <- sprintf("%s/%s_matrix.mtx.gz",args$outputdir,.x)
  io$soup.outfile <- sprintf("%s/%s_soup.tsv.gz",args$outputdir,.x)
  
  print(.x)
  print(io$sampledir)
  print(io$features_file)
  print(io$clusters_file)
  print(io$mat.outfile)
  print(io$soup.outfile)
  
  # Create output directories
  dir.create(dirname(io$mat.outfile), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(io$soup.outfile), recursive = TRUE, showWarnings = FALSE)
  
  # load gene metadata from cell ranger output folder
  gene_metadata <- fread(io$features_file, select=c(1,2)) %>% 
    setnames(c("ens_id", "gene"))
  
  # read in clustering data from cell ranger output
  clusters_dt <- fread(io$clusters_file)
  clusters <- clusters_dt[, Cluster]
  names(clusters) <- clusters_dt[, Barcode]
  
  # Load CellRanger output
  tod <- Read10X(file.path(io$sampledir, "unfiltered"), prefix=sprintf("%s_unfiltered_",.x))
  toc <- Read10X(file.path(io$sampledir, "filtered"), prefix=sprintf("%s_",.x))
  
  # Sanity check
  stopifnot(ncol(tod) > ncol(toc))
  stopifnot(colnames(toc) == names(clusters))
  
  # Run SoupX 
  sc <- SoupChannel(tod, toc) %>%
    setClusters(clusters) %>%
    autoEstCont(forceAccept = TRUE)

  # Save adjusted matrix
  writeMM(adjustCounts(sc), io$mat.outfile)
  system(sprintf("pigz %s",io$mat.outfile))
  
  # Save soup counts
  soup_counts.dt <- setDT(sc$soupProfile, keep.rownames = "gene") %>%
    merge(gene_metadata, by = "gene", all.x = TRUE) %>% 
    .[order(-rank(counts))] 
  fwrite(soup_counts.dt, io$soup.outfile, sep = "\t", na = "NA", quote = FALSE)

})


