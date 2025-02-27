library(purrr)
library(data.table)
library(furrr)

options(future.globals.maxSize = 2 * 1024^3) # 2 GiB global memory

source(here::here("settings.R"))
source(here::here("atac/SnapATAC/snapatac_settings.R"))

io$chrom_sizes <- "http://hgdownload.cse.ucsc.edu/goldenpath/mm10/bigZips/mm10.chrom.sizes"

opts$bin_sizes <- c(5000)
opts$n_cells   <- 4e3 # number of cells to subset to or FALSE

opts$min_cov   <- 1e3
opts$cores     <- 8

if (opts$cores > 1){
  plan(multisession, workers = opts$cores)
} else {
  plan(sequential)
}


# check input files
stopifnot(length(snapio$frags_files) == length(snapio$snap_files))
# check snaptools is loaded
snaptools <- system("which snaptools", intern = TRUE)
stopifnot(length(snaptools) == 1)

# read in cellRanger metrics
barcodes <- map2(snapio$metrics_files, snapio$samples, ~fread(.x)[, sample := .y]) %>%
  rbindlist()

# 10X cellRanger fragment files need to be sorted by barcode (instead of position)
# script requires lots of memory - say 100G

bedfiles <- gsub(".tsv.gz", ".bed", snapio$frags_files)
.x=1
future_map(seq_along(snapio$samples), ~{
  samp <- snapio$samples[[.x]]
  print(paste("processing", samp))
  
  # check file format
  # system(paste("zcat", snapio$frags_files[[.x]], "| head"))
  
  
  # read fragment file
  #dt <- fread(cmd = paste("zcat", .snapio$frags_files[[.x]])) # faster to use zcat than data.table
  dt <- fread(snapio$frags_files[[.x]])
  dt
  
  # filter out non-cells (to hopefully save on memory and processing time)
  cells_keep <- barcodes[sample == samp & is_cell == 1, barcode]
  
  if (opts$n_cells != FALSE && opts$n_cells < length(cells_keep)){
    cells_keep <- cells_keep[1:opts$n_cells]
  }
  
  setkey(dt, "V4")
  dt <- dt[V4 %in% cells_keep]
  # order by barcode
  setorder(dt, "V4")
  dt
  # save
  fwrite(dt, bedfiles[[.x]], sep = "\t", na = "NA", col.names = FALSE, quote = FALSE)
  paste(samp, "done")
})






chrom_sizes <- tempfile(fileext = ".txt")
download.file(io$chrom_sizes, chrom_sizes)

chrom_sizes_dt <- fread(chrom_sizes)
# remove shitty chromosomes
chrom_sizes_dt <- chrom_sizes_dt[V1 %in% paste0("chr", c(1:19, "X", "Y"))]
fwrite(chrom_sizes_dt, chrom_sizes, sep = "\t", na = "NA", col.names = FALSE, quote = FALSE)

# iterate over input/output files to run snaptools commands

future_map2(bedfiles, snapio$snap_files, ~{
  print("processing file:")
  print(.x)
  cmd1 <-  paste0("snaptools snap-pre",
                  " --input-file=", .x,
                  " --output-snap=", .y,
                  " --genome-name=mm10",
                  " --genome-size=", chrom_sizes,
                  " --min-mapq=30",
                  " --min-flen=50",
                  " --max-flen=1000",
                  " --keep-chrm=TRUE",
                  " --keep-single=FALSE",
                  " --keep-secondary=False",
                  " --overwrite=TRUE",
                  " --max-num=100000",
                  " --min-cov=", opts$min_cov,
                  " --verbose=TRUE")
  print(cmd1)
  system(cmd1)
  
  cmd2 <- paste0("snaptools snap-add-bmat",
                 " --snap-file=", .y,
                 " --bin-size-list ", paste(opts$bin_sizes, collapse = " "))
  
  print("adding bin matrix...")
  print(cmd2)
  system(cmd2)
  
  paste(.x, .y, "done")
  
})

walk(bedfiles, file.remove)

