library(SnapATAC)
library(purrr)
library(data.table)


source(here::here("settings.R"))
source(here::here("atac/SnapATAC/snapatac_settings.R"))

opts$out_prefix     <- "5k_bins"
opts$use_cell_types <- TRUE # call peaks for each cell type (or for each cluster)
opts$cores          <- 8

snap <- readRDS(snapio$rds_file)
snap
print("number of cells per cluster....")
table(snap@cluster)

if (opts$use_cell_types) {
  snap@cluster <- as.factor(snap@metaData$cell_type)
  io$out_file <- gsub(".bed", "_celltype.bed", io$out_file)
}

# switch between laptop and cluster
#snap@file <- gsub("/bi/scratch/Stephen_Clark/gastrulation_multiome_10x", io$basedir, snap@file)

snaptools <- system("which snaptools", intern = TRUE)
stopifnot(length(snaptools) == 1)
snaptools
macs2 <- system("which macs2", intern = TRUE)
stopifnot(length(macs2) == 1)
macs2



peaks <- runMACSForAll(
  obj=snap, 
  output.prefix=opts$out_prefix,
  path.to.snaptools=snaptools,
  path.to.macs=macs2,
  gsize="mm", 
  min.cells = 100,
  buffer.size=500, 
  num.cores=opts$cores,
  macs.options="--nomodel --shift 37 --ext 73 --qval 1e-2 -B --SPMR --call-summits",
  tmp.folder=tempdir()
)

# save peaks as bed file

bed <- as.data.frame(peaks)[, 1:3]

fwrite(bed, 
       snapio$peaks, 
       sep = "\t", 
       na = "NA", 
       row.names = FALSE, 
       col.names = FALSE,
       quote = FALSE)

