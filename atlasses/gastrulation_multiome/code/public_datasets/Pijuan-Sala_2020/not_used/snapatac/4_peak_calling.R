library(SnapATAC)
library(purrr)
library(data.table)


source(here::here("settings.R"))
source(here::here("public_datasets/Pijuan-Sala_2020/snapatac/snapatac_settings.R"))


opts$cores          <- 8
opts$out_prefix     <- "peaks"


print("loading rds file...")

snap <- readRDS(snapio$rds_file)
snap


snaptools <- system("which snaptools", intern = TRUE)
snaptools
macs2 <- system("which macs2", intern = TRUE)
macs2

print("running macs2...")

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

print("saving bed file...")

bed <- as.data.frame(peaks)[, 1:3]

fwrite(bed, 
       snapio$peaks, 
       sep = "\t", 
       na = "NA", 
       row.names = FALSE, 
       col.names = FALSE, 
       quote = FALSE)

