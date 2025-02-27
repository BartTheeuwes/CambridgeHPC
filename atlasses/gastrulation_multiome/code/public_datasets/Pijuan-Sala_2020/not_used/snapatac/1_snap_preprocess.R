library(SnapATAC)
library(purrr)
library(data.table)



source(here::here("settings.R"))
source(here::here("public_datasets/Pijuan-Sala_2020/snapatac/snapatac_settings.R"))


######## needs to be run with R/3.5.2 ########

# requires snaptools to be installed/loaded #



io$chrom_size_url <- "http://hgdownload.cse.ucsc.edu/goldenpath/mm10/bigZips/mm10.chrom.sizes"


opts$bin_sizes <- c(5000)
opts$cores <- 8



# first sort bam file by position

#cmd <- paste("samtools sort -o", io$sorted_bam, io$bam_file)

# convert fragments file to bed
#frags <- fread(cmd = paste("zcat <", snapio$frags_file))
frags <- fread(snapio$frags_file)
frags

# sort by barcode
setorder(frags, V4)

# re-format chromosome column
frags[, V1 := gsub("chr", "", V1) %>% paste0("chr", .)]
frags

# save as bed file without compression
bed_file <- tempfile(fileext = ".bed")
fwrite(frags, bed_file, sep = "\t", na = "NA", col.names = FALSE, row.names = FALSE, quote = FALSE)

# download chromosome sizes
chr_sizes <- tempfile()
download.file(io$chrom_size_url, chr_sizes)
fread(chr_sizes)

# run snaptools commands


cmd1 <- paste0(
  "snaptools snap-pre",
  " --input-file=", bed_file,
  " --output-snap=", snapio$snap_file,
  " --genome-name=mm10",
  " --genome-size=", chr_sizes,
  " --min-mapq=30",
  " --min-flen=50",
  " --max-flen=1000",
  " --keep-chrm=TRUE",
  " --keep-single=FALSE",
  " --keep-secondary=FALSE",
  " --overwrite=TRUE",
  " --max-num=200000",
  " --min-cov=500",
  " --verbose=TRUE"
)
cmd1
system(cmd1)

cmd2 <- paste0(
  "snaptools snap-add-bmat",
  " --snap-file ", snapio$snap_file,
  " --bin-size-list ", paste(opts$bin_sizes, collapse = " ")
)
cmd2
system(cmd2)


snap <- createSnap(snapio$snap_file, "Pijuan-Sala_2020", do.par = TRUE, num.cores = opts$cores)
snap
showBinSizes(snapio$snap_file)

head(snap@metaData)



snap <- addBmatToSnap(snap, bin.size=5000, num.cores=opts$cores)

saveRDS(snap, snapio$rds_file)
