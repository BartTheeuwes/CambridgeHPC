library(SnapATAC)
library(purrr)
library(data.table)
library(viridisLite)
library(ggplot2)
library(GenomicRanges)

source(here::here("settings.R"))
source(here::here("public_datasets/Pijuan-Sala_2020/snapatac/snapatac_settings.R"))


######## needs to be run with R/3.5.2 ########


io$blacklist_url         <- "https://www.encodeproject.org/files/ENCFF547MET/@@download/ENCFF547MET.bed.gz"
# "http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/mm10-mouse/mm10.blacklist.bed.gz"
io$gene_anno             <- "/bi/scratch/Stephen_Clark/annotations/Mmusculus_genes_BioMart.87.txt"

### options ###

opts$fragments_cutoff    <- c(5e3, 1e6)
# opts$promoter_cutoff     <- c(0.15, 0.6)

print("loading rds file...")

snap <- readRDS(snapio$rds_file)


snap@metaData$UMI <- log10(snap@metaData$UM + 1)
p <- ggplot(snap@metaData, aes(UMI)) + geom_density(alpha = 0.5)
p

print("filtering cells...")

cells_keep <- which(snap@metaData$UM >= opts$fragments_cutoff[1] & snap@metaData$UM < opts$fragments_cutoff[2])

snap
snap <- snap[cells_keep]
snap

print("binarising...")

snap <- makeBinary(snap, mat = "bmat")


print("filtering loci...")

blacklist <- paste0(tempfile(), ".gz")
download.file(io$blacklist_url, blacklist)

blacklist_gr <- fread(blacklist) %>%
  makeGRangesFromDataFrame(seqnames.field = "V1",
                           start.field = "V2",
                           end.field = "V3")

idy <- queryHits(findOverlaps(snap@feature, blacklist_gr))
if(length(idy) > 0){snap <- snap[,-idy, mat="bmat"]}
snap

# Second, we remove unwanted chromosomes

chr.exclude <- seqlevels(snap@feature)[grep("random|chrM", seqlevels(snap@feature))]
idy <- grep(paste(chr.exclude, collapse="|"), snap@feature)
if(length(idy) > 0){snap <- snap[,-idy, mat="bmat"]}
snap


# Third, the bin coverage roughly obeys a log normal distribution. 
# We remove the top 5% bins that overlap with invariant features such as 
# promoters of the house keeping genes

bin.cov <- log10(Matrix::colSums(snap@bmat)+1)

hist(bin.cov[bin.cov > 0], 
     xlab="log10(bin cov)", 
     main="log10(Bin Cov)", 
     col="lightblue", 
     xlim=c(0, 5))

bin.cutoff <- quantile(bin.cov[bin.cov > 0], 0.95)
idy        <- which(bin.cov <= bin.cutoff & bin.cov > 0)
snap       <- snap[, idy, mat="bmat"];
snap

print("saving...")

saveRDS(snap, snapio$rds_file)

