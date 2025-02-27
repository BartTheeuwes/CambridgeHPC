library(SnapATAC)
library(purrr)
library(data.table)
library(viridisLite)
library(ggplot2)
library(cowplot)
library(GenomicRanges)

source(here::here("settings.R"))
source(here::here("atac/SnapATAC/snapatac_settings.R"))

######## needs to be run with R/3.5.2 ########

# script adapted from https://github.com/r3fang/SnapATAC/blob/master/examples/10X_brain_5k/README.md

### input/output ###



io$blacklist_url         <- "https://www.encodeproject.org/files/ENCFF547MET/@@download/ENCFF547MET.bed.gz"
                          # "http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/mm10-mouse/mm10.blacklist.bed.gz"
io$gene_anno             <- "/bi/scratch/Stephen_Clark/annotations/Mmusculus_genes_BioMart.87.txt"


### options ###

opts$fragments_cutoff    <- c(1e3, 1e6)
opts$promoter_cutoff     <- c(0.15, 0.8)

opts$bin_size            <- 5e3
opts$cores               <- 8






stopifnot(length(snapio$snap_files) == length(snapio$samples))
stopifnot(length(snapio$snap_files) == length(snapio$metrics_files))


# import metadata

sample_metadata <- fread(io$metadata)

#################################
### Step 1. Barcode selection ###
#################################

# We select high-quality barcodes based on two criteria: 
# 1) number of unique fragments; 2) fragments in promoter ratio

barcodes <- map2(snapio$metrics_files, snapio$samples, ~fread(.x)[, sample := .y]) %>%
  rbindlist() %>%
  .[, promoter_ratio := (atac_TSS_fragments + 1) / (atac_fragments + 1)] %>%
  .[, UMI := log(atac_fragments + 1, 10)] %>% 
  .[, cell := paste0(sample, "_", gsub("-.*", "", barcode))] 

barcodes



# merge with metadata to get celltype
barcodes <- merge(barcodes, 
                  sample_metadata[, .(cell, cell_type = celltype.mapped, pass_rnaQC)], 
                  by = "cell", 
                  all = TRUE)
barcodes[, passQC := FALSE]

barcodes[is_cell == 1 &
           atac_fragments >= opts$fragments_cutoff[1] &
           atac_fragments <= opts$fragments_cutoff[2] &
           promoter_ratio >= opts$promoter_cutoff[1] &
           promoter_ratio <= opts$promoter_cutoff[2],
         passQC := TRUE
         ]
ggplot(barcodes[is_cell == 1], aes(sample, fill = passQC)) +
  geom_bar() +
  theme_cowplot()

p1 <- ggplot(barcodes[], aes(x= UMI, y= promoter_ratio)) + 
  geom_point(size=0.5, col="grey", alpha = 0.5) +
  theme_classic() +
  ggtitle("all barcodes") +
  ylim(0, 1) + xlim(0, 6) +
  labs(x = "log10(UMI)", y="promoter ratio") +
  geom_vline(xintercept = log10(opts$fragments_cutoff +1), linetype = "dashed", size = 0.1) +
  geom_hline(yintercept = opts$promoter_cutoff, linetype = "dashed", size = 0.1)


cowplot::save_plot(paste0(snapio$plots_out, "/qc_all_cells.pdf"), p1)


p2 <- ggplot(barcodes[is_cell == 1], aes(x= UMI, y= promoter_ratio)) + 
  geom_point(size=0.5, col="grey", alpha = 0.5) +
  theme_classic() +
  ggtitle("pass 10X QC") +
  ylim(0, 1) + xlim(0, 6) +
  labs(x = "log10(UMI)", y="promoter ratio")  +
  geom_vline(xintercept = log10(opts$fragments_cutoff + 1), linetype = "dashed", size = 0.1) +
  geom_hline(yintercept = opts$promoter_cutoff, linetype = "dashed", size = 0.1)

cowplot::save_plot(paste0(snapio$plots_out, "/qc_10X_passQC_cells.pdf"), p2)





# now load snap files

file.exists(snapio$snap_files)
snap <- createSnap(snapio$snap_files, snapio$samples)
snap

# correct any weird formatting of barcodes and add cell name to metadata

snap@barcode <- gsub('"', '', snap@barcode)

head(snap@metaData)

snap@metaData$barcode <- snap@barcode

snap@metaData$cell <- paste0(snap@sample,
                             "_",
                             gsub("-.*", "", snap@barcode) %>% gsub('"', '', .))

head(snap@metaData)


# subset snap

# (ncells <- nrow(snap@metaData))
# sub <- sample(ncells, 1000) %>% .[order(.)]
# snap_sub <- snap[sub]
# saveRDS(snap_sub, gsub(".rds", "_sub.rds", snapio$rds_file))
# snap=readRDS(gsub(".rds", "_sub.rds", snapio$rds_file))
# snap









# # filter


cells      <- barcodes[
                       cell %in% snap@metaData$cell &
                       is_cell == 1 &
                       atac_fragments >= opts$fragments_cutoff[1] &
                       atac_fragments <= opts$fragments_cutoff[2] &
                       promoter_ratio >= opts$promoter_cutoff[1] &
                       promoter_ratio <= opts$promoter_cutoff[2]
                       ]

cells[, .N, sample]

cells_keep <- which(snap@metaData$cell %in% cells$cell)
length(cells_keep)


snap <- snap[cells_keep]
snap



# update metadata

cells <- setDF(cells) %>%
  tibble::column_to_rownames("cell") %>%
  .[snap@metaData$cell,]

stopifnot(snap@barcode == cells$barcode)


snap@metaData <- cbind(snap@metaData, cells)

# ######################################
# ### Step 2. Add cell-by-bin matrix ###
# ######################################

# Next, we add the cell-by-bin matrix of 5kb resolution to the snap object.
# This function will automatically read the cell-by-bin matrix and add it to
# bmat slot of snap object.


map(snapio$snap_files, showBinSizes)
print("adding bin matrix...")
snap <- addBmatToSnap(snap,
                      bin.size  = opts$bin_size,
                      num.cores = opts$cores)


# ###################################
# ### Step 3. Matrix binarization ###
# ###################################

# We will convert the cell-by-bin count matrix to a binary matrix.
# Some items in the count matrix have abnormally high coverage perhaps due to
# the alignment errors. Therefore, we next remove 0.1% items of the highest
# coverage in the count matrix and then convert the remaining non-zero items to 1.
print("binarising...")
snap <- makeBinary(snap, mat = "bmat")

# #############################
# ### Step 4. Bin filtering ###
# #############################
#
# # First, we filter out any bins overlapping with the ENCODE blacklist to prevent
# # from potential artifacts.


blacklist <- tempfile(fileext =  ".gz")
download.file(io$blacklist_url, blacklist)

blacklist_gr <- fread(blacklist) %>%
  makeGRangesFromDataFrame(seqnames.field = "V1",
                           start.field = "V2",
                           end.field = "V3")
print("removing bins overlapping blacklist sites...")
idy <- queryHits(findOverlaps(snap@feature, blacklist_gr))
idy[1:10]
if(length(idy) > 0){snap <- snap[,-idy, mat="bmat"]}
snap
#
# # Second, we remove unwanted chromosomes
#

chr.exclude <- seqlevels(snap@feature)[grep("random|chrM", seqlevels(snap@feature))]

if (length(chr.exclude) > 0) {
  print("removing bins on unwanted chromosomes...")
  idy <- grep(paste(chr.exclude, collapse="|"), snap@feature)
  idy[1:10]
  if(length(idy) > 0){snap <- snap[,-idy, mat="bmat"]}
} else {
  print("all bins are within canonical chromosomes. nothing to filter.")
}

snap
#
#
# # Third, the bin coverage roughly obeys a log normal distribution.
# # We remove the top 5% bins that overlap with invariant features such as
# # promoters of the house keeping genes
#

print("removing bins with top/bottom 5% accessibility..")
bin.cov <- log10(Matrix::colSums(snap@bmat)+1)
#
# hist(bin.cov[bin.cov > 0],
#      xlab="log10(bin cov)",
#      main="log10(Bin Cov)",
#      col="lightblue",
#      xlim=c(0, 5))
#
bin.cutoff <- quantile(bin.cov[bin.cov > 0], 0.95)
idy        <- which(bin.cov <= bin.cutoff & bin.cov > 0)
idy[1:10]
snap       <- snap[, idy, mat="bmat"];
snap
#
# # # add gene annotation
# # genes <- fread(io$gene_anno) %>%
# #   setnames("symbol", "name") %>%
# #   makeGRangesFromDataFrame(keep.extra.columns = TRUE)
# # genes
# # snap <- createGmatFromMat(snap, genes = genes, num.cores = opts$cores)
# # snap
#snap <- makeBinary(snap, mat = "gmat")

# # save snap file
saveRDS(snap, snapio$rds_file)

# save metrics file
fwrite(barcodes, file.path(dirname(snapio$rds_file), "barcode_metrics.tsv.gz"), sep = "\t", na = "NA")
# 
# update metadata file
# sample_metadata <- fread(io$metadata)
# meta_new <- merge(sample_metadata,
#                   barcodes[, .(cell = paste0(sample, "_", gsub("-.*", "", barcode)), promoter_ratio, atac_fragments)],
#                   by = "cell",
#                   all.x = TRUE) %>%
#   .[, pass_snapQC := FALSE] %>%
#   .[cell %in% snap@metaData$cell, pass_snapQC := TRUE]
# 
# fwrite(meta_new, io$metadata)
