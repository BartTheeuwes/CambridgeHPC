library(SnapATAC)
library(BSgenome.Mmusculus.UCSC.mm10)
library(purrr)
library(data.table)
library(chromVAR)
library(SummarizedExperiment)
library(motifmatchr)
library(JASPAR2016)
library(GenomicRanges)


source(here::here("settings.R"))
source(here::here("atac/SnapATAC/snapatac_settings.R"))



snap <- readRDS(snapio$rds_file)
snap

# fix weird formatting of chromosome names and remove unmapped contigs
# peaks <- fread(snapio$peaks) 
# peaks
peaks <- as.data.table(snap@peak)

peaks[, seqnames := gsub("b'|'", "", seqnames)] 

stopifnot(ncol(snap@pmat) == nrow(peaks))


peaks_gr <- makeGRangesFromDataFrame(peaks)

snap@peak <- peaks_gr
snap@peak

real_chromosomes <- peaks$seqnames %like% "chr"
snap@peak <- snap@peak[real_chromosomes,]
snap@pmat <- snap@pmat[, real_chromosomes]
# 
# #snap@peak$name[1:10]
# class(snap@peak)




snap@peak

# run ChromVar


snap@mmat <- runChromVAR(snap, 
                         input.mat = "pmat", 
                         genome = BSgenome.Mmusculus.UCSC.mm10,
                         species = "Mus musculus")

snap

saveRDS(snap, snapio$rds_file)
