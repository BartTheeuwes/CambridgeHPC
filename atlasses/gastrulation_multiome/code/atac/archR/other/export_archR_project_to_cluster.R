suppressPackageStartupMessages(library(ArchR))

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

file <- paste0(io$archR.directory,"/Save-ArchR-Project.rds")
ArchRProj <- readRDS(file)
ArchRProj@projectMetadata$outputDirectory <- io$archR.directory

ArchRProj@peakAnnotation$Motif_cisbp$Positions <- paste0(io$archR.directory,"/Annotations/Motif_cisbp-Positions-In-Peaks.rds")
ArchRProj@peakAnnotation$Motif_cisbp$Matches <- paste0(io$archR.directory,"/Annotations/Motif_cisbp-Matches-In-Peaks.rds")
# ArchRProj@peakAnnotation$Motif_homer$Positions <- paste0(io$archR.directory,"/Annotations/Motif_homer-Positions-In-Peaks.rds")
# ArchRProj@peakAnnotation$Motif_homer$Matches <- paste0(io$archR.directory,"/Annotations/Motif_homer-Matches-In-Peaks.rds")
ArchRProj@peakAnnotation$Motif_JASPAR2020_human$Positions <- paste0(io$archR.directory,"/Annotations/Motif_JASPAR2020_human-Positions-In-Peaks.rds")
ArchRProj@peakAnnotation$Motif_JASPAR2020_human$Matches <- paste0(io$archR.directory,"/Annotations/Motif_JASPAR2020_human-Matches-In-Peaks.rds")
# ArchRProj@peakAnnotation$Motif_JASPAR2020_mouse$Positions <- paste0(io$archR.directory,"/Annotations/Motif_JASPAR2020_mouse-Positions-In-Peaks.rds")
# ArchRProj@peakAnnotation$Motif_JASPAR2020_mouse$Matches <- paste0(io$archR.directory,"/Annotations/Motif_JASPAR2020_mouse-Matches-In-Peaks.rds")


# ArchRProj@sampleColData$ArrowFiles <- sprintf("%s/ArrowFiles/%s.arrow",io$archR.directory,opts$samples)
ArchRProj@sampleColData$ArrowFiles <- sprintf("%s/ArrowFiles/%s.arrow",io$archR.directory,rownames(ArchRProj@sampleColData))
rownames(ArchRProj@sampleColData) <- opts$samples

saveRDS(ArchRProj@peakAnnotation,sprintf("%s/Annotations/peakAnnotation.rds",io$archR.directory))
saveRDS(ArchRProj, file)
# saveArchRProject(ArchRProj, io$archR.directory)
