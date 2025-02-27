library(Rsamtools)
library(GenomicAlignments)
library(rtracklayer)

source(here::here("settings.R"))



io$indir <- file.path(io$rawdata, "processed/rna/bams_by_celltype/merged")

io$outdir <- file.path(io$rawdata, "processed/rna/bigwigs")

dir.create(io$outdir, recursive = TRUE)

bams <- dir(io$indir, pattern = ".bam$", full = TRUE)
bams=bams[23:length(bams)]

walk(bams, ~{
  print(paste("converting", basename(.x), "to bigwig"))
  outfile <- paste0(io$outdir, "/", gsub(".bam", ".bw", basename(.x)))
  
  # bamfile <- system.file("extdata", .x, package="Rsamtools", mustWork=TRUE)
  # show(bamfile)
  
  alignment <- readGAlignments(.x)
  reads_coverage <- coverage(alignment)
  
  export.bw(reads_coverage, con = outfile)
})

