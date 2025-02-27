library(data.table)
library(purrr)
library(Rsamtools)


# merge bam files using Rsamtools

source(here::here("settings.R"))


io$indir <- file.path(io$rawdata, "processed/rna/bams_by_celltype")
io$outdir <- file.path(io$rawdata, "processed/rna/bams_by_celltype/merged")

dir.create(io$outdir, recursive = TRUE)

bams <- dir(io$indir, pattern = ".bam$", recursive = TRUE, full = TRUE) %>% 
  set_names(basename(.))

celltypes <- unique(names(bams))

# .x=celltypes[1]

walk(celltypes, ~{
  print(paste("merging into", .x))
  
  files <- bams[names(bams)==.x]
  
  print(files)
  
  mergeBam(
    files = files,
    destination = file.path(io$outdir, .x),
    indexDestination = TRUE
  )
  print("done")
})

print("removing old files")
file.remove(bams)




