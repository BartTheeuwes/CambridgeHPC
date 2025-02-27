library(data.table)
library(purrr)


# split bam files by celltype using sinto 

# requires sinto and samtools to be loaded

source(here::here("settings.R"))


io$script_dir <- here::here("rna/bam") # this is the path of this script. Required as sinto just spits out the files here

io$samples <- c(
  "multiome1", 
  "multiome2", 
  "rep1_L001_multiome", 
  "rep2_L002_multiome",
  "E8_0_rep1_multiome",
  "E8_0_rep2_multiome"
)

io$infiles     <- file.path(
  io$rawdata, 
  io$samples,
  "outs/possorted_genome_bam.bam"
)

# io$cellmetrics <- file.path(
#   io$rawdata, 
#   io$samples,
#   "outs/per_barcode_metrics.csv"
# )

stopifnot(all(file.exists(io$infiles)))



io$outdir <- file.path(io$rawdata, "processed/rna/bams_by_celltype")
dir.create(io$outdir, recursive = T)

filenames <- data.table(file = opts$rename.samples, sample = names(opts$rename.samples))

metadata <- fread(io$metadata) %>% 
  merge(filenames, by = "sample")


######################################
### TESTING!!!! ######################

io$samples = io$samples[6]


walk(seq_along(io$samples), ~{
  
  samp <- io$samples[.x]
  bamfile <- io$infiles[.x]
  print(paste("processing sample:", bamfile))
  
  # index bamfile
  bai <- paste0(bamfile, ".bai")
  if (!file.exists(bai)){
    print("Index for bam file not found. Generating now...")
    cmd <- paste("samtools index", bamfile)
    cmd
    system(cmd)
  }
  
  
  
  
  
  # extract celltype info with which to split on and save as tsv
  cells <- metadata[file == samp & pass_rnaQC == TRUE & !is.na(celltype.mapped), .(barcode, celltype.mapped)] 
  tmp <- tempfile(fileext = ".txt")
  fwrite(cells, tmp, sep = "\t", quote = FALSE, na = "NA", col.names = FALSE)
  
  # run sinto command
  
  cmd <- paste(
    "sinto filterbarcodes",
    "-b", bamfile,
    "-c", tmp,
   # "-o", outdir,
    "--nproc", 16
  )
  print(cmd)
  system(cmd)
  
  # files are outputed to the script directory so need to be moved
  
  
  files <- dir(io$script_dir, pattern = ".bam$", full = TRUE)
  
  outdir_samp <- paste0(io$outdir, "/", samp)
  dir.create(outdir_samp, recursive = TRUE)
  outfiles <- paste0(outdir_samp, "/", basename(files))
  
  print("moving files...")
  
  file.copy(files, outfiles, overwrite = TRUE)
  file.remove(files)
  
  print("done")
  
})


# stop()
# 
# # now merge the bams from each sample
# 
# 
# outdirs <- file.path(io$outdir, io$samples)
# bams <- dir(outdirs, pattern = ".bam$", full = TRUE)
# names(bams) <- basename(bams)
# 
# walk(unique(names(bams)), ~{
#   print(paste("merging", .x))
#   tomerge <- bams[names(bams) == .x]
#   outfile <- file.path(io$outdir, .x)
#   cmd <- paste(
#     "samtools merge",
#     outfile,
#     paste(tomerge, collapse = " ")
#   )
#   print(cmd)
#   system(cmd)
#   print("indexing...")
#   cmd <- paste(
#     "samtools index",
#     outfile
#   )
#   print("removing old files...")
#   file.remove(tomerge)
#   print("done")
# })



