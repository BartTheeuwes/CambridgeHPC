library(data.table)
library(purrr)

# requires ucsc_tools to be loaded

io <- list()
io$chrom_sizes <- "http://hgdownload.cse.ucsc.edu/goldenpath/mm10/bigZips/mm10.chrom.sizes"
io$indir       <- "/bi/scratch/Stephen_Clark/gastrulation_data/bedGraph"
io$outdir      <- "/bi/scratch/Stephen_Clark/gastrulation_data/bigwig"

dir.create(io$outdir, recursive = TRUE)

chrom_sizes <- paste0(io$outdir, "/chromsizes.txt")#tempfile(fileext = ".txt")
download.file(io$chrom_sizes, chrom_sizes)
fread(chrom_sizes)

infiles <- dir(io$indir, pattern = ".bedGraph")

walk(infiles, ~{
  
  outfile <- paste0(io$outdir, "/", gsub(".bedGraph", ".bw", basename(.x)))
  
  cmd <- paste(
    "bedGraphToBigWig",
    .x,
    chrom_sizes,
    outfile
  )
  
  cmd
  
  system(cmd)
  
})


