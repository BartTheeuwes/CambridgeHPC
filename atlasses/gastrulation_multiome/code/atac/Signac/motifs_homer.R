library(purrr)
library(data.table)

# loads difacc sites then runs homer for motif enrichement and discovery for 
# the two sets of sites (i.e. hyper-accessible sites for two groups of celltypes)


source(here::here("settings.R"))

io$difacc_file <- file.path(io$rawdata, "/processed/atac/signac/diffacc/Erythroid1_Erythroid2_Erythroid3_vs_Forebrain_Midbrain_Hindbrain_difacc.tsv.gz")
io$homer_out   <- file.path(io$rawdata, "/processed/atac/homer")
io$genome_fa   <- "/bi/scratch/Stephen_Clark/annotations/genomes/Mus_musculus.GRCm38.dna.primary_assembly.fa"
  
  
dir.create(io$homer_out, recursive = TRUE)

opts$cores     <- 8
opts$motif_length <- 6
opts$scan_size <- 200
opts$optimize_count <- 2
opts$cache <- 1e3 # cache in MB
opts$fdr_num <- 2 #5


opts$group2_as_background <- FALSE # only does enrichment on group 1 with group2 as background (instead of running each group separately with auto background)



difacc <- fread(io$difacc) 



(homer <- system("which findMotifsGenome.pl", intern = TRUE))


groups <- gsub("_difacc.tsv.gz", "", basename(io$difacc_file)) %>% 
  strsplit("_vs_") %>% 
  unlist()


peaks <- list(difacc[avg_logFC > 0 & p_val_adj < 0.05, .(locus)],
              difacc[avg_logFC < 0 & p_val_adj < 0.05, .(locus)]) %>% 
  set_names(groups) %>% 
  map(tidyr::separate, "locus", c("chr", "start", "end")) %>% 
  map(~.[, chr := gsub("chr", "", chr)])


if (opts$group2_as_background){
  peakfiles <- list(tempfile(fileext = ".txt"),
                    tempfile(fileext = ".txt"))
  
  
  walk2(peaks, peakfiles, fwrite, sep = "\t", na = "NA", col.names = FALSE, row.names = FALSE, quote = FALSE)
  outdir <- paste0(io$homer_out, "/", groups[1], "_vs_", groups[2])
  system(paste("rm -rd", outdir))
  
  cmd <- paste(
    #homer,
    "findMotifsGenome.pl",
    peakfiles[1], 
    io$genome_fa,#"mm10", 
    outdir,
    "-len", paste(opts$motif_length, collapse = ","),
    "-size", opts$scan_size,
    "-S", opts$optimize_count,
    "-p", opts$cores,
    "-cache", opts$cache,
    "-fdr", opts$fdr_num,
    "-bg", peakfiles[2]
    # "-preparsedDir", tmpdir
  )
  cmd
  system(cmd)
  
} else {
  walk(groups, ~{
    p <- peaks[[.x]]
    
    setorder(p, "chr", "start", "end")
    
    peakfile <- tempfile(fileext = ".txt")
    
    fwrite(p, peakfile, sep = "\t", na = "NA", col.names = FALSE, row.names = FALSE, quote = FALSE)
    fread(peakfile)
    #system(paste("head", io$genome_fa))
    outdir <- file.path(io$homer_out, .x)
    system(paste("rm -rd", outdir))
    
    tmpdir <- tempdir()
    
    cmd <- paste(
      #homer,
      "findMotifsGenome.pl",
      peakfile, 
      io$genome_fa,#"mm10", 
      outdir,
      "-len", paste(opts$motif_length, collapse = ","),
      "-size", opts$scan_size,
      "-S", opts$optimize_count,
      "-p", opts$cores,
      "-cache", opts$cache,
      "-fdr", opts$fdr_num
      # "-preparsedDir", tmpdir
    )
    cmd
    system(cmd)
  })
}






