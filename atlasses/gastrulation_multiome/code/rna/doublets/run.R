
################
## Define I/O ##
################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
  io$script <- "/Users/ricard/gastrulation_multiome_10x/rna/doublets/doublet_detection.R"
} else {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
  io$script <- "/homes/ricard/gastrulation_multiome_10x/rna/doublets/doublet_detection.R"
  io$tmpdir <- paste0(io$basedir,"/results/rna/doublets/tmp"); dir.create(io$tmpdir)
}
io$outdir <- paste0(io$basedir,"/results/rna/doublets")

####################
## Define options ##
####################

# opts$samples <- c(
#     "E7.5_rep1",
#     "E7.5_rep2",
#     "E8.5_rep1",
#     "E8.5_rep2"
# )

# Test mode (subsetting cells)?
opts$test_mode <- FALSE


# Number of highly variable genes
opts$hybrid_score_threshold <- 1.0

opts$memory <- 13000

#########
## Run ##
#########

for (i in opts$samples) {
  
  # Define LSF command
  if (grepl("ricard",Sys.info()['nodename'])) {
    lsf <- ""
  } else if (grepl("ebi",Sys.info()['nodename'])) {
    lsf <- sprintf("bsub -M %d -n 1 -o %s/%s_%s.txt", opts$memory, io$tmpdir,i,opts$hybrid_score_threshold)
  }
  cmd <- sprintf("%s Rscript %s --samples %s --hybrid_score_threshold %s --outdir %s", lsf, io$script, i, opts$hybrid_score_threshold, io$outdir)
  
  if (isTRUE(opts$test_mode)) cmd <- paste0(cmd, " --test")
  
  # Run
  print(cmd)
  system(cmd)
}
