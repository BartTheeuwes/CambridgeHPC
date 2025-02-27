
################
## Define I/O ##
################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
  io$script <- "/Users/ricard/gastrulation_multiome_10x/rna_atac/DORCs/DORCs.R"
} else {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
  io$script <- "/homes/ricard/gastrulation_multiome_10x/rna_atac/DORCs/DORCs.R"
  io$tmpdir <- paste0(io$basedir,"/results/rna_atac/DORCs/tmp"); dir.create(io$tmpdir)
}
io$outdir <- paste0(io$basedir,"/results/rna_atac/DORCs")

####################
## Define options ##
####################

# samples
opts$samples <- c("E7.5_rep1", "E7.5_rep2", "E8.5_rep1", "E8.5_rep2")

# Test mode (subsetting cells)?
opts$test_mode <- FALSE

# Number of highly variable genes
opts$gene_window <- c(1e4, 2.5e4, 5e4)

# Number of cores
opts$ncores <- 2

# LSF params
opts$memory <- 60000

#########
## Run ##
#########

for (i in opts$gene_window) {

  # Define LSF command
  if (grepl("ricard",Sys.info()['nodename'])) {
    lsf <- ""
  } else if (grepl("ebi",Sys.info()['nodename'])) {
    lsf <- sprintf("bsub -M %s -n %d -o %s/window%d_%dcores.txt", opts$memory, opts$ncores, io$tmpdir,i, opts$ncores)
  }
  cmd <- sprintf("%s Rscript %s --samples %s --gene_window %d --ncores %d --outdir %s", 
                 lsf, io$script, paste(opts$samples,collapse=" "), i, opts$ncores, io$outdir)

  if (isTRUE(opts$test_mode)) cmd <- paste0(cmd, " --test_mode")

  # Run
  print(cmd)
  system(cmd)
}
