#########
## I/O ##
#########

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/public_datasets/Pijuan-Sala_2020/settings.R")
  io$script <- "/Users/ricard/gastrulation_multiome_10x/public_datasets/Pijuan-Sala_2020/differential/archr_differential_accessibility_peaks.R"
} else if(grepl("ebi",Sys.info()['nodename'])){
  source("/homes/ricard/gastrulation_multiome_10x/public_datasets/Pijuan-Sala_2020/settings.R")
  io$script <- "/homes/ricard/gastrulation_multiome_10x/public_datasets/Pijuan-Sala_2020/differential/archr_differential_accessibility_peaks.R"
  io$tmpdir <- paste0(io$basedir,"/results/differential/archR/tmp"); dir.create(io$tmpdir, showWarnings = F)
} else {
  stop("Computer not recognised")
}
io$outdir <- paste0(io$basedir,"/results/differential/archR"); dir.create(io$outdir, showWarnings = F)

#############
## Options ##
#############

# Statistical test
opts$statistical.test <- "wilcoxon"

# Testing mode
# opts$test_mode <- FALSE

# Define matrix
opts$matrix <- "GeneScoreMatrix"

#########
## Run ##
#########

# opts$celltypes <- c("Erythroid","Gut")

# for (i in head(opts$celltypes,n=3)) {
for (i in opts$celltypes) {
  for (j in opts$celltypes) {
    if (i!=j) {
      outfile <- sprintf("%s/%s_%s_vs_%s.txt.gz", io$outdir,opts$matrix,i,j)
      
      # Define LSF command
      if (grepl("ricard",Sys.info()['nodename'])) {
        lsf <- ""
      } else if (grepl("ebi",Sys.info()['nodename'])) {
        lsf <- sprintf("bsub -M 3000 -n 1 -o %s/%s_%s_vs_%s.txt", io$tmpdir,opts$matrix,i,j)
      }
      cmd <- sprintf("%s Rscript %s --groupA %s --groupB %s --matrix %s --test %s --outfile %s", lsf, io$script, i, j, opts$matrix, opts$statistical.test, outfile)
      # if (isTRUE(opts$test_mode)) cmd <- paste0(cmd, " --test_mode")
      
      # Run
      print(cmd)
      system(cmd)
    }
  }
}

