
################
## Define I/O ##
################

# Load default settings
source(here::here("settings.R"))

# I/O
io$script <- here::here("rna/mapping/run/mnn/mapping_mnn.R")
io$outdir <- paste0(io$basedir,"/results/rna/mapping")

####################
## Define options ##
####################

opts$atlas_stages <- c(
  # "E6.5",
  # "E6.75",
  # "E7.0",
  # "E7.25",
  "E7.5",
  "E7.75",
  "E8.0",
  "E8.25",
  "E8.5"
  # "mixed_gastrulation"
)

opts$query_samples <- opts$samples

# Test mode (subsetting cells)?
opts$test_mode <- TRUE

if (opts$test_mode) {
  opts$memory <- 10000
} else {
  opts$memory <- 30000
}

opts$npcs <- 50
opts$n_neighbours <- 25

#########
## Run ##
#########

for (i in opts$query_samples) {
  
  # Define LSF command
  if (grepl("ricard",Sys.info()['nodename'])) {
    lsf <- ""
  } else if (grepl("ebi",Sys.info()['nodename'])) {
    lsf <- sprintf("bsub -M %d -n 1 -o %s/%s.txt", opts$memory, io$tmpdir,paste(i,collapse=" "))
  }
  cmd <- sprintf("%s Rscript %s --atlas_stages %s --query_samples %s --npcs %d --n_neighbours %d --outdir %s", 
                 lsf, io$script, paste(opts$atlas_stages,collapse=" "), paste(i,collapse=" "), opts$npcs, opts$n_neighbours, io$outdir)
  if (isTRUE(opts$test_mode)) cmd <- paste0(cmd, " --test")

  # Run
  print(cmd)
  system(cmd)
}
