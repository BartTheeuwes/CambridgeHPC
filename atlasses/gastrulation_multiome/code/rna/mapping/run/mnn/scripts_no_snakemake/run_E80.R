
################
## Define I/O ##
################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
  io$script <- "/Users/ricard/gastrulation_multiome_10x/rna/mapping/run/mnn/mapping_mnn.R"
} else {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
  io$script <- "/homes/ricard/gastrulation_multiome_10x/rna/mapping/run/mnn/mapping_mnn.R"
  io$tmpdir <- paste0(io$basedir,"/results/rna/mapping/tmp"); dir.create(io$tmpdir)
}
io$outdir <- paste0(io$basedir,"/results/rna_soupX/mapping")

####################
## Define options ##
####################

opts$atlas_stages <- c(
  # "E6.5",
  # "E6.75",
  # "E7.0",
  # "E7.25",
  # "E7.5",
  "E7.75",
  "E8.0",
  "E8.25",
  "E8.5"
  # "mixed_gastrulation"
)

opts$query_samples <- c(
  "E8.0_rep1",
  "E8.0_rep2"
)

# opts$samples <- "SIGAA6_E85_2_Dnmt3aKO_Dnmt3b_WT_L001"

# Test mode (subsetting cells)?
opts$test_mode <- FALSE

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
