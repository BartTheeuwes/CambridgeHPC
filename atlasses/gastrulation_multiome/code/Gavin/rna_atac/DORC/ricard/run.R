#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
  io$script <- "/Users/ricard/gastrulation_multiome_10x/Gavin/DORC/ricard/ComputeDORC_ArchR.R"
} else {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
  io$script <- "/homes/ricard/gastrulation_multiome_10x/Gavin/DORC/ricard/ComputeDORC_ArchR.R"
}
io$outdir <- paste0(io$basedir,"/results/rna_atac/DORCs_v2"); dir.create(io$outdir, showWarnings = F)
io$tmpdir <- paste0(io$outdir,"/tmp"); dir.create(io$tmpdir, showWarnings = F)

# samples
opts$samples <- c("E7.5_rep1", "E7.5_rep2", "E8.0_rep1", "E8.0_rep2", "E8.5_rep1", "E8.5_rep2")

# Test mode (subsetting cells)?
opts$test_mode <- FALSE

# Number of highly variable genes
opts$distance <- c(1e4)
# opts$distance <- c(1e4, 2.5e4, 5e4)

# Number of cores per job
opts$ncores <- 1

# LSF params
opts$memory <- 30000


####################################################
## Run separately per stage, only embryonic cells ##
####################################################

for (x in opts$stages) {
  samples <- opts$samples[grep(x,opts$samples)]
  for (i in opts$distance) {
    
    # Define LSF command
    if (grepl("ricard",Sys.info()['nodename'])) {
      lsf <- ""
    } else if (grepl("ebi",Sys.info()['nodename'])) {
      lsf <- sprintf("bsub -M %s -n %d -o %s/DORC_%s_distance%d_%dcores.txt", opts$memory, opts$ncores, io$tmpdir, x, i, opts$ncores)
    }
    cmd <- sprintf("%s Rscript %s --remove_ExE_celltypes --samples %s --distance %d --ncores %d --outdir %s", 
                   lsf, io$script, paste(samples,collapse=" "), i, opts$ncores, io$outdir)
    
    if (isTRUE(opts$test_mode)) cmd <- paste0(cmd, " --test_mode")
    
    # Run
    # print(cmd)
    # system(cmd)
  }
}

#############################
## Run all embryonic cells ##
#############################

for (i in opts$distance) {
  
  # Define LSF command
  if (grepl("ricard",Sys.info()['nodename'])) {
    lsf <- ""
  } else if (grepl("ebi",Sys.info()['nodename'])) {
    lsf <- sprintf("bsub -M %s -n %d -o %s/DORC_allcells_distance%d_%dcores.txt", opts$memory, opts$ncores, io$tmpdir, i, opts$ncores)
  }
  cmd <- sprintf("%s Rscript %s --remove_ExE_celltypes --samples %s --distance %d --ncores %d --outdir %s", 
                 lsf, io$script, paste(opts$samples,collapse=" "), i, opts$ncores, io$outdir)
  
  if (isTRUE(opts$test_mode)) cmd <- paste0(cmd, " --test_mode")
  
  # Run
  print(cmd)
  system(cmd)
}

#####################################
## Run all cells (embryonic + ExE) ##
#####################################

for (i in opts$distance) {
  
  # Define LSF command
  if (grepl("ricard",Sys.info()['nodename'])) {
    lsf <- ""
  } else if (grepl("ebi",Sys.info()['nodename'])) {
    lsf <- sprintf("bsub -M %s -n %d -o %s/DORC_allcells_distance%d_%dcores.txt", opts$memory, opts$ncores, io$tmpdir, i, opts$ncores)
  }
  cmd <- sprintf("%s Rscript %s --samples %s --distance %d --ncores %d --outdir %s", 
                 lsf, io$script, paste(opts$samples,collapse=" "), i, opts$ncores, io$outdir)
  
  if (isTRUE(opts$test_mode)) cmd <- paste0(cmd, " --test_mode")
  
  # Run
  # print(cmd)
  # system(cmd)
}

