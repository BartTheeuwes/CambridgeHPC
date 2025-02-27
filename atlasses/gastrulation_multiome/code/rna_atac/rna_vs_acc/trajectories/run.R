
################
## Define I/O ##
################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
  io$script <- "/Users/ricard/gastrulation_multiome_10x/rna_atac/rna_vs_acc/trajectories/run_TFexpr_vs_peakAcc_trajectories.R"
} else {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
  io$script <- "/homes/ricard/gastrulation_multiome_10x/rna_atac/rna_vs_acc/trajectories/run_TFexpr_vs_peakAcc_trajectories.R"
}
io$outdir <- paste0(io$basedir,"/results/rna_atac/rna_vs_acc/trajectories"); dir.create(io$outdir, showWarnings = F)
io$tmpdir <- paste0(io$outdir,"/tmp"); dir.create(io$tmpdir, showWarnings = F)

####################
## Define options ##
####################

# trajectories
opts$trajectories <- c("blood_trajectory", "ectoderm_trajectory", "mesoderm_trajectory", "endoderm_trajectory")
# opts$trajectories <- "blood_trajectory"

# Number of KNN for denoising
opts$knn <- c(50)

# Motif annotation
opts$motif_annotation <- c("Motif_cisbp") 

# Test mode
opts$test_mode <- FALSE

#########
## Run ##
#########

for (i in opts$trajectories) {
  trajectory <- sprintf("%s/results/rna/trajectories/%s/%s.txt.gz",io$basedir,i,i)
  for (j in opts$knn) {
    outdir <- sprintf("%s/%s_knn%s",io$outdir,i,j); dir.create(outdir, showWarnings = F)

    # Define LSF command
    if (grepl("ricard",Sys.info()['nodename'])) {
      lsf <- ""
    } else if (grepl("ebi",Sys.info()['nodename'])) {
      lsf <- sprintf("bsub -M 25000 -n 1 -o %s/rna_vs_acc_%s_knn%d.txt", io$tmpdir, i,j)
    }
    cmd <- sprintf("%s Rscript %s --trajectory %s --trajectory_name %s --denoise --knn %d --motif_annotation %s --outdir %s", 
                   lsf, io$script, trajectory, i, j, opts$motif_annotation, outdir)
    if (isTRUE(opts$test_mode)) cmd <- paste0(cmd, " --test_mode")
    
    # Run
    print(cmd)
    system(cmd)
  }
}
