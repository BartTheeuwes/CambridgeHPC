
#####################
## Define settings ##
#####################

source(here::here("settings.R"))

# I/O
io$script <- here::here("rna_atac/rna_vs_chromvar/trajectories/rna_vs_chromvar_trajectory.R")
io$outdir <- paste0(io$basedir,"/results/rna_atac/rna_vs_chromvar/trajectories"); dir.create(io$outdir, showWarnings = F)
io$tmpdir <- paste0(io$outdir,"/tmp"); dir.create(io$tmpdir, showWarnings = F)

####################
## Define options ##
####################

# samples
opts$trajectories <- list(
  "blood_trajectory" = c("Haematoendothelial_progenitors", "Blood_progenitors_1", "Blood_progenitors_2", "Erythroid1", "Erythroid2", "Erythroid3"),
  # "ectoderm_trajectory" = c("Epiblast", "Rostral_neurectoderm", "Forebrain_Midbrain_Hindbrain"),
  "mesoderm_trajectory" = c("Epiblast", "Primitive_Streak", "Nascent_mesoderm"),
  "endoderm_trajectory" = c("Epiblast", "Primitive_Streak", "Def._endoderm", "Gut")
)


# Number of KNN for denoising
opts$knn <- c(25,50)

opts$motif_annotation <- c("Motif_cisbp") # Motif_JASPAR2020_human

########################
## Run (no denoising) ##
########################

for (i in names(opts$trajectories)) {
  trajectory <- sprintf("%s/results/rna/trajectories/%s/%s.txt.gz",io$basedir,i,i)
  outdir <- sprintf("%s/%s",io$outdir,i); dir.create(outdir, showWarnings = F)
  
  # Define LSF command
  if (grepl("ricard",Sys.info()['nodename'])) {
    lsf <- ""
  } else if (grepl("ebi",Sys.info()['nodename'])) {
    lsf <- sprintf("bsub -M 20000 -n 1 -o %s/rna_vs_chromvar_%s.txt",io$tmpdir, i)
  }
  cmd <- sprintf("%s Rscript %s --celltypes %s --trajectory %s --trajectory_name %s --motif_annotation %s --outdir %s", 
                 lsf, io$script, paste(opts$trajectories[[i]],collapse=" "), trajectory, i, opts$motif_annotation, outdir)

  # Run
  print(cmd)
  system(cmd)
}


#####################
## Run (denoising) ##
#####################

for (i in names(opts$trajectories)) {
  trajectory <- sprintf("%s/results/rna/trajectories/%s/%s.txt.gz",io$basedir,i,i)
  for (j in opts$knn) {
    outdir <- sprintf("%s/%s_knn%s",io$outdir,i,j); dir.create(outdir, showWarnings = F)

    # Define LSF command
    if (grepl("ricard",Sys.info()['nodename'])) {
      lsf <- ""
    } else if (grepl("ebi",Sys.info()['nodename'])) {
      lsf <- sprintf("bsub -M 25000 -n 1 -o %s/rna_vs_chromvar_%s_knn%d.txt", io$tmpdir, i,j)
    }
    cmd <- sprintf("%s Rscript %s --celltypes %s --trajectory %s --trajectory_name %s --denoise --knn %d --motif_annotation %s --outdir %s", 
                   lsf, io$script, paste(opts$trajectories[[i]],collapse=" "), trajectory, i, j, opts$motif_annotation, outdir)

    # Run
    print(cmd)
    system(cmd)
  }
}
