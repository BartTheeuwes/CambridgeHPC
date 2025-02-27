
################
## Define I/O ##
################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
  io$script <- "/Users/ricard/gastrulation_multiome_10x/rna/dimensionality_reduction/automated/dimensionality_reduction_sce.R"
} else {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
  io$script <- "/homes/ricard/gastrulation_multiome_10x/rna/dimensionality_reduction/automated/dimensionality_reduction_sce.R"
  io$tmpdir <- paste0(io$basedir,"/results/rna/dimensionality_reduction/tmp"); dir.create(io$tmpdir)
}
io$outdir <- paste0(io$basedir,"/results/rna/dimensionality_reduction")

####################
## Define options ##
####################

opts$samples <- c(
  "E7.5_rep1",
  "E7.5_rep2",
  "E8.0_rep1",
  "E8.0_rep2",
  "E8.5_rep1",
  "E8.5_rep2"
)


# Test mode (subsetting cells)?
opts$test_mode <- FALSE

# Number of highly variable genes
# opts$features <- c(1000,2000,3000)
opts$features <- c(2500)

# Number of PCs
# opts$npcs <- c(25,50)
opts$npcs <- c(30)

# Variables to regress out
opts$vars.to.regress <- c("nFeature_RNA","mitochondrial_percent_RNA")

# Variable to do MNN batch correction on
opts$batch.correction <- c("sample")

# UMAP hyperparameters
# opts$n_neighbors <- c(20,30,40)
opts$n_neighbors <- c(25)

# opts$min_dist <- c(0.20,0.30,0.40)
opts$min_dist <- c(0.30)

# LSF params
opts$memory <- 9000


############################################################
## Run one sample at a time (no batch correction applied) ##
############################################################

opts$colour_by <- c("celltype.mapped","nFeature_RNA","hybrid_score","hybrid_call")

for (i in opts$samples) {
  outdir <- sprintf("%s/%s",io$outdir,i); dir.create(outdir, showWarnings = F)
  for (j in opts$features) {
    for (k in opts$npcs) {

      # Define LSF command
      if (grepl("ricard",Sys.info()['nodename'])) {
        lsf <- ""
      } else if (grepl("ebi",Sys.info()['nodename'])) {
        lsf <- sprintf("bsub -M %s -n 1 -o %s/%s_%d_%d.txt", opts$memory, io$tmpdir,i,j,k)
      }
      cmd <- sprintf("%s Rscript %s --samples %s --features %d --npcs %d --vars.to.regress %s --n_neighbors %s --min_dist %s --colour_by %s --outdir %s",
                     lsf, io$script, i, j, k, paste(opts$vars.to.regress,collapse=" "), paste(opts$n_neighbors,collapse=" "), paste(opts$min_dist,collapse=" "), paste(opts$colour_by,collapse=" "), outdir)

      if (isTRUE(opts$test_mode)) cmd <- paste0(cmd, " --test")

      # Run
      # print(cmd)
      # system(cmd)
    }
  }
}


############################################################
## Run one stage at a time (no batch correction applied) ##
############################################################

opts$colour_by <- c("celltype.mapped","sample","nFeature_RNA","hybrid_score","hybrid_call")

for (i in opts$stages) {
  samples <- opts$samples[grep(i,opts$samples)]
  outdir <- sprintf("%s/%s",io$outdir,i); dir.create(outdir, showWarnings = F)
  for (j in opts$features) {
    for (k in opts$npcs) {
      
      # Define LSF command
      if (grepl("ricard",Sys.info()['nodename'])) {
        lsf <- ""
      } else if (grepl("ebi",Sys.info()['nodename'])) {
        lsf <- sprintf("bsub -M %s -n 1 -o %s/%s_%d_%d.txt", opts$memory, io$tmpdir,i,j,k)
      }
      cmd <- sprintf("%s Rscript %s --samples %s --features %d --npcs %d --vars.to.regress %s --n_neighbors %s --min_dist %s --colour_by %s --outdir %s", 
                     lsf, io$script, paste(samples,collapse=" "), j, k, paste(opts$vars.to.regress,collapse=" "), paste(opts$n_neighbors,collapse=" "), paste(opts$min_dist,collapse=" "), paste(opts$colour_by,collapse=" "), outdir)
      
      if (isTRUE(opts$test_mode)) cmd <- paste0(cmd, " --test")
      
      # Run
      # print(cmd)
      # system(cmd)
    }
  }
}

#############################
## Run all samples at once ##
#############################

outdir <- sprintf("%s/all_cells",io$outdir); dir.create(outdir, showWarnings = F)
opts$colour_by <- c("celltype.mapped","sample","stage","nFeature_RNA","hybrid_score","hybrid_call")

for (j in opts$features) {
  for (k in opts$npcs) {

    # Define LSF command
    if (grepl("ricard",Sys.info()['nodename'])) {
      lsf <- ""
    } else if (grepl("ebi",Sys.info()['nodename'])) {
      lsf <- sprintf("bsub -M %s -n 1 -o %s/allsamples_%d_%d.txt", 2*opts$memory, io$tmpdir,j,k)
    }
    cmd <- sprintf("%s Rscript %s --samples %s --features %d --npcs %d --vars.to.regress %s --batch.correction %s --n_neighbors %s --min_dist %s --colour_by %s --outdir %s",
                   lsf, io$script, paste(opts$samples,collapse=" "), j, k, paste(opts$vars.to.regress,collapse=" "), opts$batch.correction, paste(opts$n_neighbors,collapse=" "), paste(opts$min_dist,collapse=" "), paste(opts$colour_by,collapse=" "), outdir)

    if (isTRUE(opts$test_mode)) cmd <- paste0(cmd, " --test")

    # Run
    # print(cmd)
    # system(cmd)
  }
}


############################################
## Run all samples at once (no ExE cells) ##
############################################

outdir <- sprintf("%s/all_cells_noExE",io$outdir); dir.create(outdir, showWarnings = F)
opts$colour_by <- c("celltype.mapped","sample","stage","nFeature_RNA","hybrid_score","hybrid_call")

for (j in opts$features) {
  for (k in opts$npcs) {

    # Define LSF command
    if (grepl("ricard",Sys.info()['nodename'])) {
      lsf <- ""
    } else if (grepl("ebi",Sys.info()['nodename'])) {
      lsf <- sprintf("bsub -M %s -n 1 -o %s/allsamples_noExE_%d_%d.txt", 2*opts$memory, io$tmpdir,j,k)
    }
    cmd <- sprintf("%s Rscript %s --samples %s --features %d --remove_ExE_cells --npcs %d --vars.to.regress %s --batch.correction %s --n_neighbors %s --min_dist %s --colour_by %s --outdir %s",
                   lsf, io$script, paste(opts$samples,collapse=" "), j, k, paste(opts$vars.to.regress,collapse=" "), opts$batch.correction, paste(opts$n_neighbors,collapse=" "), paste(opts$min_dist,collapse=" "), paste(opts$colour_by,collapse=" "), outdir)

    if (isTRUE(opts$test_mode)) cmd <- paste0(cmd, " --test")

    # Run
    print(cmd)
    system(cmd)
  }
}