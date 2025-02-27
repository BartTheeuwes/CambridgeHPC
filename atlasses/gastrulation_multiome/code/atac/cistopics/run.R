
################
## Define I/O ##
################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
  io$script <- "/Users/ricard/gastrulation_multiome_10x/atac/cistopics/run_cisTopic.R"
} else {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
  io$script <- "/homes/ricard/gastrulation_multiome_10x/atac/cistopics/run_cisTopic.R"
  io$tmpdir <- paste0(io$basedir,"/results/atac/cistopic/tmp"); dir.create(io$tmpdir)
}
io$outdir <- paste0(io$basedir,"/results/atac/cistopic")

####################
## Define options ##
####################

# samples
opts$samples <- c(
  "E7.5_rep1",
  "E7.5_rep2",
  "E8.0_rep1",
  "E8.0_rep2",
  "E8.5_rep1",
  "E8.5_rep2"
)
# opts$samples <- c("E7.5_rep1", "E7.5_rep2")

# Number of highly variable peaks
opts$nfeatures <- c(1e4)
# opts$nfeatures <- c(1e4,5e4,1e5)

# Number of topics
# opts$ntopics <- seq(from=20,to=60,by=10)
opts$ntopics <- 50

# Variable to do MNN batch correction on
# opts$batch.variable <- "stage"

# UMAP hyperparameter: number of neighbours
# opts$n_neighbors <- c(20,30,40)
opts$n_neighbors <- c(25,50)

# UMAP hyperparameter: minimum distance
# opts$min_dist <- c(0.20,0.30,0.40)
opts$min_dist <- c(0.30,0.45)

# Number of cores
opts$ncores <- 1

# Test mode (subsetting cells)?
opts$test_mode <- FALSE


##############################
## Run one sample at a time ##
##############################

# opts$colour_by <- c("celltype.predicted")
# opts$memory <- 15000

# for (i in opts$samples) {
#   outdir <- sprintf("%s/%s",io$outdir,i); dir.create(outdir, showWarnings = F)
#   for (j in opts$nfeatures) {
#     for (k in opts$ntopics) {
      
#       # Define LSF command
#       if (grepl("ricard",Sys.info()['nodename'])) {
#         lsf <- ""
#       } else if (grepl("ebi",Sys.info()['nodename'])) {
#         lsf <- sprintf("bsub -M %s -n 1 -o %s/%s_%d_%d.txt", opts$memory, io$tmpdir,i,j,k)
#       }
#       cmd <- sprintf("%s Rscript %s --samples %s --nfeatures %d --ntopics %d --n_neighbors %s --min_dist %s --colour_by %s --outdir %s",
#                      lsf, io$script, i, j, k, paste(opts$n_neighbors,collapse=" "), paste(opts$min_dist,collapse=" "), paste(opts$colour_by,collapse=" "), outdir)
      
#       # Run
#       print(cmd)
#       system(cmd)
#     }
#   }
# }



###########################################################
## Run one stage at a time (no batch correction applied) ##
###########################################################

# opts$colour_by <- c("celltype.predicted","sample")
# opts$memory <- 25000

# for (i in opts$stages) {
#   samples <- opts$samples[grep(i,opts$samples)]
#   outdir <- sprintf("%s/%s",io$outdir,i); dir.create(outdir, showWarnings = F)
#   for (j in opts$nfeatures) {
#     for (k in opts$ntopics) {

#       # Define LSF command
#       if (grepl("ricard",Sys.info()['nodename'])) {
#         lsf <- ""
#       } else if (grepl("ebi",Sys.info()['nodename'])) {
#         lsf <- sprintf("bsub -M %s -n 1 -o %s/%s_%d_%d.txt", opts$memory, io$tmpdir,i,j,k)
#       }
#       cmd <- sprintf("%s Rscript %s --samples %s --nfeatures %d --ntopics %d --n_neighbors %s --min_dist %s --colour_by %s --outdir %s",
#                      lsf, io$script, paste(samples,collapse=" "), j, k, paste(opts$n_neighbors,collapse=" "), paste(opts$min_dist,collapse=" "), paste(opts$colour_by,collapse=" "), outdir)

#       # Run
#       print(cmd)
#       system(cmd)
#     }
#   }
# }

##################################################
## Run all samples at once, no batch correction ##
##################################################

outdir <- sprintf("%s/all_cells",io$outdir); dir.create(outdir, showWarnings = F)
opts$colour_by <- c("celltype.predicted","stage")
opts$memory <- 60000


for (j in opts$nfeatures) {
  for (k in opts$ntopics) {
    
    # Define LSF command
    if (grepl("ricard",Sys.info()['nodename'])) {
      lsf <- ""
    } else if (grepl("ebi",Sys.info()['nodename'])) {
      lsf <- sprintf("bsub -M %s -n 1 -o %s/allsamples_%d_%d.txt", opts$memory, io$tmpdir,j,k)
    }
    cmd <- sprintf("%s Rscript %s --samples %s --nfeatures %d --ntopics %d --n_neighbors %s --min_dist %s --colour_by %s --outdir %s",
                   lsf, io$script, paste(opts$samples,collapse=" "), j, k, paste(opts$n_neighbors,collapse=" "), paste(opts$min_dist,collapse=" "), paste(opts$colour_by,collapse=" "), outdir)
    
    # Run
    print(cmd)
    system(cmd)
  }
}

################################################################
## Run all samples at once, no batch correction, no ExE cells ##
################################################################

outdir <- sprintf("%s/all_cells",io$outdir); dir.create(outdir, showWarnings = F)
opts$colour_by <- c("celltype.predicted","stage")
opts$memory <- 60000


for (j in opts$nfeatures) {
  for (k in opts$ntopics) {
    
    # Define LSF command
    if (grepl("ricard",Sys.info()['nodename'])) {
      lsf <- ""
    } else if (grepl("ebi",Sys.info()['nodename'])) {
      lsf <- sprintf("bsub -M %s -n 1 -o %s/allsamples_noExE_%d_%d.txt", opts$memory, io$tmpdir,j,k)
    }
    cmd <- sprintf("%s Rscript %s --remove_ExE_celltypes --samples %s --nfeatures %d --ntopics %d --n_neighbors %s --min_dist %s --colour_by %s --outdir %s",
                   lsf, io$script, paste(opts$samples,collapse=" "), j, k, paste(opts$n_neighbors,collapse=" "), paste(opts$min_dist,collapse=" "), paste(opts$colour_by,collapse=" "), outdir)
    
    # Run
    print(cmd)
    system(cmd)
  }
}


##############################################################
## Run all samples at once, MNN batch correction per sample ##
##############################################################

# outdir <- sprintf("%s/all_cells_mnn",io$outdir); dir.create(outdir, showWarnings = F)
# opts$colour_by <- c("celltype.predicted","sample","stage","log_nFrags_atac","doublet_call")

# for (j in opts$nfeatures) {
#   for (k in opts$ntopics) {

#     # Define LSF command
#     if (grepl("ricard",Sys.info()['nodename'])) {
#       lsf <- ""
#     } else if (grepl("ebi",Sys.info()['nodename'])) {
#       lsf <- sprintf("bsub -M %s -n 1 -o %s/allsamples_mnn_%d_%d.txt", opts$memory, io$tmpdir,j,k)
#     }
#     cmd <- sprintf("%s Rscript %s --samples %s --nfeatures %d --ntopics %d --batch.variable %s --batch.method MNN --n_neighbors %s --min_dist %s --colour_by %s --outdir %s",
#                    lsf, io$script, paste(opts$samples,collapse=" "), j, k, opts$batch.variable, paste(opts$n_neighbors,collapse=" "), paste(opts$min_dist,collapse=" "), paste(opts$colour_by,collapse=" "), outdir)

#     # Run
#     print(cmd)
#     system(cmd)
#   }
# }