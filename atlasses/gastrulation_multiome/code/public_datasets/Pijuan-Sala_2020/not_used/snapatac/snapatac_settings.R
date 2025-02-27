# script containing i/o for snapatac scripts
source(here::here("settings.R"))


print("loading settings....")

snapio <- list()

#snapio$frags_file <- file.path(io$basedir,  "public_datasets/Pijuan-Sala_2020/data/fragments.tsv.gz")
# for testing
snapio$frags_file <- "/bi/scratch/Stephen_Clark/gastrulation_multiome_10x/public_datasets/Pijuan-Sala_2020/bam//embryo_revision1_fragments.tsv"

snapio$snap_file  <- file.path(io$basedir,  "public_datasets/Pijuan-Sala_2020/snapatac/snap.snap")

snapio$rds_file   <- gsub(".snap$", ".rds", snapio$snap_file) #file.path(io$basedir,  "public_datasets/Pijuan-Sala_2020/snapatac/sub1000/snap.rds")
snapio$peaks      <- gsub(".snap$", "_peaks.tsv.gz", snapio$snap_file) #file.path(io$basedir,  "public_datasets/Pijuan-Sala_2020/snapatac/sub1000/peaks.tsv.gz")
snapio$difacc     <- gsub(".snap$", "_difacc.tsv.gz", snapio$snap_file) #file.path(io$basedir,  "public_datasets/Pijuan-Sala_2020/snapatac/sub1000/difacc.tsv.gz")
snapio$motifs     <- gsub(".snap$", "_motifs.tsv.gz", snapio$snap_file) #file.path(io$basedir,  "public_datasets/Pijuan-Sala_2020/snapatac/sub1000/motifs.tsv.gz")

snapio$metadata   <- file.path(io$basedir,  "public_datasets/Pijuan-Sala_2020/cell_metadata.csv")

snapio$plots_out             <- "/bi/home/clarks/plots/10X_multiome/public_datasets/Pijuan-Sala_2020/snapatac"
dir.create(snapio$plots_out, recursive = TRUE)

dir_create <- function(path){
  if (!dir.exists(path)) {
    print("creating directory")
    print(path)
    dir.create(path, recursive = TRUE)
  }
}
walk(snapio, ~dir_create(dirname(.x)))
