library(purrr)
library(data.table)


source(here::here("settings.R"))

# subset fragments file for testing purposes


# io$bam_file   <- file.path(io$basedir, "public_datasets/Pijuan-Sala_2020/embryo_revision1_sorted.bam")
# io$sorted_bam <- file.path(io$basedir, "public_datasets/Pijuan-Sala_2020/embryo_revision1_possorted.bam")

io$fragments_file <- file.path(io$basedir, "public_datasets/Pijuan-Sala_2020/bam/embryo_revision1_fragments_sorted.tsv.gz")
io$outfile <- "/bi/scratch/Stephen_Clark/fragments1000.tsv.gz"




frags <- fread(cmd = paste("zcat", io$fragments_file))


# for testing - subset to 1000 cells
 cells <- frags[, unique(V4)][1:1000]
 frags <- frags[V4 %in% cells]

 fwrite(frags, io$outfile, sep = "\t", col.names = FALSE, quote = FALSE)
 
 