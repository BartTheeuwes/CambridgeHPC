library(SnapATAC)
library(purrr)
library(data.table)
library(GenomicRanges)


source(here::here("settings.R"))
source(here::here("public_datasets/Pijuan-Sala_2020/snapatac/snapatac_settings.R"))


opts$cores          <- 8
opts$fdr_cutoff     <- 0.05
opts$min_sites      <- 100 # only run motif enrichment if a cluster has at least this many difacc sites

difacc <- fread(snapio$difacc)

snap <- readRDS(snapio$rds_file)


homer_path <- system("which findMotifsGenome.pl", intern = TRUE)



clusters <- difacc[, unique(cluster)]
.x=clusters[2]




# reformat chromosome column - otherwise homer won't work!!!
seqlevels(snap@peak) <- gsub("b'|'", "", seqlevels(snap@peak))




# iterate over clusters
motifs_by_cluster <- map(clusters, ~{
  
  sites <- difacc[FDR < opts$fdr_cutoff & cluster == .x]
  
  if (nrow(sites) < opts$min_sites) return(NULL)
  
  motifs <- runHomer(
    snap[, sites$row, "pmat"], 
    mat = "pmat",
    path.to.homer = homer_path,
    result.dir = file.path(tempdir(), .x),
    num.cores = opts$cores,
    genome = 'mm10',
    motif.length = 10,
    scan.size = 300,
    optimize.count = 2,
    background = 'automatic',
    local.background = FALSE,
    only.known = TRUE,
    only.denovo = FALSE,
    fdr.num = 5,
    cache = 100,
    overwrite = TRUE,
    keep.minimal = FALSE
  )
  motifs <- as.data.table(motifs)
  motifs[, cluster := .x]
  motifs
}) %>% 
  purrr::compact() %>%
  rbindlist()

motifs_by_cluster

fwrite(motifs_by_cluster, snapio$motifs, sep = "\t", na = "NA")




## try to re-make the runHomer function so it works with user specified genome....  

# sites <- difacc[FDR < opts$fdr_cutoff & cluster == .x]
# bed <- as.data.table(snap[, sites$row, "pmat"]@peak)[, 1:3]
# target_bed <- tempfile(fileext = ".bed")
# fwrite(bed, target_bed, sep = "\t", na = "NA", col.names = FALSE, quote = FALSE)
# result.dir <- tempdir()
# 
# 
# cmd <- paste(
#   homer_path,
#   target_bed, 
#   genome, 
#   result.dir,
#   "-len", paste0(motif.length, collapse = ","),
#   "-size", scan.size,
#   "-S", optimize.count,
#   "-p", num.cores,
#   "-cache", cache,
#   "-fdr", fdr.num
# )



