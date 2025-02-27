library(SnapATAC)
library(purrr)
library(data.table)


source(here::here("settings.R"))
source(here::here("atac/SnapATAC/snapatac_settings.R"))

opts$use_cell_types <- TRUE # run analysis per cell type (or per cluster)
opts$cores          <- 8


io$homer_out   <- file.path(io$basedir, "/processed/atac/snapatac/homer")
io$genome_fa   <- "/bi/scratch/Stephen_Clark/annotations/genomes/Mus_musculus.GRCm38.dna.primary_assembly.fa"
  
  
dir.create(io$homer_out, recursive = TRUE)

opts$cores     <- 8
opts$motif_length <- 8
opts$scan_size <- 200
opts$optimize_count <- 2
opts$cache <- 1e2 # cache in MB
opts$fdr_num <- 2 #5


if (opts$use_cell_types) {
  #snap@cluster <- as.factor(snap@metaData$cell_type)
  snapio$difacc <- gsub(".tsv", "_celltype.tsv", snapio$difacc)
}

dir.create(io$homer_out, recursive = TRUE)

#snap <- readRDS(snapio$rds_file)

difacc <- fread(snapio$difacc) %>%
  .[sig == TRUE]



(homer <- system("which findMotifsGenome.pl", intern = TRUE))


# iterate over clusters
clusters <- difacc[, unique(cluster)]

.x=clusters[[1]]

#rows <- difacc[cluster == .x, row]
# peaks <- snap@peak[rows,] %>% 
#   as.data.frame() %>% 
#   setDT() %>%
# .[, .(
#   Chr     = gsub("b|'|chr", "", seqnames), 
#   Start   = start, 
#   End     = end,
#   #extrcol = ".",
#   PeakID  = gsub("b|'| ", "", gsub(":|-", "_", name)),
#   extrcol = ".")#,
#   # Strand  = "*")
#   ]


peaks <- difacc[cluster == .x] %>%
  .[, .(chr = gsub("chr", "", chr),
        start,
        end,

        id = paste0(chr, start, end),
        strand = "+")
    ]



peaks <- peaks[chr %in% c(1:19, "X", "Y")] #peaks[chr %like% "chr"]
peaks

# generate test peak file

#peaks <- data.table(chr = paste0("chr", c(1, "X")), start = 5e6, end = 5e6 + 500)


setorder(peaks, "chr", "start", "end")
peaks[, unique(chr)]
peakfile <- tempfile(fileext = ".txt")

fwrite(peaks, peakfile, sep = "\t", na = "NA", col.names = FALSE, row.names = FALSE, quote = FALSE)
fread(peakfile)
#system(paste("head", io$genome_fa))
outdir <- file.path(io$homer_out, paste0("cluster_", .x))
system(paste("rm -rd", outdir))

tmpdir <- tempdir()

cmd <- paste(
  #homer,
  "findMotifsGenome.pl",
  peakfile, 
  io$genome_fa,#"mm10", 
  outdir,
  #"-len", paste(opts$motif_length, collapse = ","),
  "-size", opts$scan_size,
  "-S", opts$optimize_count,
  "-p", opts$cores,
  "-cache", opts$cache,
  "-fdr", opts$fdr_num
 # "-preparsedDir", tmpdir
)
cmd
system(cmd)



# dir.create(io$homer_out, recursive = TRUE)
# 
# snap <- readRDS(io$snap_file)
# 
# difacc <- fread(io$difacc_file) %>%
#   .[sig == TRUE]
# 
# 
# 
# (homer <- system("which findMotifsGenome.pl", intern = TRUE))
# homer <- "/Users/clarks/opt/miniconda3/homer/.//bin/findMotifsGenome.pl"
# 
# # iterate over clusters
# clusters <- difacc[, unique(cluster)]
# 
# .x=clusters[1]
# 
#   outdir <- file.path(io$homer_out, paste0("cluster_", .x))
#   #dir.create(outdir)
#   system(paste("rm -r -d", outdir))
#   
#   rows <- difacc[cluster == .x, row]
#   
#   motifs <- runHomer(
#     snap[, rows,"pmat"], 
#     mat = "pmat",
#     path.to.homer = homer,
#     result.dir = outdir,
#     num.cores= opts$cores,
#     genome = 'mm10',
#     motif.length = 10,
#     scan.size = 300,
#     optimize.count = 2,
#     background = 'automatic',
#     local.background = FALSE,
#     only.known = TRUE,
#     only.denovo = FALSE,
#     fdr.num = 5,
#     cache = 100,
#     overwrite = TRUE,
#     keep.minimal = FALSE
#   )
#   
#   saveRDS(motifs, paste0(outdir, "/homer.rds"))
#   
#   
# 
# 
# 
# 
