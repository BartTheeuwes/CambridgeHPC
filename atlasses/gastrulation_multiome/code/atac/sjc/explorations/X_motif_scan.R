library(purrr)
library(data.table)
# library(motifmatchr)
library(JASPAR)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)

# scan whole genome for motif positions

source(here::here("settings.R"))



io$motif_bed_out <- file.path(io$rawdata, "/processed/atac/signac/chromvar/motifbedchrY.tsv.gz")



opts$cores       <- 8
opts$testing     <- TRUE


opts$species     <- "human"

if (toupper(opts$species) == "HUMAN") opts$species <- 9606
if (toupper(opts$species) == "MOUSE") opts$species <- 10090




# extract position frequency matrices for the motifs
pfm <- getMatrixSet(
  x = JASPAR,
  opts = list(species = opts$species, all_versions = FALSE)
)



pwm <- toPWM(pfm)




pfm_names <- as.data.table(name(pfm), keep.rownames = TRUE) %>%
  setnames(c("ID", "TF")) %>%
  .[, TF := gsub("::|\\(|\\)|\\.", "_", TF) %>% gsub("_$", "", .)] %>%
  .[, rows := ID]

pfm_names




# iterate over chromosomes and over motifs

chrs <- seqnames(BSgenome.Mmusculus.UCSC.mm10)
motifs <- names(pwm)

if (opts$testing) {
  chrs <- c("chrY")
  motifs <- motifs[1:3]
}


tmpdir <- file.path(dirname(io$motif_bed_out), "tmp")

allfiles <- map(chrs, function(c) {
  print(paste("searching chromosome", c))
  
  map(motifs, function(p) {
    print(p)
    tmpfile <- file.path(tmpdir, c, p, "tmp.tsv.gz")
    # check for previously generated file
    if (file.exists(tmpfile)) return(tmpfile)
    
    dir.create(dirname(tmpfile), recursive = TRUE)
    
    dt <- searchSeq(pwm[[p]],
              BSgenome.Mmusculus.UCSC.mm10[[c]],
              mc.cores = opts$cores) %>% 
      as.data.frame() %>% 
      setDT() %>% 
      .[, .(chr = c, start, end, ID, TF)]
    
    fwrite(dt, tmpfile, sep = "\t", na = "NA", quote = FALSE)
    tmpfile
  }) 
}) 

allfiles <- unlist(allfiles)

# concatenate tmp files
positions <- map(allfiles, fread) %>% rbindlist()

fwrite(positions, io$motif_bed_out, sep = "\t", na = "NA", quote = FALSE)

# delete tmp files
system(paste("rm -rd", tmpdir))
