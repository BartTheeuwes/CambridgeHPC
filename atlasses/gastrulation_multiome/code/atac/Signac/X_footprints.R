library(Seurat)
library(Signac)
library(purrr)
library(data.table)
library(future)

library(motifmatchr)
library(JASPAR)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)



source(here::here("settings.R"))



io$signac        <- file.path(io$basedir, "/processed/atac/signac/signac.rds")
io$signac_out    <- file.path(io$basedir, "/processed/atac/signac/signacFootprints_human.rds")

io$plots_out     <- "/bi/home/clarks/plots/10X_multiome/signac/footprints"


opts$cores       <- 8
opts$mem         <- 100 # GB

opts$motifs      <- "All" # vector of motifs from JASPAR c("Sox1", "Foxo1")
opts$lineages      <- c("Gut", "ExE_mesoderm", "Forebrain_Midbrain_Hindbrain", 
                        "Erythroid1")


opts$species     <- "human"

if (toupper(opts$species) == "HUMAN") opts$species <- 9606
if (toupper(opts$species) == "MOUSE") opts$species <- 10090


plan("multiprocess", workers = opts$cores)
options(future.globals.maxSize = opts$mem * 1024 ^ 3)
plan()

signac <- readRDS(io$signac)
signac
signac[["bins"]] <- NULL
signac


# extract position frequency matrices for the motifs
pwm <- getMatrixSet(
  x = JASPAR,
  opts = list(species = opts$species, all_versions = FALSE)
)

pwm

if (toupper(opts$motifs) == "ALL"){
  opts$motifs <- name(pwm) %>%
    unname() %>%
    .[order(.)]
}

print("scanning for motifs...")

# scan the genome and record the position of each motif using motifmatchr
motif.positions <- matchMotifs(
  pwms = pwm,
  subject = granges(signac),
  out = 'positions',
  genome = 'mm10'
)

print("adding motif object to Seurat...")

# create a Motif object and add it to the assay
motif <- CreateMotifObject(
  positions = motif.positions,
  pwm = pwm
)

signac <- SetAssayData(
  object = signac,
  slot = 'motifs',
  new.data = motif
)

print("extracting footprints...")

# gather the footprinting information for sets of motifs
signac <- Footprint(
  object = signac,
  motif.name = opts$motifs,
  genome = BSgenome.Mmusculus.UCSC.mm10
)

print("saving...")

saveRDS(signac, io$signac_out)

print("plotting...")

.x=opts$motifs[1]
# plot the footprint data for each group of cells
plots <- map(opts$motifs, ~{
  print(.x)
  p <- PlotFootprint(signac, 
                     features = .x,
                     idents = opts$lineages)
  outfile <- paste0(io$plots_out, "/", .x, ".pdf")
  cowplot::save_plot(outfile, p)
})













