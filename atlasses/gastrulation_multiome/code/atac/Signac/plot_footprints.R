library(Seurat)
library(Signac)
library(purrr)
library(data.table)
library(JASPAR)
library(TFBSTools)



source(here::here("settings.R"))




io$signac        <- file.path(io$basedir, "/processed/atac/signac/signacFootprints.rds")

io$plots_out     <- "/bi/home/clarks/plots/10X_multiome/signac/footprints"


opts$motifs      <- "All" # vector of motifs from JASPAR c("Sox1", "Foxo1")

opts$lineages    <- c("Gut", "ExE_mesoderm", "Forebrain_Midbrain_Hindbrain", 
                        "Erythroid1")

opts$species     <- "mouse"

if (toupper(opts$species) == "HUMAN") opts$species <- 9606
if (toupper(opts$species) == "MOUSE") opts$species <- 10090


dir.create(io$plots_out, recursive = TRUE)

signac <- readRDS(io$signac)
signac



if (toupper(opts$motifs) == "ALL") {
  print("all motifs:")
  
  # extract position frequency matrices for the motifs
  pwm <- getMatrixSet(
    x = JASPAR,
    opts = list(species = opts$species, all_versions = FALSE)
  )
  
 
  

    opts$motifs <- name(pwm) %>%
      unname() %>%
      .[order(.)]

}

# plot the footprint data for each group of cells
plots <- map(opts$motifs, ~{
  print(.x)
  p <- PlotFootprint(signac, 
                     features = .x,
                     idents = opts$lineages)
  outfile <- paste0(io$plots_out, "/", .x, ".pdf") %>% gsub("::", "_", .)
  cowplot::save_plot(outfile, p)
})













