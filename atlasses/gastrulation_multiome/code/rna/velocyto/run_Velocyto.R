library(data.table)
library(purrr)
library(rhdf5)
library(Matrix)

# run velocyto then convert the loom file to an RDS of sparse matrices

source(here::here("settings.R"))



io$indirs     <- file.path(
  io$rawdata, c(
    "multiome1", 
    "multiome2", 
    "rep1_L001_multiome", 
    "rep2_L002_multiome"
    )
)


io$outfiles    <- file.path(
  io$indirs, 
  "processed/velocyto", 
  paste0(basename(io$indirs), ".rds")
  )

io$gtf_file   <- "/bi/scratch/Stephen_Clark/annotations/gtf/Mus_musculus.GRCm38.98.gtf"

opts$cores    <- 8
opts$mem      <- 1e4

map(dirname(io$outfiles), dir.create, recursive = TRUE)

.x=io$indirs[2]; .y=io$outfiles[2]

walk2(io$indirs, io$outfiles, ~{
  
  cmd <- paste("velocyto run10x",
               "--samtools-threads", opts$cores,
               "--samtools-memory", opts$mem,
               .x,
               io$gtf_file)
  cmd
  system(cmd)
  
  loom <- dir(file.path(.x, "velocyto"), pattern = ".loom$", full = TRUE)
  
  h5ls(loom)
  cells <- h5read(loom, "col_attrs")
  genes <- h5read(loom, "row_attrs")
  mats <- h5read(loom, "layers")
  
  cells <- gsub(":", "_", cells[[1]]) %>% gsub("x", "", .)
  head(cells)
  
  genes <- genes$Accession
  head(genes)
  
  
  # matrices do not have dim names and need transposing to genes x cells
  mats <- map(mats, function(mat){
    mat <- t(mat)
    colnames(mat) <- cells
    rownames(mat) <- genes
    Matrix(mat, sparse = TRUE)
  })
  
  saveRDS(mats, .y)
})