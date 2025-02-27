library(SoupX)
library(purrr)
library(data.table)

source(here::here("settings.R"))

# script to run the SoupX method to quantify and remove contamination from ambient RNA
# takes as input 10X CellRanger output files


# output from 10x CellRanger can be found in the following folders:

io$indirs <- file.path(
  io$rawdata, 
  c("multiome1",
    "multiome2",
    "rep1_L001_multiome",
    "rep2_L002_multiome"
    ), 
  "outs"
  )



io$outdirs <- gsub("outs", "soupX", io$indirs)
map(io$outdirs, dir.create)
.x=io$indirs[[1]]
.y=io$outdirs[[1]]

gene_metadata <- fread(io$gene_metadata) %>% 
  setnames("symbol", "gene")

# iterate over the different 10X samples 

walk2(io$indirs, io$outdirs, ~{
  
  clusters_dt <- file.path(.x, "analysis/clustering/gex/graphclust/clusters.csv") %>% 
    fread()
  clusters <- clusters_dt[, Cluster]
  names(clusters) <- clusters_dt[, Barcode]
  
  # run the SoupX commands
  sc <- load10X(.x)
  sc <- setClusters(sc, clusters)
  sc <- autoEstCont(sc, forceAccept = TRUE)
  adj <- adjustCounts(sc)
  
  # save adjusted data  
  Matrix::writeMM(adj, paste0(.y, "/soupX_adjusted_matrix.mtx.gz"))
  
  soup <- setDT(sc$soupProfile, keep.rownames = "gene") %>%
    merge(gene_metadata, by = "gene") %>% 
    .[order(-rank(counts))] 
  
  fwrite(soup, paste0(.y, "/soup.tsv.gz"), sep = "\t", na = "NA")
  
})
