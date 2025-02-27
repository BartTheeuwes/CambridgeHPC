suppressPackageStartupMessages({
  library(SoupX)
  library(purrr)
  library(data.table)
  library(Matrix)
  library(furrr)
})

## run soupX pipeline from snakemake ##


#print(snakemake)


cores <- snakemake@threads[[1]]


if (cores > 1){
  plan(multisession, workers = cores)
}



# iterate over samples
n <- seq_along(snakemake@input[["inputdir"]])

future_map(n, ~{
  
  # load snakemake i/o
  indir         <- snakemake@input[["inputdir"]][[.x]]
  features_file <- snakemake@input[["features"]][[.x]]
  clusters_file <- snakemake@input[["clusters"]][[.x]]
  matfile       <- snakemake@output[["matrix"]][[.x]]
  soupfile      <- snakemake@output[["soup"]][[.x]]
  
  print(indir);print(matfile);print(soupfile)
  
  dir.create(dirname(matfile), recursive = TRUE)
  dir.create(dirname(soupfile), recursive = TRUE)
  
  # load gene metadata from cell ranger output folder
  gene_metadata <- fread(features_file) %>% 
    setnames(c("ens_id", "gene", "assay", "chr", "start", "end"))
  
  # read in clustering data from cell ranger output
  clusters_dt <- fread(clusters_file)
  clusters <- clusters_dt[, Cluster]
  names(clusters) <- clusters_dt[, Barcode]
  
  # run the SoupX commands
  sc <- load10X(indir) %>% 
    setClusters(clusters) %>% 
    autoEstCont(forceAccept = TRUE)
  
  # adjust counts (to remove soup)
  adj <- adjustCounts(sc)
  
  #save adjusted data  
  writeMM(adj, matfile)
  
  # save soup counts as tsv
  soup <- setDT(sc$soupProfile, keep.rownames = "gene") %>%
    merge(gene_metadata, by = "gene", all.x = TRUE) %>% 
    .[order(-rank(counts))] 
  
  fwrite(soup, soupfile, sep = "\t", na = "NA", quote = FALSE)

})


