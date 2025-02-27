library(SnapATAC)
library(purrr)
library(data.table)


source(here::here("settings.R"))
source(here::here("atac/SnapATAC/snapatac_settings.R"))

opts$use_cell_types <- TRUE # run analysis per cell type (or per cluster)
opts$cores          <- 8



opts$fdr_cutoff <- 0.05

snap <- readRDS(snapio$rds_file)

if (opts$use_cell_types) {
  snap@cluster <- as.factor(snap@metaData$cell_type)
  snapio$difacc <- gsub(".tsv", "_celltype.tsv", snapio$difacc)
}


# iterate over clusters
clusters <- levels(snap@cluster)

# do differential testing
difacc <- map(clusters, ~{
  DARs <- findDAR(obj=snap,
                  input.mat="pmat",
                  cluster.pos=.x,
                  cluster.neg.method="knn",
                  test.method="exactTest",
                  bcv=0.1, #0.4 for human, 0.1 for mouse
                  seed.use=10)
  
  setDT(DARs)
  
  DARs[, cluster := .x]
  DARs[, row := .I]
  DARs[, FDR := p.adjust(DARs$PValue, method="BH")]
  setkey(DARs, FDR)
  DARs[, sig := FALSE]
  DARs[FDR < opts$fdr_cutoff & logFC > 0, sig := TRUE]
  
  DARs

}) %>% 
  rbindlist()


covs <- Matrix::rowSums(snap@pmat)

# plot in umap

walk(clusters, ~{
  idy <- difacc[cluster == .x & sig == TRUE, row]
  
  if (length(idy) < 2000L){
    idy <- difacc[cluster == .x & logFC > 0] %>%
      setorder("FDR") %>%
      .[1:2000, row]
  }
  
  (title   <- paste0("Differentially accessible peaks in cluster ", .x))
  (outfile <- paste0(snapio$plots_out, "/umap_difacc_cluster", .x, ".pdf"))
  
  
  
  vals        <- Matrix::rowSums(snap@pmat[,idy]) / covs
  vals.zscore <- (vals - mean(vals)) / sd(vals)
  
  plotFeatureSingle(
    obj=snap,
    feature.value=vals.zscore,
    method="umap", 
    main=title,
    point.size=0.1, 
    point.shape=19, 
    down.sample=10000,
    quantiles=c(0.01, 0.99),
    pdf.file.name = outfile
  )
})


# save diffacc sites

peaks <- setDT(as.data.frame(snap@peak)) %>%
  .[, .(chr = gsub("b'|'", "", seqnames),
        start,
        end)] %>%
  .[, row := .I]

difacc <- merge(difacc, peaks, by = "row") %>% 
  setorder("FDR")

fwrite(difacc, snapio$difacc, sep = "\t", na = "NA")