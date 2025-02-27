library(data.table)
library(purrr)
library(ggplot2)
library(cowplot)
library(Matrix)

source(here::here("settings.R"))

io$soup_files <- c(
  file.path(io$rawdata, "/multiome1/soupX/soup.tsv.gz"),         
  file.path(io$rawdata, "/multiome2/soupX/soup.tsv.gz"),        
  file.path(io$rawdata, "/rep1_L001_multiome/soupX/soup.tsv.gz"),
  file.path(io$rawdata, "/rep2_L002_multiome/soupX/soup.tsv.gz")
  )

io$matrix_files <- c(
  file.path(io$rawdata, "/multiome1/outs/filtered_feature_bc_matrix/matrix.mtx.gz"),
  file.path(io$rawdata, "/multiome2/outs/filtered_feature_bc_matrix/matrix.mtx.gz"),
  file.path(io$rawdata, "/rep1_L001_multiome//outs/filtered_feature_bc_matrix/matrix.mtx.gz"),
  file.path(io$rawdata, "/rep2_L002_multiome//outs/filtered_feature_bc_matrix/matrix.mtx.gz")
  
  
)
file.exists(io$matrix_files)

io$rna_features <- c(
  file.path(io$rawdata, "multiome1/outs/filtered_feature_bc_matrix/features.tsv.gz"),
  file.path(io$rawdata, "multiome2/outs/filtered_feature_bc_matrix/features.tsv.gz"),
  file.path(io$rawdata, "rep1_L001_multiome//outs/filtered_feature_bc_matrix/features.tsv.gz"),
  file.path(io$rawdata, "rep2_L002_multiome//outs/filtered_feature_bc_matrix/features.tsv.gz")
  
)
file.exists(io$rna_features)

io$merged_metrics_file    <- file.path(io$basedir, "/processed/atac/signac/cell_metrics.tsv")
.x=io$soup_files[[1]]

soup <- map(io$soup_files, ~{
  sample <- basename(dirname(dirname(.x)))
  fread(.x) %>% 
    .[, sample := sample]
  
}) %>% 
  rbindlist()


soup[, totals := sum(counts), sample]
soup[, counts_by_chr := sum(counts), .(chr, sample)]

mt <- soup[chr == "chrMT"][, .(mt_fraction = counts_by_chr / totals), sample] %>% 
  unique()
mt
ggplot(mt, aes(sample, mt_fraction, colour = sample, fill = sample)) +
  geom_bar(alpha = 0.5, stat = "identity") +
  theme_cowplot() +
  guides(fill = FALSE, colour = FALSE) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle("MT reads in the soup")
mt


metrics <- fread(io$merged_metrics_file) %>% 
  .[, mt_fraction := atac_mitochondrial_reads / atac_raw_reads]
.x=io$matrix_files[[1]]
.y=io$rna_features[[1]]

rna <- map2(io$matrix_files, io$rna_features, ~{
  features <- fread(.y) %>% 
    .[, row := .I]
  rows <- features[V2 %like% "mt-", row]
  mat <- readMM(.x)
  totals <- colSums(mat)
  mito <- colSums(mat[rows, ])
  
  data.table(mt_rna = mean(mito/totals),
             sample = basename(dirname(dirname(dirname(.x))))
})

rna

