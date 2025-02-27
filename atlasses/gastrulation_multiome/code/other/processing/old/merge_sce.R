library(SingleCellExperiment)
library(scran)
library(batchelor)
library(scater)

a <- readRDS("/Users/ricard/data/gastrulation_multiome_10x/multiome1/processed/rna/SingleCellExperiment.rds")
colnames(a) <- colnames(a) %>% stringr::str_replace_all(.,"-1","") %>% paste0("multiome1_",.)

b <- readRDS("/Users/ricard/data/gastrulation_multiome_10x/multiome2/processed/rna/SingleCellExperiment.rds")

sample_metadata <- fread("/Users/ricard/data/gastrulation_multiome_10x/sample_metadata.csv")

a <- a[intersect(rownames(a), rownames(b)),]
b <- b[intersect(rownames(a), rownames(b)),]
# rowdata <- rowData(a)[,c("symbol","description")]

# b <- b[,colnames(b)%in%sample_metadata$id_rna]

# rowData(a) <- NULL
# colData(a) <- NULL
# rowData(b) <- NULL
# colData(b) <- NULL

ab <- SingleCellExperiment::SingleCellExperiment(
  list(counts=Matrix::Matrix(cbind(counts(a),counts(b)),sparse=TRUE)))

batch <- substr(colnames(ab), 1,9)
ab <- multiBatchNorm(ab, batch=batch)

coldata <- sample_metadata %>%
  .[cell %in% colnames(ab)] %>% 
  setkey(cell) %>% .[colnames(ab)] %>%
  tibble::column_to_rownames("cell")
stopifnot(all(colnames(ab)==rownames(coldata)))
colData(ab) <- DataFrame(coldata)

saveRDS(ab, "/Users/ricard/data/gastrulation_multiome_10x/processed/rna/SingleCellExperiment.rds")
