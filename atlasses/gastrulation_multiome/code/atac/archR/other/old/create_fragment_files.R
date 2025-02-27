features <- fread("/Users/ricard/data/gastrulation_multiome_10x/multiome2/original/filtered_feature_bc_matrix/features.tsv.gz") %>%
  .[,id:=1:.N]

rna.features <- features[V3=="Gene Expression",id]
atac.features <- features[V3=="Peaks",id]

barcodes <- fread("/Users/ricard/data/gastrulation_multiome_10x/multiome2/original/filtered_feature_bc_matrix/barcodes.tsv.gz") %>%
  .[,id:=1:.N]

matrix <- fread("/Users/ricard/data/gastrulation_multiome_10x/multiome2/original/filtered_feature_bc_matrix/matrix.mtx.gz", skip=4, quote="", sep=" ") %>% 
  setnames(c("feature","barcode","value"))

matrix.atac <- matrix[feature>=min(atac.features) & feature<=max(atac.features)]
matrix.rna <- matrix[feature>=min(rna.features) & feature<=max(rna.features)]


length(unique(matrix$feature))
length(unique(matrix$barcode))

# Save fragments files for ATAC. Format:
# chr1    10079   10209   TAGGAGGGTTAACCGT-1      1
# chr1    10090   10278   CCTTGGTAGTTCTCCC-1      5