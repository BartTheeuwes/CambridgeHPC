library(data.table)
library(purrr)

metadata <- fread("/Users/argelagr/data/gastrulation_multiome_10x/test/results/atac/archR/celltype_assignment/sample_metadata_after_celltype_assignment_fixed.txt.gz") %>%
  .[pass_rnaQC==T & pass_atacQC==TRUE & doublet_call==FALSE & !is.na(celltype.mapped_mnn)] %>%
  .[,c("cell", "sample", "stage", "nFeature_RNA", "nFrags_atac", "pass_atacQC", "pass_rnaQC",  "celltype.mapped_mnn")] %>%
  setnames("celltype.mapped_mnn","celltype")


fwrite(metadata, "~/Downloads/cell_metadata.txt.gz", quote=F, na="NA", sep="\t")
