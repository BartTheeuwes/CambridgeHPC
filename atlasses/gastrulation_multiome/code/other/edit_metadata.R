a <- fread("/Users/ricard/data/gastrulation_multiome_10x/results/rna/mapping/sample_metadata_after_mapping.txt.gz") %>%
  .[,c("cell","celltype.mapped","celltype.score", "closest.cell")]
b <- fread("/Users/ricard/data/gastrulation_multiome_10x/results/rna/doublets/sample_metadata_after_doublets.txt.gz") %>%
  .[,c("celltype.mapped","celltype.score", "closest.cell"):=NULL]

ab <- merge(a,b,by="cell")

fwrite(ab, "/Users/ricard/data/gastrulation_multiome_10x/results/rna/doublets/sample_metadata_after_doublets.txt.gz", sep="\t", quote=F, na="NA")
