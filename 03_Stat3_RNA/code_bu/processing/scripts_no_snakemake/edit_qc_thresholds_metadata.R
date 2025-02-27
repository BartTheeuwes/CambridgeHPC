io$metadata <- "/Users/ricard/data/gastrulation_multiome_10x/processed/atac/archR/sample_metadata_after_archR.txt.gz"

# load metadata
sample_metadata <- fread(io$metadata)
sample_metadata[,mean(pass_atacQC==TRUE),by="sample"]

# QC thresholds
opts$min.TSSEnrichment <- 9
opts$min.nFrags <- 4000
opts$max.BlacklistRatio <- 0.05

# edit metadata
sample_metadata %>%
  .[,pass_atacQC:=TSSEnrichment_atac>=opts$min.TSSEnrichment & nFrags_atac>=opts$min.nFrags & BlacklistRatio_atac<=opts$max.BlacklistRatio] %>%
  .[is.na(pass_atacQC),pass_atacQC:=FALSE]
sample_metadata[,(mean(pass_atacQC)),by="sample"]

# save
fwrite(sample_metadata, io$metadata, quote=F, na="NA", sep="\t")
