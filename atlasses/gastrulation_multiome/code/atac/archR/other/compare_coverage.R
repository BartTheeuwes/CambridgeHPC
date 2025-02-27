sample_metadata.old <- fread("/Users/argelagr/data/gastrulation_multiome_10x/results/atac/archR/qc/sample_metadata_after_qc.txt.gz")
sample_metadata.new <- fread("/Users/argelagr/data/gastrulation_multiome_10x/test/results/atac/archR/qc/sample_metadata_after_qc.txt.gz")

foo <- merge(
  sample_metadata.old[pass_atacQC==TRUE,c("cell","nFrags_atac","sample")],
  sample_metadata.new[pass_atacQC==TRUE,c("cell","nFrags_atac","sample")],
  suffixes=c(".old",".new"),
  by=c("cell","sample")
) 

to.plot <- foo %>% melt(id.vars=c("cell","sample"))
ggscatter(foo[nFrags_atac.new<=10000][sample.int(.N,size=1000)], x="nFrags_atac.old", y="nFrags_atac.new", size=0.5) +
  geom_abline(intercept=0, slope=1)
