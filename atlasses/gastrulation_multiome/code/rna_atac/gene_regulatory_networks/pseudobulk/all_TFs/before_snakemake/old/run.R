#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
  io$script <- "/Users/ricard/gastrulation_multiome_10x/rna_atac/gene_regulatory_networks/pseudobulk/gene_regulatory_network_test.R"
} else {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
  io$script <- "/homes/ricard/gastrulation_multiome_10x/rna_atac/gene_regulatory_networks/pseudobulk/gene_regulatory_network_test.R"
}
io$outdir <- paste0(io$basedir,"/results/rna_atac/gene_regulatory_networks/pseudobulk"); dir.create(io$outdir, showWarnings = F)
io$tmpdir <- paste0(io$outdir,"/tmp"); dir.create(io$tmpdir, showWarnings = F)

# Cell types
# opts$celltypes <- ...

# genomic distance
opts$distance <- 1e5

# TFs to plot
# opts$genes <- c("Foxa2","Gata1")
opts$genes <- fread(paste0(io$basedir,"/results/rna_atac/gene_regulatory_networks/pseudobulk/TFs.txt"))[[1]]# %>% head(n=3)

# LSF params
opts$memory <- 5000

###################
## Run stringent ##
###################

for (i in opts$genes) {
  
  outdir <- sprintf("%s/%s",io$outdir,i); dir.create(outdir, showWarnings = F)
  # Define LSF command
  if (grepl("ricard",Sys.info()['nodename'])) {
    lsf <- ""
  } else if (grepl("ebi",Sys.info()['nodename'])) {
    lsf <- sprintf("bsub -M %s -n 1 -o %s/gene_regulatory_network_pseudobulk_gene%s_stringent.txt", opts$memory, io$tmpdir, i)
  }
  cmd <- sprintf("%s Rscript %s --genes %s --distance %d --stringent --outdir %s", lsf, io$script, i, opts$distance, outdir)
  
  # Run
  print(cmd)
  system(cmd)
}

#################
## Run lenient ##
#################

for (i in opts$genes) {
  
  outdir <- sprintf("%s/%s",io$outdir,i); dir.create(outdir, showWarnings = F)
  # Define LSF command
  if (grepl("ricard",Sys.info()['nodename'])) {
    lsf <- ""
  } else if (grepl("ebi",Sys.info()['nodename'])) {
    lsf <- sprintf("bsub -M %s -n 1 -o %s/gene_regulatory_network_pseudobulk_gene%s_lenient.txt", opts$memory, io$tmpdir, i)
  }
  cmd <- sprintf("%s Rscript %s --genes %s --distance %d --outdir %s", lsf, io$script, i, opts$distance, outdir)
  
  # Run
  print(cmd)
  system(cmd)
}