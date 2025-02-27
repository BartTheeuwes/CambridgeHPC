#####################
## Define settings ##
#####################

source(here::here("settings.R"))
io$script <- here::here("rna_atac/gene_regulatory_networks/pseudobulk/per_celltype/gene_regulatory_network_per_celltype.R")
io$outdir <- file.path(io$basedir,"results/rna_atac/gene_regulatory_networks/pseudobulk/per_celltype"); dir.create(io$outdir, showWarnings = F)
io$tmpdir <- file.path(io$outdir,"tmp"); dir.create(io$tmpdir, showWarnings = F)

# genomic distance
# opts$distance <- 1e5

# celltypes to plot
# opts$celltypes <- c("Gut")

# LSF params
opts$memory <- 5000

#########
## Run ##
#########

for (i in opts$celltypes) {
  
  outdir <- sprintf("%s/%s",io$outdir,paste(i,collapse="_")); dir.create(outdir, showWarnings = F)
  # Define LSF command
  if (grepl("BI2404M",Sys.info()['nodename'])) {
    lsf <- ""
  } else if (grepl("babraham",Sys.info()['nodename'])) {
    lsf <- sprintf("bsub -M %s -n 1 -o %s/gene_regulatory_network_pseudobulk_%s.txt", opts$memory, io$tmpdir, paste(i,collapse=" "))
  }
  cmd <- sprintf("%s Rscript %s --celltype %s --outdir %s", lsf, io$script, i, outdir)
  
  # Run
  print(cmd)
  system(cmd)
}
