here::i_am("atac/archR/chromvar_chip/pseudobulk/differential/run_differential_chromvar_chip_pseudobulk.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--motif_annotation',  type="character",              help='Motif annotation') 
p$add_argument('--chromvar_chip',  type="character",              help='') 
p$add_argument('--outdir',  type="character",              help='Motif annotation') 
args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
io$basedir <- file.path(io$basedir,"test")
args$motif_annotation <- "JASPAR"
args$chromvar_chip <- file.path(io$basedir,sprintf("results/atac/archR/chromvar_chip/pseudobulk/chromVAR_chip_%s_archr.rds",args$motif_annotation))
args$outdir <- file.path(io$basedir,sprintf("results/atac/archR/chromvar_chip/pseudobulk/differential/celltype/%s",args$motif_annotation))
## END TEST ##

#####################
## Define settings ##
#####################

# I/O
dir.create(args$outdir, showWarnings=F, recursive=T)

#####################################
## Load pseudobulk chromVAR scores ##
#####################################

atac_chromvar_pseudobulk.se <- readRDS(args$chromvar_chip)
atac_chromvar_pseudobulk.se

######################################
## Differential motif accessibility ##
######################################

celltypes <- colnames(atac_chromvar_pseudobulk.se)

# i <- 1; j <- 2
for (i in 1:length(celltypes)) {
  for (j in i:length(celltypes)) {
    if (i!=j) {
      foo <- assay(atac_chromvar_pseudobulk.se[,celltypes[[j]]],"deviations")[,1]
      bar <- assay(atac_chromvar_pseudobulk.se[,celltypes[[i]]],"deviations")[,1]
      
      diff.dt <- data.table(
        gene = names(foo), 
        diff = round(foo-bar,2) 
        # groupA = celltypes[[i]], 
        # groupB = celltypes[[j]]
      ) %>% sort.abs("diff") 
      
      # save      
      outfile <- file.path(args$outdir,sprintf("chromVAR_%s_vs_%s.txt.gz",celltypes[[i]],celltypes[[j]]))
      fwrite(diff.dt, outfile, sep="\t")
    }
  }
}

# Completion token
file.create(file.path(args$outdir,"completed.txt"))
