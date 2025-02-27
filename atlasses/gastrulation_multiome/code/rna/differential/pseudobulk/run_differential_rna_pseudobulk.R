#####################
## Define settings ##
#####################

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

# Options

# I/O
io$outdir <- file.path(io$basedir,"results_new/rna/differential/pseudobulk/test"); dir.create(io$outdir, showWarnings = F)

##############################
## Load pseudobulk RNA data ##
##############################

rna.sce <- readRDS(io$rna.pseudobulk.sce)[,opts$celltypes]

#############################
## Differential expression ##
#############################

for (i in 1:length(opts$celltypes)) {
  for (j in i:length(opts$celltypes)) {
    if (i!=j) {
      
      tmp <- data.table(
        gene = rownames(rna.sce),
        expr_A = logcounts(rna.sce[,opts$celltypes[[i]]])[,1] %>% round(2),
        expr_B = logcounts(rna.sce[,opts$celltypes[[j]]])[,1] %>% round(2)
      ) %>% .[,diff:=round(expr_B-expr_A,2)] %>% sort.abs("diff") 

      # save      
      outfile <- file.path(io$outdir,sprintf("%s_vs_%s.txt.gz",opts$celltypes[[i]],opts$celltypes[[j]]))
      fwrite(tmp, outfile, sep="\t")
    }
  }
}
