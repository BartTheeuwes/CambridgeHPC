here::i_am("rna/differential/celltypes/pseudobulk/run_differential_rna_pseudobulk.R")

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# Options

# I/O
io$outdir <- file.path(io$basedir,"results/rna/differential/celltype.mapped/pseudobulk"); dir.create(io$outdir, showWarnings = F)
io$rna.pseudobulk.sce <- file.path(io$basedir,"results/rna/pseudobulk/celltype.mapped/SingleCellExperiment_pseudobulk.rds")

##############################
## Load pseudobulk RNA data ##
##############################

rna.sce <- readRDS(io$rna.pseudobulk.sce)

celltypes <- colnames(rna.sce)

#############################
## Differential expression ##
#############################

for (i in 1:length(celltypes)) {
  for (j in i:length(celltypes)) {
    if (i!=j) {
      
      tmp <- data.table(
        gene = rownames(rna.sce),
        expr_A = logcounts(rna.sce[,celltypes[[i]]])[,1] %>% round(2),
        expr_B = logcounts(rna.sce[,celltypes[[j]]])[,1] %>% round(2)
      ) %>% .[,diff:=round(expr_B-expr_A,2)] %>% sort.abs("diff") 

      # save      
      outfile <- file.path(io$outdir,sprintf("%s_vs_%s.txt.gz",celltypes[[i]],celltypes[[j]]))
      fwrite(tmp, outfile, sep="\t")
    }
  }
}


##########
## TEST ##
##########

tmp <- data.table(
  gene = rownames(rna.sce),
  expr_A = logcounts(rna.sce[,"Endothelium"])[,1] %>% round(2),
  expr_B = logcounts(rna.sce[,"Epiblast"])[,1] %>% round(2)
) %>% .[,diff:=round(expr_B-expr_A,2)] %>% sort.abs("diff") 
