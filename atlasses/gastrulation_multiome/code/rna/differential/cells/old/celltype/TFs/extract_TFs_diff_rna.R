# here::i_am("rna/differential/celltypes/pseudobulk/run_differential_rna_pseudobulk.R")

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# Options

# I/O
io$indir <- file.path(io$basedir,"results/rna/differential/celltype.mapped")
io$outdir <- file.path(io$basedir,"results/rna/differential/celltype.mapped/TFs"); dir.create(io$outdir, showWarnings = F)
io$TFs <- file.path(io$basedir,"processed/atac/archR/Annotations/CISBP_TFs.txt.gz")

##############
## Load TFs ##
##############

TFs <- fread(io$TFs)[["gene"]]

################################################
## Load differential expression and fetch TFs ##
################################################

celltypes <- opts$celltypes

diff_tf.dt <- celltypes %>% map(function(i) {
  celltypes %>% map(function(j) {
    if (i!=j) {
      file <- file.path(io$indir,sprintf("%s_vs_%s.txt.gz",i,j))
      if (file.exists(file)) {
        fread(file) %>% 
          .[,c("gene","logFC","padj_fdr")] %>%
          .[,gene:=toupper(gene)] %>% .[gene%in%TFs] %>% 
          .[,c("celltypeA","celltypeB"):=list(i,j)]
      }
    }
  }) %>% rbindlist
}) %>% rbindlist

##########
## Save ##
##########

fwrite(diff_tf.dt, file.path(io$outdir,"diff_expr_TFs.txt.gz"), sep="\t", quote=F,na="NA")
