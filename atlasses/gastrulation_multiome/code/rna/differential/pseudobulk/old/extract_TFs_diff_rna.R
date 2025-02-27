here::i_am("rna/differential/pseudobulk/with_replicates/extract_TFs_diff_rna_metacells.R")

# Load default settings
source(here::here("settings.R"))
# source(here::here("utils.R"))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--TFs',             type="character",     help='Cell metadata file')
p$add_argument('--diff_results_dir',   type="character",     help='File')
p$add_argument('--outfile',             type="character",     help='File')
args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# io$basedir <- file.path(io$basedir,"test")
# args <- list()
# args$diff_results_dir <- file.path(io$basedir,"results/rna/differential/metacells/celltype")
# args$outfile <- file.path(io$basedir,"results/rna/differential/metacells/celltype/TFs/diff_expr_celltype_tfs_metacells.txt.gz")
# args$TFs <- "/Users/argelagr/data/mm10_regulation/TFs/TFs.txt"
## END TEST ##

dir.create(dirname(args$outfile), showWarnings = F)

##############
## Load TFs ##
##############

# TFs <- fread(args$TFs)[["gene"]]
TFs <- fread(args$TFs)[[1]] %>% str_to_title

################################################
## Load differential expression and fetch TFs ##
################################################

diff_tf.dt <- opts$celltypes %>% map(function(i) {
  opts$celltypes %>% map(function(j) {
    if (i!=j) {
      file <- file.path(args$diff_results_dir,sprintf("%s_vs_%s.txt.gz",i,j))
      if (file.exists(file)) {
        fread(file) %>% 
          setnames(c("gene", "logFC", "padj_fdr")) %>%
          .[gene%in%TFs] %>% .[,c("celltypeA","celltypeB"):=list(i,j)]
      }
    }
  }) %>% rbindlist
}) %>% rbindlist %>% .[,gene:=toupper(gene)]

print(sprintf("Number of TFs in the differential expression results: %s",length(unique(diff_tf.dt$gene))))

##########
## Save ##
##########

fwrite(diff_tf.dt, args$outfile, sep="\t", quote=F, na="NA")

