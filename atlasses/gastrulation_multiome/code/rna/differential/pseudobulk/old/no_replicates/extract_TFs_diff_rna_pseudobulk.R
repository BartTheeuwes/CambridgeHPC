here::i_am("rna/differential/pseudobulk/no_replicates/extract_TFs_diff_rna_pseudobulk.R")

# Load default settings
source(here::here("settings.R"))
# source(here::here("utils.R"))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--TFs',        			type="character",     help='Cell metadata file')
p$add_argument('--diff_results_file',   type="character",     help='File')
p$add_argument('--outfile',             type="character",     help='File')

args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# args <- list()
# args$TFs <- file.path(io$basedir,"processed/atac/archR/Annotations/CISBP_TFs.txt.gz")
# args$diff_results_file <- file.path(io$basedir,"results/rna/differential/pseudobulk/celltype/diff_expr_celltype_pseudobulk.txt.gz")
# args$outfile <- file.path(io$basedir,"results/rna/differential/pseudobulk/celltype/TFs/diff_expr_celltype_tfs_pseudobulk.txt.gz")
## END TEST ##

dir.create(dirname(args$outfile), showWarnings = F)

##############
## Load TFs ##
##############

TFs <- fread(args$TFs)[["gene"]]

################################################
## Load differential expression and fetch TFs ##
################################################

# celltypes <- opts$celltypes

# diff_tf.dt <- celltypes %>% map(function(i) {
#   celltypes %>% map(function(j) {
#     if (i!=j) {
#       file <- file.path(args$indir,sprintf("%s_vs_%s.txt.gz",i,j))
#       if (file.exists(file)) {
#         fread(file) %>% .[,gene:=toupper(gene)] %>% .[gene%in%TFs] %>% .[,c("celltypeA","celltypeB"):=list(i,j)]
#       }
#     }
#   }) %>% rbindlist
# }) %>% rbindlist

diff_tf.dt <- fread(args$diff_results_file) %>% 
	.[,gene:=toupper(gene)] %>% .[gene%in%TFs]

print(sprintf("Number of TFs: %s",length(TFs)))
print(sprintf("Number of TFs in the differential expression results: %s",length(unique(diff_tf.dt$gene))))

##########
## Save ##
##########

fwrite(diff_tf.dt, args$outfile, sep="\t", quote=F, na="NA")
