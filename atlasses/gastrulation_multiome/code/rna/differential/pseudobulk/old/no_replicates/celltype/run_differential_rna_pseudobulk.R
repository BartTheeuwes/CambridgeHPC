here::i_am("rna/differential/pseudobulk/celltype/run_differential_rna_pseudobulk.R")

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--sce',             type="character",     help='SingleCellExperiment file')
p$add_argument('--outfile',          type="character",     help='Output directory')

args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# args <- list()
# args$sce <- file.path(io$basedir,"results/rna/pseudobulk/celltype.mapped/SingleCellExperiment_pseudobulk.rds")
# args$outfile <- file.path(io$basedir,"results/rna/differential/pseudobulk/celltype/diff_expr_celltype_pseudobulk.txt.gz")
## END TEST ##

#####################
## Define settings ##
#####################

##############################
## Load pseudobulk RNA data ##
##############################

rna.sce <- readRDS(args$sce)

#############################
## Differential expression ##
#############################

celltypes <- colnames(rna.sce)

# for (i in 1:length(celltypes)) {
#   for (j in i:length(celltypes)) {
#     if (i!=j) {
      
#       tmp <- data.table(
#         gene = rownames(rna.sce),
#         expr_A = logcounts(rna.sce[,celltypes[[i]]])[,1] %>% round(2),
#         expr_B = logcounts(rna.sce[,celltypes[[j]]])[,1] %>% round(2)
#       ) %>% .[,diff:=round(expr_B-expr_A,2)] %>% sort.abs("diff") 

#       # save      
#       outfile <- file.path(args$outdir,sprintf("%s_vs_%s.txt.gz",celltypes[[i]],celltypes[[j]]))
#       fwrite(tmp, outfile, sep="\t")
#     }
#   }
# }

diff.dt <- celltypes %>% map(function(i) {
  celltypes %>% map(function(j) {
    if (i!=j) {
      tmp <- data.table(
        gene = rownames(rna.sce),
        expr_A = logcounts(rna.sce[,i])[,1] %>% round(2),
        expr_B = logcounts(rna.sce[,j])[,1] %>% round(2)
      ) %>% .[,diff:=round(expr_B-expr_A,2)] %>% sort.abs("diff") %>%
        .[,c("celltypeA","celltypeB"):=list(i,j)]
    }
  }) %>% rbindlist
}) %>% rbindlist

# Save
fwrite(diff.dt, file.path(args$outfile), sep="\t", na="NA", quote=F)

##########
## TEST ##
##########

# tmp <- data.table(
#   gene = rownames(rna.sce),
#   expr_A = logcounts(rna.sce[,"Endothelium"])[,1] %>% round(2),
#   expr_B = logcounts(rna.sce[,"Epiblast"])[,1] %>% round(2)
# ) %>% .[,diff:=round(expr_B-expr_A,2)] %>% sort.abs("diff") 
