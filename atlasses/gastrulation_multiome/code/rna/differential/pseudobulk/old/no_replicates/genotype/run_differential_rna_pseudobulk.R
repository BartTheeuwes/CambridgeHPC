here::i_am("rna/differential/pseudobulk/genotype/run_differential_rna_pseudobulk.R")

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--sce',             type="character",     help='SingleCellExperiment file')
p$add_argument('--min_cells',             type="integer",     help='Minimum number of cells per cell type')
p$add_argument('--pseudobulk_stats',             type="character",     help='Pseudobulk stats')
p$add_argument('--outfile',          type="character",     help='Output directory')

args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# args <- list()
# args$sce <- file.path(io$basedir,"results/rna/pseudobulk/celltype_genotype/SingleCellExperiment_pseudobulk.rds")
# args$pseudobulk_stats <- file.path(io$basedir,"results/rna/pseudobulk/celltype_genotype/stats.txt")
# args$min_cells <- 30
# args$outfile <- file.path(io$basedir,"results/rna/differential/pseudobulk/genotype/diff_expr_genotype_pseudobulk.txt.gz")
## END TEST ##

dir.create(dirname(args$outfile), showWarnings=F, recursive=T)

###############
## Load data ##
###############

# Load SingleCellExperiment
rna.sce <- readRDS(args$sce)

# Load stats
rna_pseudobulk_stats.dt <- fread(args$pseudobulk_stats)

# Sanity checks
stopifnot(colnames(rna.sce)%in%rna_pseudobulk_stats.dt$group)
print(rna_pseudobulk_stats.dt)

#############################
## Differential expression ##
#############################

tmp <- table(strsplit(rna_pseudobulk_stats.dt[N>=args$min_cells,group], split = "-") %>% map_chr(1))
celltypes.to.use <- tmp[tmp==2] %>% names
cat(sprintf("Using the following cell types that have more than %d cells:\n- %s\n",args$min_cells,paste(celltypes.to.use, collapse="\n- ")))

# i <- "NMP"
diff.dt <- celltypes.to.use %>% map(function(i) {
  
    foo <- logcounts(rna.sce[,paste0(i,"-WT")])[,1] %>% round(2)
    bar <- logcounts(rna.sce[,paste0(i,"-T_KO")])[,1] %>% round(2)
    tmp <- data.table(
      gene = names(foo), 
      expr_WT = foo,
      expr_T_KO = bar,
      diff = round(bar-foo,2), 
      celltype = i
    ) %>% sort.abs("diff") %>% return
}) %>% rbindlist

# Save
fwrite(diff.dt, args$outfile, sep="\t")
