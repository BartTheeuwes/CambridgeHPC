here::i_am("atac/archR/differential/metacells/genotype/run_diff_acc_genotype_metacells.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--atac_pseudobulk_file',        type="character",                               help='')
p$add_argument('--matrix',          type="character",  default="PeakMatrix",   help='Matrix to use')
p$add_argument('--outdir',          type="character",                               help='Output directory')
p$add_argument('--test_mode',    action="store_true",             help='Test mode? subset data')

args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# io$basedir <- file.path(io$basedir,"test")
# args <- list()
# args$matrix <- "PeakMatrix" # "GeneScoreMatrix_TSS"
# args$atac_pseudobulk_file <- file.path(io$basedir,sprintf("results/atac/archR/pseudobulk/celltype_genotype/%s/pseudobulk_%s_summarized_experiment.rds",args$matrix,args$matrix))
# args$outdir <- file.path(io$basedir,sprintf("results/atac/archR/differential/pseudobulk/genotype/%s",args$matrix))
# args$test_mode <- TRUE
## END TEST ##

# I/O
dir.create(args$outdir, showWarnings=F, recursive = T)

#################################
## Load pseudobulk ATAC Matrix ##
#################################

print(sprintf("Fetching pseudobulk ATAC matrix '%s'  grouped by '%s'...", args$matrix,args$group_variable))

atac_pseudobulk.se <- readRDS(args$atac_pseudobulk_file)

# Normalise
assay(atac_pseudobulk.se) <- log2(1e6*(sweep(assay(atac_pseudobulk.se),2,colSums(assay(atac_pseudobulk.se),na.rm=T),"/"))+0.5)

#########
## Run ##
#########

tmp <- colnames(atac_pseudobulk.se) %>% strsplit("-") %>% map_chr(1) %>% table
celltypes <- which(tmp==2) %>% names
# celltypes <- opts$celltypes

if (args$test_mode) {
  print("Test mode activated, running only a few comparisons...")
  celltypes <- celltypes %>% head(n=3)
}

for (i in celltypes) {
  foo <- assay(atac_pseudobulk.se[,paste0(i,"-WT")])[,1]
  bar <- assay(atac_pseudobulk.se[,paste0(i,"-T_KO")])[,1]
  
  atac_diff.dt <- data.table(
    idx = names(foo), 
    diff = round(bar-foo,2) 
    # groupA = celltypes[[i]], 
    # groupB = celltypes[[j]]
  ) %>% sort.abs("diff") 
  
  # save      
  outfile <- file.path(args$outdir,sprintf("%s_WT_vs_KO_pseudobulk.txt.gz",i))
  fwrite(atac_diff.dt, outfile, sep="\t")
}

# Completion token
file.create(file.path(args$outdir,"completed.txt"))