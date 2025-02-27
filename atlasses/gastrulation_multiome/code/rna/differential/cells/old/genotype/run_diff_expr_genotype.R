here::i_am("rna/differential/cells/genotype/run_diff_expr_genotype.R")

source(here::here("settings.R"))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--metadata',        type="character",     help='Cell metadata file')
p$add_argument('--sce',             type="character",     help='SingleCellExperiment file')
p$add_argument('--group_variable',  type="character",     help='Group variable')
p$add_argument('--wt_class',  type="character",     help='WT class')
p$add_argument('--ko_class',  type="character",     help='KO class')
p$add_argument('--min_cells',       type="integer",       default=50,      help='Minimum number of cells per cell type')
p$add_argument('--outdir',          type="character",     help='Output directory')
p$add_argument('--test_mode',       action="store_true",  help='Test mode? subset data')

args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# args <- list()
# args$metadata <- file.path(io$basedir,"results/rna/mapping/sample_metadata_after_mapping.txt.gz")
# args$sce <- io$rna.sce
# args$group_variable <- "genotype"
# args$ko_class <- "T_KO"
# args$wt_class <- "WT"
# args$min_cells <- 50
# args$outdir <- file.path(io$basedir,"results/rna/differential/genotype")
# args$test_mode <- FALSE
## END TEST ##

#####################
## Define settings ##
#####################

io$script <- here::here("rna/differential/cells/differential.R")
dir.create(args$outdir, showWarnings=FALSE, recursive=TRUE)

opts$samples <- c(
  "E8.5_CRISPR_T_KO",
  "E8.5_CRISPR_T_WT"
)

# opts$rename_celltypes <- c(
#   "Erythroid3" = "Erythroid",
#   "Erythroid2" = "Erythroid",
#   "Erythroid1" = "Erythroid",
#   "Blood_progenitors_1" = "Blood_progenitors",
#   "Blood_progenitors_2" = "Blood_progenitors"
#   # "Intermediate_mesoderm" = "Mixed_mesoderm",
#   # "Paraxial_mesoderm" = "Mixed_mesoderm",
#   # "Nascent_mesoderm" = "Mixed_mesoderm",
#   # "Pharyngeal_mesoderm" = "Mixed_mesoderm"
#   # "Visceral_endoderm" = "ExE_endoderm"
# )

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(args$metadata) %>%
  .[pass_rnaQC==TRUE & doublet_call==FALSE & sample%in%opts$samples]# %>%
  # .[,celltype:=stringr::str_replace_all(celltype,opts$rename_celltypes)]

stopifnot(args$group_variable%in%colnames(sample_metadata))

# Consider cell types with sufficient observations in WT cells
celltypes.to.use <- sample_metadata[genotype==args$wt_class,.(N=.N),by="celltype"] %>% .[N>=args$min_cells,celltype]
sample_metadata <- sample_metadata[celltype%in%celltypes.to.use]

# print stats
celltype_genotype_stats.dt <- table(sample_metadata$celltype,sample_metadata$genotype)
print(celltype_genotype_stats.dt)

#########
## Run ##
#########

# Define cell types to use 
celltypes.to.use <- sample_metadata %>% .[genotype==args$ko_class,.N,by="celltype"] %>% .[N>=args$min_cells,celltype]

if (args$test_mode) {
  print("Test mode activated, running only a few comparisons...")
  celltypes.to.use <- celltypes.to.use %>% head(n=3)
}

stats.dt <- data.table(celltype=as.character(NA), N_WT=as.integer(NA), N_KO=as.integer(NA))

# j <- "Blood_progenitors"
for (i in celltypes.to.use) {
  outfile <- sprintf("%s/%s_%s_vs_%s.txt.gz", args$outdir,i,args$wt_class,args$ko_class); dir.create(dirname(outfile), showWarnings = F)
  if (!file.exists(outfile)) {
    
    # Define LSF command
    if (grepl("BI",Sys.info()['nodename'])) {
      lsf <- ""
    } else if (grepl("pebble|headstone", Sys.info()['nodename'])) {
      lsf <- sprintf("sbatch -n 1 --mem 10G --wrap")
    }
    cmd <- sprintf("%s 'Rscript %s --metadata %s --sce %s --samples %s --celltypes %s --groupA %s --groupB %s --group_variable %s --outfile %s'", 
      lsf, io$script, args$metadata, args$sce, paste(opts$samples,collapse=" "), i, args$wt_class, args$ko_class, args$group_variable, outfile)

  # save stats
    stats.dt <- rbind(stats.dt, data.table(celltype=i, N_WT=celltype_genotype_stats.dt[i,"WT"], N_KO=celltype_genotype_stats.dt[i,"T_KO"]))

    # Run
    print(cmd)
    system(cmd)
  }
}


# Save stats
fwrite(stats.dt[-1], file.path(args$outdir,"diff_stats.txt"), sep="\t", quote=F)

# Completion token
file.create(file.path(args$outdir,"completed.txt"))