here::i_am("rna/differential/cells/celltype/run_diff_expr_celltype.R")

source(here::here("settings.R"))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--metadata',        type="character",     help='Cell metadata file')
p$add_argument('--sce',             type="character",     help='SingleCellExperiment file')
p$add_argument('--group_variable',  type="character",     help='')
p$add_argument('--min_cells',       type="integer",       default=50,      help='Minimum number of cells per cell type')
p$add_argument('--outdir',          type="character",     help='Output directory')
p$add_argument('--ignore_small_celltypes',       action="store_true",  help='Ignore cell types with a small number of cells')
p$add_argument('--test_mode',       action="store_true",  help='Test mode? subset data')

args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# args <- list()
# args$metadata <- file.path(io$basedir,"results/rna/mapping/sample_metadata_after_mapping.txt.gz")
# args$sce <- io$rna.sce
# args$group_variable <- "celltype.mapped"
# args$min_cells <- 50
# args$outdir <- file.path(io$basedir,"results/rna/differential/celltype.mapped")
# args$test_mode <- TRUE
## END TEST ##

#####################
## Define settings ##
#####################

io$script <- here::here("rna/differential/cells/differential.R")
dir.create(args$outdir, showWarnings=FALSE, recursive=TRUE)

opts$samples <- c(
  "E7.5_rep1",
  "E7.5_rep2",
  "E7.75_rep1",
  "E8.0_rep1",
  "E8.0_rep2",
  "E8.5_rep1",
  "E8.5_rep2",
  "E8.75_rep1",
  "E8.75_rep2",
  # "E8.5_CRISPR_T_KO",
  "E8.5_CRISPR_T_WT"
)

opts$celltypes <- c(
  "Epiblast",
  "Primitive_Streak",
  "Caudal_epiblast",
  "PGC",
  "Anterior_Primitive_Streak",
  "Notochord",
  "Def._endoderm",
  "Gut",
  "Nascent_mesoderm",
  "Mixed_mesoderm",
  "Intermediate_mesoderm",
  "Caudal_Mesoderm",
  "Paraxial_mesoderm",
  "Somitic_mesoderm",
  "Pharyngeal_mesoderm",
  "Cardiomyocytes",
  "Allantois",
  "ExE_mesoderm",
  "Mesenchyme",
  "Haematoendothelial_progenitors",
  "Endothelium",
  "Blood_progenitors_1",
  "Blood_progenitors_2",
  "Erythroid1",
  "Erythroid2",
  "Erythroid3",
  "NMP",
  "Rostral_neurectoderm",
  "Caudal_neurectoderm",
  "Neural_crest",
  "Forebrain_Midbrain_Hindbrain",
  "Spinal_cord",
  "Surface_ectoderm",
  "Visceral_endoderm",
  "ExE_endoderm",
  "ExE_ectoderm",
  "Parietal_endoderm"
)

########################
## Load cell metadata ##
########################

sample_metadata <- fread(args$metadata) %>%
  .[pass_rnaQC==TRUE & doublet_call==FALSE & sample%in%opts$samples]

stopifnot(args$group_variable%in%colnames(sample_metadata))

sample_metadata <- sample_metadata %>%
  .[,celltype:=eval(as.name(args$group_variable))] %>%
  .[celltype%in%opts$celltypes]

# subset celltypes with sufficient number of cells
if (args$ignore_small_celltypes) {
  sample_metadata <- sample_metadata %>%
    .[,N:=.N,by=c("celltype")] %>% .[N>args$min_cells] %>% .[,N:=NULL]
}
celltypes.to.use <- unique(sample_metadata$celltype)# %>% head(n=3)

# print stats
celltype.stats <- table(sample_metadata$celltype)
print(celltype.stats)

#########
## Run ##
#########

if (args$test_mode) {
  print("Test mode activated, running only a few comparisons...")
  celltypes.to.use <- celltypes.to.use %>% head(n=3)
}

stats.dt <- data.table(groupA=as.character(NA), groupB=as.character(NA), N_A=as.integer(NA), N_B=as.integer(NA))

for (i in 1:length(celltypes.to.use)) {
  for (j in i:length(celltypes.to.use)) {
    if (i!=j) {
      groupA <- celltypes.to.use[[i]]
      groupB <- celltypes.to.use[[j]]
      # print(sprintf("i=%s (%s), j=%s (%s)",i,groupA,j,groupB))
      
      outfile <- sprintf("%s/%s_vs_%s.txt.gz", args$outdir,groupA,groupB)
      if (!file.exists(outfile)) {

        # Define LSF command
        if (grepl("BI",Sys.info()['nodename'])) {
          lsf <- ""
        } else if (grepl("pebble|headstone", Sys.info()['nodename'])) {
          lsf <- sprintf("sbatch -n 1 --mem 12G --wrap")
        }
        cmd <- sprintf("%s 'Rscript %s --metadata %s --sce %s --samples %s --groupA %s --groupB %s --group_variable %s --outfile %s'", 
          lsf, io$script, args$metadata, args$sce, paste(opts$samples,collapse=" "), groupA, groupB, args$group_variable, outfile)
        # if (isTRUE(opts$test_mode)) cmd <- paste0(cmd, " --test_mode")

        # Run
        print(cmd)
        system(cmd)
        
        # save stats
        stats.dt <- rbind(stats.dt, data.table(groupA=groupA, groupB=groupB, N_A=celltype.stats[[groupA]], N_B=celltype.stats[[groupB]]))
        
      } else {
        print(sprintf("%s already exists...",outfile))
      }
    }
  }
}


# Save stats
fwrite(stats.dt[-1], file.path(args$outdir,"diff_stats.txt"), sep="\t", quote=F)

# Completion token
file.create(file.path(args$outdir,"completed.txt"))