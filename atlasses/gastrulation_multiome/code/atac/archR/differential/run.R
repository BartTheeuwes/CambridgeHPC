#####################
## Define settings ##
#####################

# Load default settings
source(here::here("settings.R"))

# I/O
io$script <- here::here("atac/archR/differential/archr_differential_accessibility.R")
io$tmpdir <- file.path(io$basedir,"results_new/atac/archR/differential/tmp"); dir.create(io$tmpdir, showWarnings=F)

# Statistical test
opts$statistical.test <- "wilcoxon"

# Celltype label
opts$celltype_label <- "celltype.mapped_mnn"

# Define matrix
opts$matrix <- "GeneScoreMatrix_TSS" # "PeakMatrix"

opts$ignore_small_celltypes <- TRUE
opts$min.cells <- 100

# Cell types
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

sample_metadata <- fread(io$metadata) %>%
  .[pass_atacQC==TRUE & doublet_call==FALSE]

stopifnot(opts$celltype_label%in%colnames(sample_metadata))

sample_metadata <- sample_metadata %>%
  .[,celltype:=eval(as.name(opts$celltype_label))] %>%
  .[celltype%in%opts$celltypes]

# subset celltypes with sufficient number of cells
if (opts$ignore_small_celltypes) {
  sample_metadata <- sample_metadata %>%
    .[,N:=.N,by=c("celltype")] %>% .[N>opts$min.cells] %>% .[,N:=NULL]
}
opts$celltypes <- unique(sample_metadata$celltype)# %>% head(n=3)

table(sample_metadata$celltype)

#########
## Run ##
#########

io$outdir <- file.path(io$basedir,sprintf("results_new/atac/archR/differential/%s",opts$matrix)); dir.create(io$outdir, showWarnings=F)

for (i in 1:length(opts$celltypes)) {
  for (j in i:length(opts$celltypes)) {
    if (i!=j) {
      groupA <- opts$celltypes[[i]]
      groupB <- opts$celltypes[[j]]
      # print(sprintf("i=%s (%s), j=%s (%s)",i,groupA,j,groupB))
      
      outfile <- sprintf("%s/%s_%s_vs_%s.txt.gz", io$outdir,opts$matrix,groupA,groupB)
      if (!file.exists(outfile)) {

        # Define LSF command
        if (grepl("BI",Sys.info()['nodename'])) {
          lsf <- ""
        } else if (grepl("ebi",Sys.info()['nodename'])) {
          lsf <- sprintf("bsub -M 7000 -n 1 -o %s/%s_%s_vs_%s.txt", io$tmpdir,opts$matrix,groupA,groupB)
        } else if (grepl("pebble|headstone", Sys.info()['nodename'])) {
          lsf <- sprintf("sbatch -n 1 --mem 7G --wrap")
        }
        cmd <- sprintf("%s 'Rscript %s --groupA %s --groupB %s --matrix %s --celltype_label %s --test %s --outfile %s'", lsf, io$script, groupA, groupB, opts$matrix, opts$celltype_label, opts$statistical.test, outfile)
        # if (isTRUE(opts$test_mode)) cmd <- paste0(cmd, " --test_mode")


        # Run
        print(cmd)
        system(cmd)
      }
    }
  }
}

