#####################
## Define settings ##
#####################

source(here::here("settings.R"))

# io$Rscript <- "/Library/Frameworks/R.framework/Versions/Current/Resources/bin/Rscript"
# io$Rscript <- "Rscript"

io$script <- here::here("rna/differential/TFs/run_differential_rna_TFs.R")

# io$tmpdir <- file.path(io$basedir,"results_new/rna/differential/TFs/tmp"); dir.create(io$tmpdir, showWarnings=F)
io$outdir <- file.path(io$basedir,"results_new/rna/differential/TFs"); dir.create(io$outdir, showWarnings=F)

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(io$metadata) %>% 
  .[pass_rnaQC==TRUE & doublet_call==FALSE & !is.na(celltype.mapped_mnn)]

#############
## Options ##
#############

# Testing mode
opts$test_mode <- FALSE

# Define cell types
# opts$groups <- c("Epiblast", "Primitive_Streak", "Gut")
opts$groups <- names(which(table(sample_metadata[stage%in%opts$stages,celltype.mapped_mnn])>=50))

opts$group_label <- "celltype.mapped_mnn"

###################################
## Run all pair-wise comparisons ##
###################################

for (i in 1:length(opts$groups)) {
  groupA <- opts$groups[[i]]
  for (j in i:length(opts$groups)) {
    if (i!=j) {
      groupB <- opts$groups[[j]]
      outfile <- sprintf("%s/%s_vs_%s.txt.gz", io$outdir,groupA,groupB)
      
      if (!file.exists(outfile)) {
        
        # Define LSF command
        if (grepl("BI",Sys.info()['nodename'])) {
          cmd <- sprintf("Rscript %s --groupA %s --groupB %s --group_label %s --outfile %s", io$script, groupA, groupB, opts$group_label, outfile)
        } else if (grepl("pebble|headstone", Sys.info()['nodename'])) {
          cmd <- sprintf("sbatch -n 1 --mem 10G --wrap 'Rscript %s --groupA %s --groupB %s --group_label %s --outfile %s'", io$script, groupA, groupB, opts$group_label, outfile)
        }
        if (isTRUE(opts$test_mode)) cmd <- paste0(cmd, " --test_mode")
        
        # Run
        print(cmd)
        system(cmd)
      }
    }
  }
}

