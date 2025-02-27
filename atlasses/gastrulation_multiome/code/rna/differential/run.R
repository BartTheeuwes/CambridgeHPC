#####################
## Define settings ##
#####################

source(here::here("settings.R"))

# io$Rscript <- "/Library/Frameworks/R.framework/Versions/Current/Resources/bin/Rscript"
# io$Rscript <- "Rscript"

io$script <- here::here("rna/differential/differential.R")
io$tmpdir <- file.path(io$basedir,"results_new/rna/differential/tmp"); dir.create(io$tmpdir, showWarnings=F)
io$outdir <- file.path(io$basedir,"results_new/rna/differential/test"); dir.create(io$outdir, showWarnings=F)

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
          lsf <- ""
        } else if (grepl("pebble|headstone", Sys.info()['nodename'])) {
          lsf <- sprintf("sbatch -n 1 --mem 10G --wrap")
        }
        cmd <- sprintf("%s 'Rscript %s --groupA %s --groupB %s --group_label %s --outfile %s'", lsf, io$script, groupA, groupB, opts$group_label, outfile)
        if (isTRUE(opts$test_mode)) cmd <- paste0(cmd, " --test_mode")
        
        # Run
        print(cmd)
        system(cmd)
      }
    }
  }
}


##############################
## Run selected comparisons ##
##############################

# opts$comparisons <- list(
#   c("groupA"="Mixed_mesoderm",          "groupB"="ExE_ectoderm")
#   # c("groupA"="Blood_progenitors_2", "groupB"="Caudal_Mesoderm"),
#   # c("groupA"="Allantois",           "groupB"="Haematoendothelial_progenitors"),
#   # c("groupA"="Blood_progenitors_2", "groupB"="NMP"),
#   # c("groupA"="Erythroid1",          "groupB"="Anterior_Primitive_Streak"),
#   # c("groupA"="Erythroid1",          "groupB"="Haematoendothelial_progenitors")
# )
# 
# for (comparison in opts$comparisons) {
#   groupA <- comparison[["groupA"]]; groupB <- comparison[["groupB"]]
#   for (test in opts$statistical.test) {
#     outfile <- sprintf("%s/%s_vs_%s.txt.gz", io$outdir,groupA,groupB)
#     
#     # Define LSF command
#     if (grepl("ricard",Sys.info()['nodename'])) {
#       lsf <- ""
#     } else if (grepl("ebi",Sys.info()['nodename'])) {
#       lsf <- sprintf("bsub -M 15000 -n 1 -q research-rh74 -o %s/%s_vs_%s.txt", io$tmpdir,groupA,groupB)
#     }
#     cmd <- sprintf("%s Rscript %s --stages %s --groupA %s --groupB %s --test %s --outfile %s", lsf, io$script, paste(opts$stages, collapse=" "), groupA, groupB, test, outfile)
#     if (isTRUE(opts$test_mode)) cmd <- paste0(cmd, " --test_mode")
#     
#     # Run
#     print(cmd)
#     system(cmd)
#   }
# }
# 
