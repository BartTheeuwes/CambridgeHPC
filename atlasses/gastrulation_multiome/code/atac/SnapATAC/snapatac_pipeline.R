library(SnapATAC)
library(purrr)
library(data.table)



# requires the following packages to be installed / loaded:
# snaptools, macs2, bedtools

source(here::here("settings.R"))
source(here::here("atac/SnapATAC/snapatac_settings.R"))

logfile <- gsub(".rds", paste0("_", Sys.Date() ,"_logfile.txt"), snapio$rds_file)
logfile
dir.create(dirname(logfile), recursive = TRUE)

scripts <- c(
  #"1_fragments_to_snap.R",
  #"2_load_snap_and_QC.R",
  #"3_dim_reduction.R",
  "4_peak_calling.R",
  "5_cell_by_peak_matrix.R",
  "6_diffacc.R",
  "7_ChromVar.R",
  "8_plot_motifs.R" 
)

run_script <- function(script, logfile){
  path <- here::here("atac/SnapATAC", script)
  
  print(script)
  
  log <- paste0(Sys.time(), ": running ", script)
  cat(log, 
      file = logfile, 
      sep = "\n",
      append = TRUE)
  
  source(path)
  
  log <- paste0(Sys.time(), ": ", script,  " script completed")
  cat(log, 
      file = logfile, 
      sep = "\n",
      append = TRUE)
}

lapply(scripts, run_script, logfile = logfile)


