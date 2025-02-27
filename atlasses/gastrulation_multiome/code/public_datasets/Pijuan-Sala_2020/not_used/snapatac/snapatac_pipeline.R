# requires the following packages to be installed / loaded:
# snaptools, macs2, bedtools


scripts <- c("1_snap_preprocess.R",
             "2_snapQC.R",
             "3_dim_reduction.R",
             "4_peak_calling.R",
             "5_cell_by_peak_matrix.R",
             "6_diffacc.R")

run_script <- function(script){
  path <- here::here("public_datasets/Pijuan-Sala_2020/snapatac/", script)
  print(script)
  source(path)
}

lapply(scripts, run_script)

