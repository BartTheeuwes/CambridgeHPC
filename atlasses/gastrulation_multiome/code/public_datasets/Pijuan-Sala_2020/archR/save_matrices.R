if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/public_datasets/Pijuan-Sala_2020/archR/load_archR_project.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/public_datasets/Pijuan-Sala_2020/archR/load_archR_project.R")
} else {
  stop("Computer not recognised")
}

io$outfile <- sprintf("%s/data/processed/ArchR/Matrices/GeneScoreMatrix.rds",io$basedir)
foo <- getMatrixFromProject(ArchRProject, "GeneScoreMatrix")
saveRDS(foo, io$outfile)

io$outfile <- sprintf("%s/data/processed/ArchR/Matrices/TileMatrix.rds",io$basedir)
foo <- getMatrixFromProject(ArchRProject, "TileMatrix")
saveRDS(foo, io$outfile)