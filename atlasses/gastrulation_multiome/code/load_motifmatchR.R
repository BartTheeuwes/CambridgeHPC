stopifnot(!is.null(opts$motif_annotation))

# io$motifmatcher.se <- sprintf("%s/Annotations/%s-Scores.rds",io$archR.directory,opts$motif_annotation)
io$motifmatcher.se <- sprintf("%s/Annotations/%s-Scores.rds",io$archR.directory,opts$motif_annotation)

# Load
motifmatcher.se <- readRDS(io$motifmatcher.se)