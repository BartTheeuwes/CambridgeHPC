if (!is.null(opts$motif_annotation)) {
  io$motifmatcher_positions.se <- sprintf("%s/Annotations/%s-Positions-In-Peaks.rds",io$archR.directory,opts$motif_annotation)
}

# Load
motifmatcher_positions.se <- readRDS(io$motifmatcher_positions.se)

# Rename motifs
names(motifmatcher_positions.se) <- names(motifmatcher_positions.se) %>% toupper %>% stringr::str_split(.,"_") %>% map_chr(1)
names(motifmatcher_positions.se) <- gsub("TCFAP","TFAP",names(motifmatcher_positions.se))
names(motifmatcher_positions.se) <- gsub("NKX2","NKX2-",names(motifmatcher_positions.se))
names(motifmatcher_positions.se) <- gsub("NKX3","NKX3-",names(motifmatcher_positions.se))
names(motifmatcher_positions.se) <- gsub("NKX6","NKX6-",names(motifmatcher_positions.se))