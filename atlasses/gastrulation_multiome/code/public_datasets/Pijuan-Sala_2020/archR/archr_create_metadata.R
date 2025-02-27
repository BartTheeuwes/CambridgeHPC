###################
## Save metadata ##
###################

# c("Sample", "TSSEnrichment", "ReadsInTSS", "ReadsInPromoter", 
#   "ReadsInBlacklist", "PromoterRatio", "PassQC", "NucleosomeRatio", 
#   "nMultiFrags", "nMonoFrags", "nFrags", "nDiFrags", "BlacklistRatio", 
#   "Clusters")

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  addTileMat = FALSE,
  addGeneScoreMat = FALSE,
  excludeChr = c("chrM", "chrY"),
  # QC metrics
  removeFilteredCells = FALSE,
  filterFrags = 1000,          # The minimum number of fragments per cell
  filterTSS = 4,               # The minimum TSS enrichment score per cell
  minFrags = 500,              # The minimum fragments per cell to be filtered before any QC calculations (such as TSS Enrichment Score).
  maxFrags = 2e+05             # The maximum fragments per cell to be filtered before any QC calculations (such as TSS Enrichment Score).
)

ArchRProject <- ArchRProject(
  ArrowFiles,  
  outputDirectory = "/Users/ricard/test/archR", 
  copyArrows = FALSE
)

dt <- getCellColData(ArchRProject, select = c("Sample","nFrags","TSSEnrichment","NucleosomeRatio","PassQC")) %>%
  as.data.table(keep.rownames = T) %>% setnames("rn","cell") %>%
  .[,cell:=stringr::str_replace_all(cell,"#","_") %>% stringr::str_replace_all(.,"-1","")] %>%
  .[,Sample:=NULL] %>%
  setnames("PassQC","pass_atacQC") %>%
  .[,pass_atacQC:=ifelse(pass_atacQC==1,TRUE,FALSE)]

# Sanity checks
# foo <- ArchRProject$cellNames %>% stringr::str_replace_all("#","_") %>% stringr::str_replace_all("-1","")
mean(dt$cell %in% sample_metadata$cell)
# mean(sample_metadata$cell %in% foo)
# sample_metadata[!cell%in%foo] %>% View
dt[!cell%in%sample_metadata$cell] %>% View

# Update sample metadata 
# sample_metadata <- fread(io$metadata) %>%
#   merge(dt, by=c("cell"), all.x=TRUE)
# sample_metadata[is.na(pass_atacQC),pass_atacQC:=FALSE]
# fwrite(sample_metadata, io$metadata, sep="\t", na="NA", quote=F)