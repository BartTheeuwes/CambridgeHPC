########################
## Load ArchR project ##
########################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/atac/archR/load_archR_project.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/atac/archR/load_archR_project.R")
} else {
  stop("Computer not recognised")
}

############################################
## Merge archR metadata with RNA metadata ##
############################################


# fetch archR's metadata (first time)
# note that QC is done later in the QC/qc.R script
archr_metadata <- getCellColData(ArchRProject) %>%
  .[,c("cell", "barcode", "sample", "TSSEnrichment", "ReadsInTSS", "PromoterRatio", "NucleosomeRatio", "nFrags",  "BlacklistRatio")] %>%
  as.data.table(keep.rownames = T) %>%
  setnames("rn","archR_cell") %>%
  .[is.na(cell),cell:=stringr::str_replace_all(archR_cell,"#","_")] %>%
  .[is.na(sample),sample:=strsplit(archR_cell,"#") %>% map_chr(1)] %>%
  .[is.na(barcode),barcode:=strsplit(archR_cell,"#") %>% map_chr(2)]

cols.to.rename <- c("TSSEnrichment","ReadsInTSS","PromoterRatio","NucleosomeRatio","nFrags","BlacklistRatio")
idx.cols.to.rename <- which(colnames(archr_metadata)%in%cols.to.rename)
colnames(archr_metadata)[idx.cols.to.rename] <- paste0(colnames(archr_metadata)[idx.cols.to.rename], "_atac")

# fetch pre-computed archR's metadata
# io$archr.metadata <- paste0(io$basedir,"/processed/atac/archR/sample_metadata_after_archR.txt.gz")
# archr_metadata <- fread(io$archr.metadata)
# cols.to.rename <- c("TSSEnrichment","ReadsInTSS","PromoterRatio","NucleosomeRatio","nFrags","BlacklistRatio")
# idx.cols.to.rename <- which(colnames(archr_metadata)%in%cols.to.rename)
# colnames(archr_metadata)[idx.cols.to.rename] <- paste0(colnames(archr_metadata)[idx.cols.to.rename], "_atac")

# fetch RNA metadata
io$rna.metadata <- paste0(io$basedir,"/results/rna/celltype_denoising/sample_metadata_after_celltype_denoising.txt.gz")
rna_metadata <- fread(io$rna.metadata)
rna_metadata[,log_nCount_RNA:=NULL]
rna_metadata %>% setnames("pass_QC","pass_rnaQC")
rna_metadata %>% setnames("percent.mt","mtFraction_RNA")

# remove overlapping columns (rna_metadata is more updated)
overlaping.columns <- intersect(colnames(rna_metadata),colnames(archr_metadata))
overlaping.columns <- overlaping.columns[!overlaping.columns%in%c("sample","cell","barcode")]
archr_metadata <- archr_metadata[,which(!colnames(archr_metadata)%in%overlaping.columns),with=F]

# merge
foo <- rna_metadata %>% merge(archr_metadata, by=c("cell","sample","barcode"), all=TRUE)
# foo[,archR_cell:=sprintf("%s#%s",sample,barcode)]

#############################
## Update ArchR's metadata ##
#############################

bar <- foo %>% 
  .[archR_cell%in%rownames(ArchRProject)] %>% setkey(archR_cell) %>% .[rownames(ArchRProject)] %>%
  as.data.frame() %>% tibble::column_to_rownames("archR_cell")

all(bar$TSSEnrichment_atac == getCellColData(ArchRProject, "TSSEnrichment")[[1]])

for (i in colnames(bar)) {
  ArchRProject <- addCellColData(
    ArchRProject,
    data = bar[[i]], 
    name = i,
    cells = rownames(bar),
    force = TRUE
  )
}

head(getCellColData(ArchRProject))

##########
## Save ##
##########

io$metadata.out <- paste0(io$basedir,"/processed/atac/archR/sample_metadata_after_archR.txt.gz")
fwrite(foo, io$metadata.out, sep="\t", na="NA", quote=F)

saveArchRProject(ArchRProject)
