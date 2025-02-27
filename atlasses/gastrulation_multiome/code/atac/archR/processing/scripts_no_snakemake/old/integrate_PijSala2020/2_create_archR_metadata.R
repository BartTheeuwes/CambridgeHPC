########################
## Load ArchR project ##
########################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/atac/archR/load_integrated_archR_project.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/atac/archR/load_integrated_archR_project.R")
} else {
  stop("Computer not recognised")
}

io$multiome.metadata <- io$metadata
io$pijuansala.metadata <- paste0(io$pijuansala.basedir,"/cell_metadata.csv.gz")
io$metadata.out <- paste0(io$archR.directory,"/sample_metadata.txt.gz")

############################
## Load multiome metadata ##
############################

metadata.multiome <- fread(io$multiome.metadata) %>%
  .[,c("cell", "sample", "stage", "barcode", "archR_cell", "nFeature_RNA", "nCount_RNA", "mtFraction_RNA", "pass_rnaQC", 
       "celltype.mapped", "celltype.score", "closest.cell", 
       "cxds_score", "cxds_call", "bcds_score", "bcds_call", "hybrid_score", "hybrid_call", "doublet_call")] %>%
       # "TSSEnrichment_atac", "ReadsInTSS_atac", "PromoterRatio_atac", "NucleosomeRatio_atac", "nFrags_atac", "BlacklistRatio_atac", "pass_atacQC")
  .[,dataset:="Multiome"] %>%
  setnames("celltype.mapped","celltype")

##################################
## Load PijuanSala metadata ##
##################################

metadata.pijuansala <- fread(io$pijuansala.metadata) %>%
  .[,c("cell", "nuclei_type", "celltype")] %>%
  setnames("cell","barcode") %>%
  .[,stage:="E8.25"] %>%
  .[,sample:="E8.25_PijuanSala"] %>%
  .[,dataset:="E8.25_PijuanSala"] %>%
  .[,archR_cell:=sprintf("%s#%s",sample,barcode)] %>%
  .[,cell:=NA]

######################
## Load ArchR stats ##
######################
  
# fetch archR's metadata (first time)
# note that QC is done later in the QC/qc.R script
archR_metadata <- getCellColData(ArchRProject) %>%
  as.data.table(keep.rownames = T) %>% setnames("rn","archR_cell") %>%
  .[,c("archR_cell", "Sample", "TSSEnrichment", "ReadsInTSS", "PromoterRatio", "NucleosomeRatio", "nFrags",  "BlacklistRatio")]# %>%
  # setnames("Sample","sample") %>%
  # .[,cell:=stringr::str_replace_all(archR_cell,"#","_")] %>%
  # .[,sample:=strsplit(archR_cell,"#") %>% map_chr(1)] %>%
  # .[,barcode:=strsplit(archR_cell,"#") %>% map_chr(2)]

cols.to.rename <- c("TSSEnrichment","ReadsInTSS","PromoterRatio","NucleosomeRatio","nFrags","BlacklistRatio")
idx.cols.to.rename <- which(colnames(archR_metadata)%in%cols.to.rename)
colnames(archR_metadata)[idx.cols.to.rename] <- paste0(colnames(archR_metadata)[idx.cols.to.rename], "_atac")

###########
## Merge ##
###########

metadata <-  plyr::rbind.fill(metadata.multiome, metadata.pijuansala) %>%
  as.data.table %>%
  merge(archR_metadata,by="archR_cell") %>%
  .[,Sample:=NULL]

head(metadata,n=3)
tail(metadata,n=3)
table(metadata$sample)
table(metadata$stage)

#############################
## Update ArchR's metadata ##
#############################

metadata.to.archR <- metadata %>% 
  .[archR_cell%in%rownames(ArchRProject)] %>% setkey(archR_cell) %>% .[rownames(ArchRProject)] %>%
  as.data.frame() %>% tibble::column_to_rownames("archR_cell")

stopifnot(all(metadata.to.archR$TSSEnrichment_atac == getCellColData(ArchRProject, "TSSEnrichment")[[1]]))

for (i in colnames(metadata.to.archR)) {
  ArchRProject <- addCellColData(
    ArchRProject,
    data = metadata.to.archR[[i]], 
    name = i,
    cells = rownames(metadata.to.archR),
    force = TRUE
  )
}

head(getCellColData(ArchRProject))

##########
## Save ##
##########

fwrite(metadata, io$metadata.out, sep="\t", na="NA", quote=F)

saveArchRProject(ArchRProject)
