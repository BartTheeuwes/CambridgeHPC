here::i_am("atac/archR/processing/2_create_archR_metadata.R")

source(here::here("settings.R"))

suppressPackageStartupMessages(library(ArchR))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--archr_directory',    type="character",    help='ArchR directory')
p$add_argument('--metadata',    type="character",    help='metadata file')
p$add_argument('--outfile',     type="character",    help='Output file')

args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# args$archr_directory <- file.path(io$basedir,"processed/atac/archR")
# args$metadata <- file.path(io$basedir,"results/rna/mapping/sample_metadata_after_mapping.txt.gz")
# args$outfile <- file.path(io$basedir,"processed/atac/archR/sample_metadata_after_archR.txt.gz")
## END TEST ##

###################
## Load metadata ##
###################

# Note that this metadata file is derived from the RNA pipeline
sample_metadata <- fread(args$metadata)

########################
## Load ArchR project ##
########################

# source(here::here("atac/archR/load_archR_project.R"))

setwd(args$archr_directory)

addArchRGenome("mm10")
addArchRThreads(threads = 1)

ArchRProject <- loadArchRProject(args$archr_directory)#[sample_metadata$cell]

# print stats
print(sprintf("Number of cells that are derived from the ATAC pipeline (before QC): %s",length(rownames(ArchRProject))))
print(sprintf("Number of cells that are derived from the RNA pipeline (before QC): %s",nrow(sample_metadata)))
print(sprintf("Number of cells that are derived from the RNA pipeline but are not present in the ArchR project: %s",sum(!sample_metadata$cell%in%rownames(ArchRProject))))
print(table(sample_metadata[!cell%in%rownames(ArchRProject),sample]))
print(rownames(ArchRProject@cellColData)[!rownames(ArchRProject@cellColData)%in%sample_metadata$cell] %>% strsplit("#") %>% map_chr(1) %>% table)
print(sprintf("Number of cells that are derived from the ArchR project but are not present in the RNA pipeline: %s",sum(!rownames(ArchRProject@cellColData)%in%sample_metadata$cell)))
print(rownames(ArchRProject@cellColData)[!rownames(ArchRProject@cellColData)%in%sample_metadata$cell] %>% strsplit("#") %>% map_chr(1) %>% table)

nrow(sample_metadata[sample=="rv_eo_deg_day3_5_control"])
nrow(sample_metadata[pass_rnaQC==T & sample=="rv_eo_deg_day3_5_control"])

######################
## Load ArchR stats ##
######################
  
# fetch archR's metadata
archR_metadata <- getCellColData(ArchRProject) %>%
  as.data.table(keep.rownames = T) %>% setnames("rn","cell") %>%
  .[,c("cell", "TSSEnrichment", "ReadsInTSS", "PromoterRatio", "NucleosomeRatio", "nFrags",  "BlacklistRatio")]

cols.to.rename <- c("TSSEnrichment","ReadsInTSS","PromoterRatio","NucleosomeRatio","nFrags","BlacklistRatio")
idx.cols.to.rename <- which(colnames(archR_metadata)%in%cols.to.rename)
colnames(archR_metadata)[idx.cols.to.rename] <- paste0(colnames(archR_metadata)[idx.cols.to.rename], "_atac")

###########
## Merge ##
###########

sample_metadata_tosave <- sample_metadata %>% 
  merge(archR_metadata,by="cell", all=TRUE) 

# Fill missing entries
sample_metadata_tosave %>%
  .[is.na(sample),sample:=strsplit(cell,"#") %>% map_chr(1)] %>%
  .[is.na(barcode),barcode:=strsplit(cell,"#") %>% map_chr(2)]

# round
sample_metadata_tosave[,c("TSSEnrichment_atac","NucleosomeRatio_atac","PromoterRatio_atac","BlacklistRatio_atac"):=list(round(TSSEnrichment_atac,2),round(NucleosomeRatio_atac,2),round(PromoterRatio_atac,2),round(BlacklistRatio_atac,2))]
sample_metadata_tosave[,c("ribosomal_percent_RNA","mitochondrial_percent_RNA"):=list(round(ribosomal_percent_RNA,2),round(mitochondrial_percent_RNA,2))]

# sanity checks
print("Number of cells that pass atac QC per sample:")
table(sample_metadata$sample)
stopifnot(all(!is.na(sample_metadata_tosave$sample)))
stopifnot(all(!is.na(sample_metadata_tosave$barcode)))

#############################
## Update ArchR's metadata ##
#############################

metadata.to.archR <- sample_metadata_tosave %>% 
  .[cell%in%rownames(ArchRProject)] %>% setkey(cell) %>% .[rownames(ArchRProject@cellColData)] %>%
  as.data.frame() %>% tibble::column_to_rownames("cell")

# stopifnot(all(metadata.to.archR$TSSEnrichment_atac == getCellColData(ArchRProject,"TSSEnrichment")[[1]]))

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

fwrite(sample_metadata_tosave, args$outfile, sep="\t", na="NA", quote=F)

saveArchRProject(ArchRProject)
