
####################
## Load libraries ##
####################

#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
  source("/Users/ricard/gastrulation_multiome_10x/utils.R")
  source("/Users/ricard/gastrulation_multiome_10x/rna_atac/gene_regulatory_networks/DORC/pseudobulk/utils.R")
  source("/Users/ricard/gastrulation_multiome_10x/rna_atac/gene_regulatory_networks/DORC/pseudobulk/plot_utils.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
  source("/homes/ricard/gastrulation_multiome_10x/utils.R")
  source("/homes/ricard/gastrulation_multiome_10x/rna_atac/gene_regulatory_networks/DORC/pseudobulk/utils.R")
  source("/homes/ricard/gastrulation_multiome_10x/rna_atac/gene_regulatory_networks/DORC/pseudobulk/plot_utils.R")
} else if (grepl("Workstation",Sys.info()['nodename'])){
  source("/home/lijingyu/gastrulation/gastrulation_multiome_10x/settings.R")
  source("/home/lijingyu/gastrulation/gastrulation_multiome_10x/utils.R")
  source("/home/lijingyu/gastrulation/gastrulation_multiome_10x/rna_atac/gene_regulatory_networks/DORC/pseudobulk/utils.R")
  source("/home/lijingyu/gastrulation/gastrulation_multiome_10x/rna_atac/gene_regulatory_networks/DORC/pseudobulk/plot_utils.R")
} else {
  stop("Computer not recognised")
}

# Define I/O
io$rna.sce <- paste0(io$basedir, '/processed/rna/pseudobulk/SingleCellExperiment.rds')
io$dorc.pseudobulk <- paste0(io$basedir, '/results/rna_atac/DORCs/pseudobulk/DORCs_samples_distance10000.rds')
io$atac.peaks.se <- paste0(io$basedir,"/processed/atac/archR/pseudobulk/pseudobulk_PeakMatrix_summarized_experiment.rds")
io$bigwig <-  paste0(io$basedir,'/processed/atac/archR/GroupBigWigs')
io$outdir <- paste0(io$basedir,"/results/rna_atac/DORCs/pseudobulk/pdf")

# Options
opts$remove_ExE_celltypes <- FALSE

celltypes.to.plot <- opts$celltypes[1:3]


#############################
## Load pre-computed DORCs ##
#############################

DORC_output <- readRDS(io$dorc.pseudobulk)
# sum(is.na(DORC_output))

##########################
## Load atac.peaks data ##
##########################

atac.peaks.se <- readRDS(io$atac.peaks.se)
peak <- atac.peaks.se@elementMetadata
rownames(atac.peaks.se) <- paste0(peak$seqnames,':',peak$start,'-',peak$end)

########################
## Load ArchR Project ##
########################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/atac/archR/load_archR_project.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/atac/archR/load_archR_project.R")
} else if(grepl('Workstation',Sys.info()['nodename'])){
  # source("/home/lijingyu/gastrulation/gastrulation_multiome_10x/atac/archR/load_archR_project.R")
  ArchRProject=readRDS("/home/lijingyu/gastrulation/data/gastrulation_multiome_10x/processed/atac/archR/Save-ArchR-Project.rds")
}else{
  stop("Computer not recognised")
}

##############################
## Load RNA expression data ##
##############################

# Load SingleCellExperiment
sce <- readRDS(io$rna.pseudobulk.sce)
if (opts$remove_ExE_celltypes) {
  sce<- sce[,!(colnames(sce)%in% c("Visceral_endoderm","ExE_endoderm","ExE_ectoderm","Parietal_endoderm"))]
}

# Filter genes
sce <- sce[rownames(sce)%in%unique(DORC_output$gene),]

##########################
## Load gene annotation ##
##########################

gene_annotation <- getGeneAnnotation(ArchRProject)[["genes"]]
colnames(elementMetadata(gene_annotation)) <- c("gene_id","gene_name")

gene_annotation <- gene_annotation[!is.na(gene_annotation$gene_name)]

##########################################################
# Define bigwig files with pseudobulk ATAC measurements ##
##########################################################

bigwig <- list.files(io$bigwig,full.names = T)
names(bigwig) <- sapply(strsplit(list.files(io$bigwig),'-'), FUN = `[[`, 1)
bigwig <- as.list(bigwig)
bigwig <- bigwig[names(bigwig)%in%celltypes.to.plot]


##########
## Plot ##
##########



p<-CoveragePlot(bigwig = bigwig,window = 300,
                annotation = gene_annotation,Links = DORC_output,region  =c('Sox17'),
                peak.use=rownames(atac.peaks.se) ,
                idents =idents,extend.upstream = 5000,extend.downstream=5000)

##########
## Test ##
##########

# 
# metadata = sample_metadata.filt
# region = "Sox17"
# annotation = gene_annotation
# Links = DORC_output
# expression.data = logcounts(sce)
# bigwig = bigwig
# group.by = "celltype.predicted"
# # fragment.path
# features = NULL
# assay = NULL
# show.bulk = TRUE
# peak.use
# anno= TRUE
# peaks = TRUE
# peaks.group.by = NULL
# ranges = NULL
# ranges.group.by = NULL
# ranges.title = "Ranges"
# links = TRUE
# tile = FALSE
# tile.size = 100
# tile.cells = 100
# window = 300
# extend.upstream = 5000
# extend.downstream = 5000
# ymax = NULL
# scale.factor = NULL
# cells = NULL
# idents = NULL
# sep = c("-", "-")
# heights = NULL
# bigwig.type='coverage'
# max.downsample = 3000
# downsample.rate = 0.1
# add_explot=F