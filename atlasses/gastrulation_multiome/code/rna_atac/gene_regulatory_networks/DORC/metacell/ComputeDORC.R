
####################
## Load libraries ##
####################
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(future))
suppressPackageStartupMessages(library(pbapply))
suppressPackageStartupMessages(library(future.apply))
suppressPackageStartupMessages(library(Matrix))
library(ArchR)
#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
  source("/Users/ricard/gastrulation_multiome_10x/utils.R")
  source("/Users/ricard/gastrulation_multiome_10x/rna_atac/gene_regulatory_networks/DORC/metacell/utils.R")
  source("/Users/ricard/gastrulation_multiome_10x/rna_atac/gene_regulatory_networks/DORC/metacell/plot_utils.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
  source("/homes/ricard/gastrulation_multiome_10x/utils.R")
  source("/homes/ricard/gastrulation_multiome_10x/rna_atac/gene_regulatory_networks/DORC/metacell/utils.R")
  source("/homes/ricard/gastrulation_multiome_10x/rna_atac/gene_regulatory_networks/DORC/metacell/plot_utils.R")
} else if (grepl("Workstation",Sys.info()['nodename'])){
  source("/home/lijingyu/gastrulation/gastrulation_multiome_10x/settings.R")
  source("/home/lijingyu/gastrulation/gastrulation_multiome_10x/utils.R")
  source("/home/lijingyu/gastrulation/gastrulation_multiome_10x/rna_atac/gene_regulatory_networks/DORC/metacell/utils.R")
  source("/home/lijingyu/gastrulation/gastrulation_multiome_10x/rna_atac/gene_regulatory_networks/DORC/metacell/plot_utils.R")
} else {
  stop("Computer not recognised")
}

# Define I/O
io$sce <- paste0(io$basedir, '/processed/rna/metacell/SingleCellExperiment_2500metacells.rds')
io$sce <- paste0(io$basedir, '/processed/rna/metacell/SingleCellExperiment_1000metacells.rds')
io$atac.peaks.se <- paste0(io$basedir,"/processed/atac/archR/atac_SummarizedExperiment.rds")
io$bigwig <-  paste0(io$basedir,'/processed/atac/archR/GroupBigWigs')


#####################
## Load metadata ##
#####################

sample_metadata <- fread(io$metadata) %>%
  .[pass_rnaQC==TRUE & doublet_call==FALSE] %>%
  .[celltype.mapped%in%opts$celltypes & sample%in%opts$samples]

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
# p$add_argument('--samples',      type="character",    nargs="+",  help='Samples')
# p$add_argument('--celltypes',  type="character",    nargs="+",  help='Cell type')
p$add_argument('--distance',      type="character",         default=1e4  ,       help='Distance')
p$add_argument('--ncores',       type="integer",      default=1,  help='Number of cores')
# p$add_argument('--denoised',     action="store_true",             help='Use denoised ATAC data?')
p$add_argument('--remove_ExE_celltypes', action="store_true",  help='Remove ExE cell types?')
p$add_argument('--test_mode',    action="store_true",            help='Test mode? subset number of genes')
p$add_argument('--outdir',       type="character",      default=paste0(io$basedir,"/results/rna_atac/GRN/DORC/metacell") ,        help='Output file')
args <- p$parse_args(commandArgs(TRUE))

## START TEST
args$samples <- c("E7.5_rep1", "E7.5_rep2", "E8.0_rep1", "E8.0_rep2", "E8.5_rep1", "E8.5_rep2")
args$distance <- 1e4
args$ncores <- 1
args$remove_ExE_celltypes <- FALSE
args$test_mode <- T
args$outdir <- paste0(io$basedir,"/results/rna_atac/GRN/DORC/metacell")
## END TEST

# Define multiprocessing
options(future.globals.maxSize = 10 * 1024^3)
if (args$ncores>1) {
  plan(multiprocess, workers=args$ncores)
} else {
  plan(sequential)
}
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
sce <- load_SingleCellExperiment(io$sce, normalise = TRUE, remove_non_expressed_genes = TRUE)
if (args$remove_ExE_celltypes) {
  sce<- sce[,!(colnames(sce)%in% c("Visceral_endoderm","ExE_endoderm","ExE_ectoderm","Parietal_endoderm"))]
}
# Load cell type markers
marker_genes.dt <- fread(io$rna.atlas.marker_genes) 

# Filter genes
sce <- sce[!grepl("Rik|Gm",rownames(sce)),]
sce <- sce[rownames(sce)%in%unique(marker_genes.dt$gene),]

##########################
## Load atac.peaks data ##
##########################


atac.peaks.se=readRDS(io$atac.peaks.se)
atac.peaks.se=atac.peaks.se[,colnames(atac.peaks.se)%in%colnames(sce)]
sce<-sce[,colnames(atac.peaks.se)]

##########################
## Load gene annotation ##
##########################

## run if use gene_annotation of EnsDb
# library(EnsDb.Mmusculus.v79)
# suppressWarnings(gene_annotation <- GetGRangesFromEnsDb(EnsDb.Mmusculus.v79))
# seqlevelsStyle(gene_annotation) <- "UCSC"
# names(gene_annotation)<-NULL

gene_annotation <- getGeneAnnotation(ArchRProject)[["genes"]]
colnames(elementMetadata(gene_annotation)) <- c("gene_id","gene_name")
gene_annotation=gene_annotation[!is.na(gene_annotation$gene_name)]
if (args$test_mode) {
  gene_annotation <- gene_annotation[1:1000]
}
# names(gene_annotation)=NULL

####################
## Define peakSet ##
####################

peakSet.gr <- getPeakSet(ArchRProject)

####################
## Load Bgdpeaks  ##
####################
background.peaks.se <- readRDS(io$archR.bgdPeaks)
# background.peaks.se <- getBgdPeaks(ArchRProject)

#########
## Run ##
#########

# plan(multisession, gc = TRUE, workers = 1)

DORC_output <- ComputeDORC(
  atac.peaks.se=atac.peaks.se,
  background.peaks.se=background.peaks.se,
  expression.data  = logcounts(sce),
  regions = peakSet.gr,
  gene.coords = gene_annotation,
  distance = args$distance
  # min.cells = 0,
  # method = "pearson",
  # n_sample = 50,
  # pvalue_cutoff =0.05,
  # score_cutoff = 0.05,
  # sep=c(":", "-"),
  # verbose = TRUE
)
# Save (TO-DO: SAVE OPTIONS )
if (args$remove_ExE_celltypes) {
  outfile <- sprintf("%s/DORCs_samples%s_distance%d_noExE.rds", args$outdir, paste(args$samples,collapse="-"), args$distance)
} else {
  outfile <- sprintf("%s/DORCs_samples%s_distance%d.rds", args$outdir, paste(args$samples,collapse="-"), args$distance)
}
saveRDS(DORC_output, outfile)

#########
## Plot##
#########

# # choose the celltypes to plot
# idents=opts$celltypes[1:7]
# # idents=c('Epiblast',"Visceral_endoderm","Parietal_endoderm")
# 
# # load bigwig
# bigwig=list.files(io$bigwig,full.names = T)
# names(bigwig)=sapply(strsplit(list.files(io$bigwig),'-'), FUN = `[[`, 1)
# bigwig=as.list(bigwig)
# bigwig=bigwig[names(bigwig)%in%idents]
# 
# 
# # plot 
# p=CoveragePlot(metadata = sample_metadata,group.by = 'celltype.mapped',bigwig = bigwig,window = 300,
#              annotation = gene_annotation,Links = DORC_output,region  =c('Pou5f1'),
#              peak.use=rownames(atac.peaks.se) ,expression.data = logcounts(sce),
#              idents =idents,extend.upstream = 5000,extend.downstream=5000)
