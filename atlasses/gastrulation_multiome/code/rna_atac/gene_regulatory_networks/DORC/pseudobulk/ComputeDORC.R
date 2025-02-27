
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
io$sce <- paste0(io$basedir, '/processed/rna/pseudobulk/SingleCellExperiment.rds')
io$atac.peaks.se <- paste0(io$basedir,"/processed/atac/archR/pseudobulk/pseudobulk_PeakMatrix_summarized_experiment.rds")
io$bigwig <-  paste0(io$basedir,'/processed/atac/archR/GroupBigWigs')


################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
# p$add_argument('--samples',      type="character",    nargs="+",  help='Samples')
# p$add_argument('--celltypes',  type="character",    nargs="+",  help='Cell type')
p$add_argument('--distance',      type="integer",         default=1e4  ,       help='Distance')
p$add_argument('--ncores',       type="integer",      default=1,  help='Number of cores')
# p$add_argument('--denoised',     action="store_true",             help='Use denoised ATAC data?')
p$add_argument('--remove_ExE_celltypes', action="store_true",  help='Remove ExE cell types?')
p$add_argument('--test_mode',    action="store_true",            help='Test mode? subset number of genes')
p$add_argument('--outdir',       type="character",      default=paste0(io$basedir,"/results/rna_atac/GRN/DORC") ,        help='Output file')
args <- p$parse_args(commandArgs(TRUE))

## START TEST
# args$samples <- c("E7.5_rep1", "E7.5_rep2", "E8.0_rep1", "E8.0_rep2", "E8.5_rep1", "E8.5_rep2")
# args$distance <- 1e4
# args$ncores <- 1
# args$remove_ExE_celltypes <- FALSE
# args$test_mode <- F
# args$outdir <- paste0(io$basedir,"/results/rna_atac/GRN/DORC")
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

##########################
## Load atac.peaks data ##
##########################


atac.peaks.se<-readRDS(io$atac.peaks.se)
peak<-atac.peaks.se@elementMetadata
rownames(atac.peaks.se)<-paste0(peak$seqnames,':',peak$start,'-',peak$end)

##############################
## Load RNA expression data ##
##############################

# Load SingleCellExperiment
sce <- readRDS(io$sce)
if (args$remove_ExE_celltypes) {
  sce<- sce[,!(colnames(sce)%in% c("Visceral_endoderm","ExE_endoderm","ExE_ectoderm","Parietal_endoderm"))]
}
# Load cell type markers
marker_genes.dt <- fread(io$rna.atlas.marker_genes) 

# Filter genes
sce <- sce[!grepl("Rik|Gm",rownames(sce)),]
sce <- sce[rownames(sce)%in%unique(marker_genes.dt$gene),]

##########################
## Load gene annotation ##
##########################
# define gene annotation

## run if use gene_annotation of EnsDb
# library(EnsDb.Mmusculus.v79)
# suppressWarnings(annotation <- GetGRangesFromEnsDb(EnsDb.Mmusculus.v79))
# seqlevelsStyle(annotation) <- "UCSC"
# names(annotation)=NULL

gene_annotation <- getGeneAnnotation(ArchRProject)[["genes"]]
colnames(elementMetadata(gene_annotation)) <- c("gene_id","gene_name")
gene_annotation <- gene_annotation[!(is.na(gene_annotation$gene_name))]
if (args$test_mode) {
  gene_annotation <- gene_annotation[1:90]
}


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
  distance = args$distance,
  score_cutoff = 0.01,
  pvalue_cutoff = 1
  # min.cells = 0,
  # method = "pearson",
  # n_sample = 50,
  # pvalue_cutoff =0.05,
  # score_cutoff = 0.25,
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

# choose the celltypes to plot
idents<-opts$celltypes[1:7]
# idents<-c('Epiblast',"Visceral_endoderm","Parietal_endoderm")

# load bigwig
bigwig<-list.files(io$bigwig,full.names = T)
names(bigwig)<-sapply(strsplit(list.files(io$bigwig),'-'), FUN = `[[`, 1)
bigwig<-as.list(bigwig)
bigwig<-bigwig[names(bigwig)%in%idents]


# plot
p<-CoveragePlot(bigwig = bigwig,window = 300,
                annotation = gene_annotation,Links = DORC_output,region  =c('Sox17'),
                peak.use=rownames(atac.peaks.se) ,
                idents =idents,extend.upstream = 5000,extend.downstream=5000)

