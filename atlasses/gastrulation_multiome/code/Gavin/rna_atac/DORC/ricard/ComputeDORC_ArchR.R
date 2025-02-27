suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(future))
suppressPackageStartupMessages(library(pbapply))
suppressPackageStartupMessages(library(future.apply))
suppressPackageStartupMessages(library(Matrix))

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--samples',      type="character",    nargs="+",  help='Samples')
# p$add_argument('--celltypes',  type="character",    nargs="+",  help='Cell type')
p$add_argument('--distance',      type="integer",                  help='Distance')
p$add_argument('--ncores',       type="integer",      default=1,  help='Number of cores')
# p$add_argument('--denoised',     action="store_true",             help='Use denoised ATAC data?')
p$add_argument('--remove_ExE_celltypes', action="store_true",   help='Remove ExE cell types?')
p$add_argument('--test_mode',    action="store_true",             help='Test mode? subset number of cells')
p$add_argument('--outdir',       type="character",                help='Output file')
args <- p$parse_args(commandArgs(TRUE))

## START TEST
# args$samples <- c("E7.5_rep1", "E7.5_rep2", "E8.0_rep1", "E8.0_rep2", "E8.5_rep1", "E8.5_rep2")
# args$distance <- 1e4
# args$ncores <- 1
# args$remove_ExE_celltypes <- FALSE
# args$test_mode <- TRUE
# args$outdir <- "/Users/ricard/data/gastrulation_multiome_10x/results/rna_atac/DORCs_v2"
## END TEST

########################
## Load ArchR Project ##
########################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/atac/archR/load_archR_project.R")
  source("/Users/ricard/gastrulation_multiome_10x/Gavin/DORC/ricard/utils.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/atac/archR/load_archR_project.R")
  source("/homes/ricard/gastrulation_multiome_10x/Gavin/DORC/ricard/utils.R")
} else {
  stop("Computer not recognised")
}

# Define multiprocessing
options(future.globals.maxSize = 10 * 1024^3)
if (args$ncores>1) {
  plan(multiprocess, workers=args$ncores)
} else {
  plan(sequential)
}

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(io$metadata) %>%
  .[pass_atacQC==TRUE & pass_rnaQC==TRUE & doublet_call==FALSE] %>%
  .[sample%in%args$samples]

if (args$remove_ExE_celltypes) {
  sample_metadata <- sample_metadata %>%
    .[!celltype.mapped%in%c("Visceral_endoderm","ExE_endoderm","ExE_ectoderm","Parietal_endoderm")]
}

if (args$test_mode) {
  sample_metadata <- sample_metadata %>% head(n=100)
}

table(sample_metadata$sample)

##################
## Subset ArchR ##
##################

ArchRProject.filt <- ArchRProject[sample_metadata$cell,]

# Updatea peakSet
# ArchRProject.filt@peakSet <- readRDS(io$archR.peakSet.granges)

##############################
## Load RNA expression data ##
##############################

# Load SingleCellExperiment
sce <- load_SingleCellExperiment(io$rna.sce, cells = sample_metadata$cell, normalise = TRUE, remove_non_expressed_genes = TRUE)

# Load cell type markers
marker_genes.dt <- fread(io$rna.atlas.marker_genes) 

# Filter genes
sce <- sce[!grepl("Rik|Gm",rownames(sce)),]
sce <- sce[rownames(sce)%in%unique(marker_genes.dt$gene),]

##########################
## Load gene annotation ##
##########################

gene_annotation <- getGeneAnnotation(ArchRProject.filt)[["genes"]]
colnames(elementMetadata(gene_annotation)) <- c("gene_id","gene_name")

if (args$test_mode) {
  gene_annotation <- gene_annotation[1:100]
}

##################
## Define peakSet ##
##################

peakSet.gr <- getPeakSet(ArchRProject.filt)

#########
## Run ##
#########

# plan(multisession, gc = TRUE, workers = 1)

DORC_output <- ComputeDORC_ArchR(
  ArchRProject = ArchRProject.filt,
  expression.data = logcounts(sce),
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
