suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(ggplot2))

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--gene',       type="character",  help='')
p$add_argument('--samples',      type="character",    nargs="+",  help='Samples')
# p$add_argument('--celltypes',      type="character",    nargs="+",  help='')
# p$add_argument('--celltypes',  type="character",    nargs="+",  help='Cell type')
p$add_argument('--remove_ExE_celltypes', action="store_true",   help='Remove ExE cell types?')
p$add_argument('--distance',      type="integer",                  help='Distance')
p$add_argument('--test_mode',    action="store_true",             help='Test mode? subset number of cells')
p$add_argument('--outdir',       type="character",                help='Output file')
args <- p$parse_args(commandArgs(TRUE))

## START TEST
args$gene <- "Hoxa9"
args$samples <- c("E7.5_rep1", "E7.5_rep2", "E8.0_rep1", "E8.0_rep2", "E8.5_rep1", "E8.5_rep2")
args$distance <- 1e4
args$remove_ExE_celltypes <- FALSE
args$test_mode <- TRUE
args$outdir <- "/Users/ricard/data/gastrulation_multiome_10x/results/rna_atac/DORCs_v2/pdf"
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

##############################
## Load RNA expression data ##
##############################

# Load SingleCellExperiment
sce <- load_SingleCellExperiment(io$rna.sce, cells = sample_metadata$cell, normalise = TRUE, remove_non_expressed_genes = TRUE)

# Filter genes
sce <- sce[!grepl("Rik|Gm",rownames(sce)),]

###########################
## Load ATAC peak matrix ##
###########################

# Fetch peak matrix
atac.peaks.se <- getMatrixFromProject(ArchRProject.filt, binarize = TRUE, useMatrix = "PeakMatrix")

# Rename peaks
rownames(atac.peaks.se) <- paste0(seqnames(rowRanges(atac.peaks.se)),":",ranges(rowRanges(atac.peaks.se)))

##########################
## Load gene annotation ##
##########################

gene_annotation <- getGeneAnnotation(ArchRProject.filt)[["genes"]]
colnames(elementMetadata(gene_annotation)) <- c("gene_id","gene_name")

if (args$test_mode) {
  gene_annotation <- gene_annotation[1:100]
}

##################
## Load peakSet ##
##################

peakSet.gr <- getPeakSet(ArchRProject.filt)

#######################
## Load DORC results ##
#######################

io$dorc.file <- paste0(io$basedir,"/results/rna_atac/DORCs_v2/DORCs_samplesE7.5_rep1-E7.5_rep2-E8.0_rep1-E8.0_rep2-E8.5_rep1-E8.5_rep2_distance10000_noExE.rds")
DORCs.gr <- readRDS(io$dorc.file)

##########
## Plot ##
##########

io$atac_fragments.path = "/hps/nobackup2/research/stegle/users/ricard/gastrulation_multiome_10x/original/fragments.tsv.gz"

CoveragePlot(
  metadata = sample_metadata,
  region = "Sox17",
  annotation = gene_annotation,
  peak.data = assay(atac.peaks.se),
  expression.data = logcounts(sce),
  group.by = "celltype.predicted",
  Links = DORCs.gr,
  window = 300,
  extend.upstream = 10000,
  extend.downstream = 10000,
  features = NULL,
  assay = NULL,
  show.bulk = TRUE,
  anno= TRUE,
  peaks = TRUE,
  peaks.group.by = NULL,
  ranges = NULL,
  ranges.group.by = NULL,
  ranges.title = "Ranges",
  links = TRUE,
  tile = FALSE,
  tile.size = 100,
  tile.cells = 100,
  ymax = NULL,
  scale.factor = NULL,
  cells = NULL,
  idents = NULL,
  sep = c("-", "-"),
  heights = NULL,
  max.downsample = 3000,
  downsample.rate = 0.1
)


##########
## Test ##
##########


metadata = sample_metadata %>% tibble::column_to_rownames("cell")
region = "Sox17"
annotation = gene_annotation
peak.data = assay(atac.peaks.se)
expression.data = logcounts(sce)
group.by = "celltype.predicted"
Links = DORCs.gr
fragment.path = io$atac_fragments.path
window = 300
extend.upstream = 10000
extend.downstream = 10000
features = NULL
assay = NULL
show.bulk = TRUE
anno= TRUE
peaks = TRUE
peaks.group.by = NULL
ranges = NULL
ranges.group.by = NULL
ranges.title = "Ranges"
links = TRUE
tile = FALSE
tile.size = 100
tile.cells = 100
ymax = NULL
scale.factor = NULL
cells = NULL
idents = NULL
sep = c("-", "-")
heights = NULL
max.downsample = 3000
downsample.rate = 0.1

fragments <- sprintf("%s/%s_atac_fragments.tsv.gz",io$atac_fragments.path,opts$samples)
