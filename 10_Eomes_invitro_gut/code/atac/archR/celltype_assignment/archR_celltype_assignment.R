# Cell type assignments are done by mapping the RNA expression profiles to the PijuanSala2019 atlas.
# there are some cells which did not pass QC for RNA expression but did pass QC on the ATAC modality
# This script predicts the cell type label for these cells by using a kNN prediction

suppressPackageStartupMessages(library(argparse))
# suppressPackageStartupMessages(library(ArchR))

here::i_am("atac/archR/celltype_assignment/archR_celltype_assignment.R")

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--metadata',        type="character",                               help='Cell metadata file')
p$add_argument('--lsi',             type="character",                               help='File with the LSI coordinates')
p$add_argument('--umap',             type="character",                               help='File with the UMAP coordinates')
p$add_argument('--input_celltype_column',             type="character",             help='Input column for the celltype information')
p$add_argument('--output_celltype_column',             type="character",            help='Output column for the celltype information')
p$add_argument('--k',             type="integer",  default=25,                      help='Number of kNN for the celltype prediction')
p$add_argument('--outdir',          type="character",                               help='Output directory')

args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

source(here::here("settings.R"))
source(here::here("utils.R"))

## START TEST ##
# args$metadata <- file.path(io$basedir,"results/atac/archR/qc/sample_metadata_after_qc.txt.gz") # io$metadata
# args$lsi <- file.path(io$basedir,"results/atac/archR/dimensionality_reduction/lsi_PeakMatrix_nfeatures5000_ndims30.txt.gz")
# args$umap <- file.path(io$basedir,"results/atac/archR/dimensionality_reduction/umap_PeakMatrix_nfeatures5000_ndims30_neigh25_dist0.3.txt.gz")
# args$input_celltype_column <- "celltype.mapped_mnn"
# args$output_celltype_column <- "celltype.predicted"
# args$k <- 25
# args$outdir <- file.path(io$basedir,"results/atac/archR/celltype_assignment")
## END TEST ##

dir.create(args$outdir, showWarnings=F)

########################
## Load ArchR project ##
########################

# source(here::here("atac/archR/load_archR_project.R"))
# addArchRThreads(threads=args$threads) 

########################
## Load cell metadata ##
########################

sample_metadata <- fread(args$metadata)
sample_metadata.filt <- sample_metadata[pass_atacQC==TRUE & doublet_call==FALSE] %>% setnames(args$input_celltype_column,"celltype")
# sample_metadata.filt <- sample_metadata[pass_atacQC==TRUE] %>% setnames(args$input_celltype_column,"celltype")

if (mean(is.na(sample_metadata.filt$celltype))==0) {
  stop("No cells missing cell type assignments")
}

print(paste0("Fraction of cells missing cell type assignment: ",round(mean(is.na(sample_metadata.filt$celltype)),2)))

################################################
## Load pre-computed dimensionality reduction ##
################################################

# LSI
lsi.dt <- fread(args$lsi) 
stopifnot(lsi.dt$cell%in%sample_metadata$cell)
stopifnot(sample_metadata.filt$cell%in%lsi.dt$cell)
lsi.mtx <- lsi.dt %>% matrix.please

# UMAP
umap.dt <- fread(args$umap)
stopifnot(umap.dt$cell%in%sample_metadata$cell)
stopifnot(sample_metadata.filt$cell%in%umap.dt$cell)

###########################################################################
## Create neighbourhood graph for cells that are missing cell type label ##
###########################################################################

df.observed <- sample_metadata.filt[!is.na(celltype),c("cell","celltype")]
cells.missing.label <- sample_metadata.filt[is.na(celltype),cell]

X.observed <- lsi.mtx[df.observed$cell,]
X.missing <- lsi.mtx[cells.missing.label,]

# returns a list of 2 (N,K) matrices.:
# - nn.idx: 1-indexed indices
# - nn.dists: distances
knnObj <- nabor::knn(data = X.observed, query = X.missing, k = args$k)

####################
## kNN prediction ##
####################

neighbour.celltypes <- apply(knnObj$nn.idx, 1, function(x) df.observed$celltype[x]) %>% t
predicted.celltype <- apply(neighbour.celltypes, 1, function(x) getmode(x, 1:length(x)))

############################
## Update sample metadata ##
############################

dt.predictions <- data.table(
  cell = cells.missing.label,
  celltype.predicted = predicted.celltype
)

sample_metadata.updated <- fread(args$metadata)
if (args$output_celltype_column%in%colnames(sample_metadata.updated)) sample_metadata.updated[[args$output_celltype_column]] <- NULL
sample_metadata.updated <- sample_metadata.updated %>% 
  merge(dt.predictions,by="cell", all.x=TRUE) %>%
  .[is.na(celltype.predicted),celltype.predicted:=eval(as.name(args$input_celltype_column))]

# Save metadata
fwrite(sample_metadata.updated, file.path(args$outdir,"sample_metadata_after_celltype_assignment.txt.gz"), sep="\t", na="NA", quote=F)

##########
## Plot ##
##########

to.plot <- umap.dt %>%
  merge(sample_metadata.updated, by="cell") %>%
  setnames(args$input_celltype_column,"celltype")

p1 <- ggplot(to.plot[!is.na(celltype)], aes(x=umap1, y=umap2)) +
  geom_point(aes(fill=celltype), size=1.25, shape=21, color="black", stroke=0.05) +
  scale_fill_manual(values=opts$celltype.colors) +
  theme_classic() +
  ggplot_theme_NoAxes() +
  labs(title=sprintf("Cells that have cell type assignment (%s, N=%s)",args$input_celltype_column,to.plot[!is.na(celltype),.N])) +
  theme(legend.position = "none", plot.title = element_text(hjust=0.5, size=rel(0.9)))

p2 <- ggplot(to.plot, aes(x=umap1, y=umap2)) +
  geom_point(aes(fill=is.na(celltype), size=is.na(celltype), alpha=is.na(celltype)), shape=21, color="black", stroke=0.05) +
  scale_size_manual(values=c("TRUE"=1, "FALSE"=0.5)) +
  scale_fill_manual(values=c("TRUE"="red", "FALSE"="gray60")) +
  scale_alpha_manual(values=c("TRUE"=0.75, "FALSE"=0.25)) +
  theme_classic() +
  ggplot_theme_NoAxes() +
  labs(title=sprintf("Highlighting cells that do not have cell type assignment (in red, N=%s)",to.plot[is.na(celltype),.N])) +
  theme(legend.position = "none", plot.title = element_text(hjust=0.5, size=rel(0.9)))
  
p3 <- ggplot(to.plot[is.na(celltype)], aes(x=umap1, y=umap2)) +
  geom_point(aes(fill=celltype.predicted), size=1.25, shape=21, color="black", stroke=0.05) +
  scale_fill_manual(values=opts$celltype.colors) +
  theme_classic() +
  ggplot_theme_NoAxes() +
  labs(title=sprintf("Subset of cells after celltype prediction (%s, N=%s)",args$output_celltype_column,to.plot[is.na(celltype),.N])) +
  theme(legend.position = "none", plot.title = element_text(hjust=0.5, size=rel(0.9)))

p <- cowplot::plot_grid(plotlist=list(p1,p2,p3), nrow = 1)

png(paste0(args$outdir,"/celltype_assignment.png"), width=1400, height=500)
print(p)
dev.off()


