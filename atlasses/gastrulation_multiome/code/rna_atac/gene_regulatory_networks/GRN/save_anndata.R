#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
  source("/Users/ricard/gastrulation_multiome_10x/utils.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
  source("/homes/ricard/gastrulation_multiome_10x/utils.R")
} else if (grepl("Workstation",Sys.info()['nodename'])){
  source("/home/lijingyu/gastrulation/gastrulation_multiome_10x/settings.R")
  source("/home/lijingyu/gastrulation/gastrulation_multiome_10x/utils.R")
} else {
  stop("Computer not recognised")
}
# Options

# opts$celltypes = c(
#   # "Mixed_mesoderm",
#   "Haematoendothelial_progenitors",
#   "Blood_progenitors_1",
#   "Blood_progenitors_2",
#   "Erythroid1",
#   "Erythroid2",
#   "Erythroid3"
# )

# mesoderm
opts$celltypes = c(
 "Epiblast",
 "Primitive_Streak",
 "Nascent_mesoderm"
)


io$outdir <- paste0(io$basedir,"/results/rna_atac/rna_vs_acc/trajectories/mesoderm_trajectory")
io$trajectory <- paste0(io$basedir,"/results/rna/trajectories/mesoderm_trajectory/mesoderm_trajectory.txt.gz")
io$metadata <- paste0(io$basedir, "/sample_metadata.txt.gz")
io$anndata <- paste0(io$outdir,'/anndata.h5ad')
########################
## Load cell metadata ##
########################


sample_metadata <- fread(io$metadata) %>%
  .[(pass_atacQC==TRUE & pass_rnaQC==TRUE) & doublet_call==FALSE] %>%
  .[sample%in%opts$samples & celltype.mapped%in%opts$celltypes]
trajectory.dt <- fread(io$trajectory)

# Filter cells
trajectory.dt <- trajectory.dt[cell%in%sample_metadata$cell]
sample_metadata <- sample_metadata[cell%in%trajectory.dt$cell] %>% setkey(cell) %>% .[trajectory.dt$cell]
sample_metadata <- merge(sample_metadata,trajectory.dt)
sample_metadata <- sample_metadata[trajectory.dt$cell,]
stopifnot(!is.na(sample_metadata$sample))
stopifnot(sample_metadata$cell==trajectory.dt$cell)

#################
## save anndata##
#################

sce <- load_SingleCellExperiment(
  file = io$rna.sce, 
  cells = sample_metadata$cell,
  normalise = TRUE, 
  remove_non_expressed_genes = FALSE
)
sce$celltype.mapped <- sample_metadata$celltype.mapped
sce$PC1 <- sample_metadata$PC1
sce$DC1 <- sample_metadata$DC1
marker_genes.dt <- fread(io$rna.atlas.marker_genes) 
sce <- sce[rownames(sce)%in%unique(marker_genes.dt$gene),]
ad <- AnnData(
  X = t(counts(sce)),
  obs = as.data.frame(sce@colData),
  uns = list(
    celltype.mapped_colors=opts$celltype.colors[sort(unique(sce$celltype.mapped))] )
)


ad$write_h5ad(io$anndata)


