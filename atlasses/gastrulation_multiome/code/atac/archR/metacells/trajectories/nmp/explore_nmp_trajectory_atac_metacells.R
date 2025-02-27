
# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# opts$celltypes = c(
#   "Epiblast",
#   "Primitive_Streak",
#   "Caudal_mesoderm",
#   "Somitic_mesoderm",
#   "NMP",
#   "Spinal_cord"
# )

# I/O
io$cell_metadata <- io$metadata
io$metacell_metadata <- file.path(io$basedir,"results/rna/trajectories/nmp/nmp_sample_metadata.txt.gz")
io$trajectory <- file.path(io$basedir,"results/rna/trajectories/nmp/nmp_trajectory.txt.gz")
io$rna_metacells_sce <- file.path(io$basedir,"results/rna/metacells/trajectories/nmp/SingleCellExperiment_metacells.rds")
io$atac_peak_matrix_metacells <- file.path(io$basedir,"results/atac/archR/metacells/trajectories/nmp/PeakMatrix_summarized_experiment_metacells.rds")
io$atac_peak_diff_dir <- file.path(io$basedir,"results/atac/archR/differential/celltype.mapped/PeakMatrix")
# io$outdir <- paste0(io$basedir,"/results/atac/archR/celltype_hierarchies")

# options
opts$TFs <- fread(io$TFs_file)[[1]]

###################
## Load metadata ##
###################

cell_metadata.dt <- fread(io$cell_metadata) %>%
  .[pass_rnaQC==TRUE & pass_atacQC==TRUE & doublet_call==FALSE] %>%
  setnames("celltype.mapped","celltype")

# metacell_metadata.dt <- fread(io$metacell_metadata) %>%
#   .[,c("cell","celltype.mapped")] %>% setnames("celltype.mapped","celltype")

#####################
## Load trajectory ##
#####################

trajectory.dt <- fread(io$trajectory) %>% setnames(c("cell","V1","V2")) %>%
  merge(cell_metadata.dt[,c("cell","celltype")])

# Filter common cells
cells <- intersect(cell_metadata.dt$cell, trajectory.dt$cell)
cell_metadata.dt <- cell_metadata.dt[cell%in%cells]
trajectory.dt <- trajectory.dt[cell%in%cells]

###################
## Load RNA data ##
###################

# Load SingleCellExperiment for metacells
# rna_metacells.sce <- load_SingleCellExperiment(io$rna_metacells_sce, cells=metacell_metadata.dt$cell, normalise = TRUE)
rna_metacells.sce <- load_SingleCellExperiment(io$rna_metacells_sce, normalise = TRUE)

# Load SingleCellExperiment for cells
# sce <- load_SingleCellExperiment(io$rna.sce, cells=cell_metadata.dt$cell, normalise = TRUE)

# Add sample metadata as colData
# colData(sce) <- metacell_metadata.dt %>% tibble::column_to_rownames("cell") %>% DataFrame

####################
## Load ATAC data ##
####################

atac_peak_matrix_metacells.se <- readRDS(io$atac_peak_matrix_metacells)

# normalise
assay(atac_peak_matrix_metacells.se) <- 1e6*(sweep(assay(atac_peak_matrix_metacells.se),2,colSums(assay(atac_peak_matrix_metacells.se),na.rm=T),"/"))
assay(atac_peak_matrix_metacells.se) <- log2(assay(atac_peak_matrix_metacells.se)+0.5)

# colnames(atac_peak_matrix_metacells.se)

intersect(colnames(rna_metacells.sce), colnames(atac_peak_matrix_metacells.se))

#####################
## Define DE genes ##
#####################

io$rna.differential <- file.path(io$basedir,"results_new/rna/differential")
rna_diff.dt <- fread(paste0(io$rna.differential,"/Somitic_Mesoderm_vs_Spinal_cord.txt.gz")) %>%
  .[,sig:=padj_fdr<0.01 & abs(logFC)>2] %>% 
  .[,sign:="Up in Somitic_Mesoderm"] %>% 
  .[logFC<0,sign:=c("Up in Spinal_cord")]

genes.to.plot <- rna_diff.dt[sig==T,gene]

#####################
## Define DA peaks ##
#####################

atac.diff <- fread(paste0(io$atac_peak_diff_dir,"/PeakMatrix_Somitic_Mesoderm_vs_Spinal_cord.txt.gz")) %>%
  .[,sig:=FDR<0.01 & abs(MeanDiff)>=0.20] %>% 
  .[,sign:="Up in Somitic_Mesoderm"] %>% 
  .[MeanDiff<0,sign:=c("Up in Spinal_cord")]

peaks.to.plot <- atac.diff[sig==T,idx]

somitic_mesoderm_peaks <- atac.diff[sig==T & MeanDiff>=0.2,idx]
spinal_cord_peaks <- atac.diff[sig==T & MeanDiff<=0.2,idx]

#############
## Explore ##
#############

to.plot <- data.table(
  cell = colnames(atac_peak_matrix_metacells.se),
  somitic_mesoderm_peaks_acc = colMeans(assay(atac_peak_matrix_metacells.se[somitic_mesoderm_peaks,])),
  spinal_cord_peaks_acc = colMeans(assay(atac_peak_matrix_metacells.se[spinal_cord_peaks,]))
) %>% merge(trajectory.dt, by="cell")

ggplot(to.plot, aes_string(x="V1", y="V2", color ="celltype")) +
  geom_point(size=2.5) +
  scale_color_manual(values=opts$celltype.colors) +
  theme_classic() +
  ggplot_theme_NoAxes() +
  theme(
    legend.position = "none"
  )

ggplot(to.plot, aes_string(x="V1", y="V2", color ="somitic_mesoderm_peaks_acc")) +
  geom_point(size=2.5) +
  scale_color_distiller(palette = "YlOrRd", direction=1) +
  theme_classic() +
  ggplot_theme_NoAxes() +
  theme(
    legend.position = "none"
  )

ggplot(to.plot, aes_string(x="spinal_cord_peaks_acc", y="somitic_mesoderm_peaks_acc", color ="celltype")) +
  geom_point(size=2.5) +
  scale_color_manual(values=opts$celltype.colors) +
  theme_classic() +
  theme(
    legend.position = "none"
  )

