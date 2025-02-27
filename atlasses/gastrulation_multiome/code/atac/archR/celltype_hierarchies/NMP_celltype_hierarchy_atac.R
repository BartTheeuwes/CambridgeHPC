
# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

opts$celltypes = c(
  "Epiblast",
  "Primitive_Streak",
  # "Caudal_epiblast",
  # "Caudal_Mesoderm",
  "NMP",
  "Somitic_mesoderm",
  "Spinal_cord"
)

# I/O

io$rna_sce_pseudobulk_file <- file.path(io$basedir,"results_new/rna/pseudobulk/SingleCellExperiment_pseudobulk_celltype.mapped_mnn.rds")
io$rna_sce_tf_pseudobulk_file <- file.path(io$basedir,"results_new/rna/pseudobulk/SingleCellExperiment_CISBP_pseudobulk_celltype.mapped_mnn.rds")
io$diff_atac_peak_matrix <- file.path(io$basedir,"/results_new/atac/archR/differential/PeakMatrix/PeakMatrix_Somitic_Mesoderm_vs_Spinal_cord.txt.gz")
io$atac_chromvar_chip_pseudobulk_file <- file.path(io$basedir,"results_new/atac/archR/chromvar_chip/pseudobulk/chromVAR_deviations_CISBP_archr_chip.rds")
io$atac_peak_matrix <- file.path(io$basedir,"results_new/atac/archR/pseudobulk/celltype.mapped_mnn/pseudobulk_PeakMatrix_summarized_experiment.rds")
# io$outdir <- paste0(io$basedir,"/results_new/atac/archR/celltype_hierarchies")

########################
## Load cell metadata ##
########################

sample_metadata <- fread(io$metadata) %>%
  .[pass_rnaQC==TRUE & pass_atacQC==TRUE & doublet_call==FALSE & celltype.mapped_mnn%in%opts$celltypes] %>%
  setnames("celltype.mapped_mnn","celltype")
# .[pass_rnaQC==TRUE]

###############################
## Load single-cell RNA data ##
###############################

# Load RNA expression data as SingleCellExperiment object
# sce <- load_SingleCellExperiment(io$rna.sce, cells=sample_metadata$cell, normalise = TRUE)

# Add sample metadata as colData
# colData(sce) <- sample_metadata %>% tibble::column_to_rownames("cell") %>% DataFrame

#######################################
## Load pseudobulk RNA and ATAC data ##
#######################################

# Load pseudobulk TF RNA expression
rna_pseudobulk.sce <- readRDS(io$rna_sce_pseudobulk_file)[,opts$celltypes]
rna_pseudobulk_tf.sce <- readRDS(io$rna_sce_tf_pseudobulk_file)[,opts$celltypes]

# Load pseudobulk ATAC chromVAR-ChIP matrix
atac_pseudobulk_chromvar.se <- readRDS(io$atac_chromvar_chip_pseudobulk_file)[,opts$celltypes]

# Load pseudobulk ATAC peak matrix
atac_pseudobulk_PeakMatrix.se <- readRDS(io$atac_peak_matrix)[,opts$celltypes]

########################
## Prepare data table ##
########################

# atac_chromvar_pseudobulk.dt <- assay(atac_pseudobulk_chromvar.se) %>% t %>%
#   as.data.table(keep.rownames = T) %>%
#   setnames("rn","celltype") %>%
#   melt(id.vars=c("celltype"), variable.name="gene", value.name="chromvar_zscore")
# 
# rna_tf_pseudobulk.dt <- logcounts(rna_pseudobulk_tf.sce) %>%
#   as.data.table(keep.rownames = T) %>%
#   setnames("rn","gene") %>%
#   data.table::melt(id.vars="gene", variable.name="celltype", value.name="expr")

##################
## Define genes ##
##################

rna_diff.dt <- fread(paste0(io$rna.differential,"/Somitic_Mesoderm_vs_Spinal_cord.txt.gz")) %>%
  .[,sig:=padj_fdr<0.01 & abs(logFC)>2] %>% 
  .[,sign:="Up in Somitic_Mesoderm"] %>% 
  .[logFC<0,sign:=c("Up in Spinal_cord")]

genes.to.plot <- rna_diff.dt[sig==T,gene]

# Create data.table
rna_pseudobulk.dt <- logcounts(rna_pseudobulk.sce[genes.to.plot,]) %>%
  as.data.table(keep.rownames = T) %>%
  setnames("rn","gene") %>%
  data.table::melt(id.vars="gene", variable.name="celltype", value.name="expr")

rna_tf_pseudobulk.dt <- logcounts(rna_pseudobulk_tf.sce) %>%
  as.data.table(keep.rownames = T) %>%
  setnames("rn","gene") %>%
  data.table::melt(id.vars="gene", variable.name="celltype", value.name="expr")

# rna_cells.dt <- as.matrix(logcounts(sce[genes.to.plot])) %>%
#   as.data.table(keep.rownames = T) %>%
#   setnames("rn","gene") %>%
#   data.table::melt(id.vars="gene", variable.name="cell", value.name="expr")

##################
## Define peaks ##
##################

atac.diff <- fread(io$diff_atac_peak_matrix) %>%
  .[,sig:=FDR<0.01 & abs(MeanDiff)>=0.25] %>% 
  .[,sign:="Up in Somitic_Mesoderm"] %>% 
  .[MeanDiff<0,sign:=c("Up in Spinal_cord")]

stopifnot(atac.diff$idx%in%rownames(atac_pseudobulk_PeakMatrix.se))

peaks.to.plot <- atac.diff[sig==T,idx]

# Create data.table
atac_chromvar_pseudobulk.dt <- assay(atac_pseudobulk_PeakMatrix.se[peaks.to.plot,]) %>% t %>%
  as.data.table(keep.rownames = T) %>%
  setnames("rn","celltype") %>%
  melt(id.vars=c("celltype"), variable.name="idx", value.name="value")

##########
## Plot ##
##########

to.plot <- merge(
  atac.diff[idx%in%peaks.to.plot,c("idx","sign")], 
  atac_chromvar_pseudobulk.dt, 
  by = "idx", 
  allow.cartesian = TRUE
)

order.celltypes <- to.plot %>% 
  .[,.(value=mean(value)),by=c("celltype","sign")] %>% 
  dcast(celltype~sign) %>% 
  .[,diff:=`Up in Somitic_Mesoderm` - `Up in Spinal_cord`] %>%
  setorder(-diff) %>% .$celltype

# to.plot[,celltype:=factor(celltype,levels=order.celltypes)]
to.plot[,celltype:=factor(celltype,levels=opts$celltypes)]

ggboxplot(to.plot, x="celltype", y="value", fill="sign", outlier.shape=NA) +
  coord_cartesian(ylim=c(0,1)) +
  theme_classic() +
  labs(x="", y="Average accessibility") +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    # axis.text.x = element_text(color="black", angle=30, hjust=1),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black"),
  )



##########
## Plot ##
##########

to.plot <- merge(
  rna_diff.dt[gene%in%genes.to.plot,c("gene","sign")], 
  rna_pseudobulk.dt, 
  by = "gene", 
  allow.cartesian = TRUE
)

order.celltypes <- to.plot %>% 
  .[,.(expr=mean(expr)),by=c("celltype","sign")] %>% 
  dcast(celltype~sign) %>% 
  .[,diff:=`Up in Somitic_Mesoderm` - `Up in Spinal_cord`] %>%
  setorder(-diff) %>% .$celltype

to.plot[,celltype:=factor(celltype,levels=order.celltypes)]

ggboxplot(to.plot, x="celltype", y="expr", fill="sign", outlier.shape=NA) +
  # coord_cartesian(ylim=c(0,1)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(color="black", angle=30, hjust=1),
  )


##########
## TEST ##
##########

corr.mtx <- psych::corr.test(t(as.matrix(logcounts(sce[genes.to.plot]))),t(as.matrix(logcounts(sce[genes.to.plot]))), ci = F)$r
diag(corr.mtx) <- NA
corr.mtx["T","Sox2"]

to.plot <- rna_cells.dt %>% merge(sample_metadata[,c("cell","celltype")])

to.plot2 <- to.plot[gene%in%c("T","Sox2")] %>% dcast(cell+celltype~gene,value.var="expr")

ggscatter(to.plot2, x="Sox2", y="T", fill="celltype", size=2, shape=21) +
  scale_fill_manual(values=opts$celltype.colors) +
  theme(
    legend.position = "none"
  )


foo <- fread("/Users/argelagr/data/gastrulation_multiome_10x/results_new/rna/differential/pseudobulk/TFs/Somitic_mesoderm_vs_Spinal_cord.txt.gz") %>% .[abs(diff)>=3]
to.plot <- rna_tf_pseudobulk.dt[gene%in%foo$gene] %>% dcast(celltype~gene,value.var="expr")

ggscatter(to.plot, x="SOX2", y="T", fill="celltype", size=5, shape=21) +
  scale_fill_manual(values=opts$celltype.colors) +
  theme(
    legend.position = "none"
  )

for (i in foo$gene) {
  to.plot <- rna_tf_pseudobulk.dt[gene==i]
  
  p <- ggbarplot(to.plot, x="celltype", y="expr", fill="celltype") +
    scale_fill_manual(values=opts$celltype.colors) +
    labs(x="", title=i) +
    theme(
      axis.text.x = element_blank(),
      legend.position = "none"
    )
  print(p)
}
