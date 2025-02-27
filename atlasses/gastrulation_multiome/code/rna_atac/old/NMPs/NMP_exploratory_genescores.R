########################
## Load ArchR Project ##
########################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/atac/archR/load_archR_project.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/atac/archR/load_archR_project.R")
} else {
  stop("Computer not recognised")
}

#####################
## Define settings ##
#####################

# I/O
io$outdir <- paste0(io$basedir,"/results/rna_atac/NMPs")

# Options
opts$samples <- c(
  # "E7.5_rep1",
  # "E7.5_rep2",
  "E8.5_rep1",
  "E8.5_rep2"
)

opts$celltypes = c(
  # "Intermediate_mesoderm",
  # "Caudal_Mesoderm",
  # "Paraxial_mesoderm",
  "Somitic_mesoderm",
  "NMP",
  # "Forebrain_Midbrain_Hindbrain",
  "Spinal_cord"
)

#####################
## Update metadata ##
#####################

sample_metadata <- fread(io$metadata) %>%
  .[pass_atacQC==TRUE & pass_rnaQC==TRUE] %>%
  .[sample%in%opts$samples & celltype.predicted%in%opts$celltypes]
stopifnot(sample_metadata$archR_cell %in% rownames(ArchRProject))

##################
## Subset ArchR ##
##################

ArchRProject.filt <- ArchRProject[sample_metadata$archR_cell,]

###############
## Load data ##
###############

# Load RNA-based trajectory
io$nmp.pseudotime <- "/Users/ricard/data/gastrulation_multiome_10x/results/rna/NMP_trajectory/NMP_trajectory.txt.gz"
nmp.trajectory <- fread(io$nmp.pseudotime)# %>%
  # merge(sample_metadata[,c("cell","archR_cell")]) %>% 
  # .[,cell:=NULL] %>% setnames("archR_cell","cell")

# Load highly variable along the NMP trajectory
io$nmp.hvgs <- "/Users/ricard/data/gastrulation_multiome_10x/results/rna/NMP_trajectory/hvgs.rds"
hvgs <- readRDS(io$nmp.hvgs)

# Load SingleCellExperiment
sce <- load_SingleCellExperiment(io$sce, cells = sample_metadata$cell, genes=hvgs, normalise = TRUE, remove_non_expressed_genes = TRUE)

# Load gene accessibility scores
gene.score.se <- getMatrixFromProject(ArchRProject.filt, useMatrix = "GeneScoreMatrix_nodistal")
gene.score.se.filt <- gene.score.se %>%
  .[,sample_metadata$archR_cell] %>%
  .[Matrix::rowMeans(assay(.))>0.5,] 
gene.score.matrix <- assay(gene.score.se.filt) %>% as.matrix %>% t
colnames(gene.score.matrix) <- rowData(gene.score.se.filt)$name

# Subset genes
gene.score.matrix <- gene.score.matrix[,colnames(gene.score.matrix)  %in% hvgs]
dim(gene.score.matrix)

# chromVAR scores
# atac.deviation.mtx <- readRDS(io$archR.deviations.se)[,sample_metadata$archR_cell] %>% assay(.,"z")
# atac.deviation.mtx <- atac.deviation.mtx[apply(atac.deviation.mtx,1,var)>5,]

################
## parse data ##
################

acc_dt <- gene.score.matrix %>% as.data.table(keep.rownames = T) %>%
  setnames("rn","archR_cell") %>%
  melt(id.vars=c("archR_cell"), variable.name="gene", value.name="accessibility")

# chromvar_dt <- atac.deviation.mtx %>% t %>% as.data.table(keep.rownames = T) %>%
#   setnames("rn","cell") %>%
#   merge(nmp.trajectory,by="cell") %>%
#   melt(id.vars=c("cell","V1","V2"), variable.name="motif")

rna_dt <- logcounts(sce) %>% as.matrix %>% 
  as.data.table(keep.rownames = T) %>%
  setnames("rn","gene") %>%
  melt(id.vars="gene", variable.name="cell", value.name="expr") %>%
  merge(sample_metadata[,c("cell","celltype.predicted")],by="cell") %>%
  setnames("celltype.predicted","celltype")

##########
## Plot ##
##########

# For each gene, plot RNA expression vs chromatin accessibility along the trajectory

accrna_dt <- merge(
  rna_dt,
  acc_dt %>% merge(sample_metadata[,c("cell","archR_cell")],by="archR_cell"),
  by = c("cell","gene")
) %>% .[,archR_cell:=NULL] %>%
  merge(nmp.trajectory,by="cell")


genes.to.plot <- unique(to.plot$gene)

for (i in genes.to.plot) {
  
  to.plot <- accrna_dt[gene==i] %>% 
    .[accessibility>7,accessibility:=7] %>%
    melt(id.vars=c("cell","V1","V2","celltype"), measure.vars=c("accessibility","expr"), variable.name="modality") %>%
    .[,value_scaled:=value/max(value),by="modality"]
    
  p <- ggplot(to.plot, aes(x=V1, y=value_scaled)) +
    geom_point(aes(fill=modality), size=1.5, shape=21, stroke=0.1) +
    stat_smooth(aes(fill=modality), method="loess", color="black", alpha=0.75, span=0.5) +
    geom_rug(aes(color=celltype), sides="b") +
    scale_color_manual(values=opts$celltype.colors) +
    guides(color=F) +
    # scale_fill_manual(values=opts$celltype.colors) +
    # scale_fill_brewer(palette="Dark2") +
    labs(x="Pseudotime", y=sprintf("Scaled levels (%s)",i)) +
    theme_classic() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.title = element_blank(),
      legend.position="top"
    )
  
  pdf(sprintf("%s/%s_rna_acc_vs_pseudotime.pdf",io$outdir,i), width=8, height=4)
  print(p)
  dev.off()
}
