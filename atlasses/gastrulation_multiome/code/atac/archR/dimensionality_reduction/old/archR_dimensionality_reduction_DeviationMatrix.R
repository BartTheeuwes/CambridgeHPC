
########################
## Load ArchR project ##
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

# io$metadata <- paste0(io$basedir,"/processed/atac/archR/sample_metadata_after_archR.txt.gz")
# io$metadata <- paste0(io$basedir,"/sample_metadata.txt.gz")
# io$metadata <- paste0(io$basedir,"/results/atac/archR/celltype_assignment/sample_metadata_after_archR.txt.gz")
io$outdir <- paste0(io$basedir,"/results/atac/archR/dimensionality_reduction/test")

opts$samples <- c(
  "E7.5_rep1",
  "E7.5_rep2",
  "E8.0_rep1",
  "E8.0_rep2",
  "E8.5_rep1",
  "E8.5_rep2"
)

opts$celltypes <- c(
  "Epiblast",
  "Primitive_Streak",
  "Caudal_epiblast",
  # "PGC",
  "Anterior_Primitive_Streak",
  "Notochord",
  "Def._endoderm",
  "Gut",
  "Nascent_mesoderm",
  "Mixed_mesoderm",
  "Intermediate_mesoderm",
  "Caudal_Mesoderm",
  "Paraxial_mesoderm",
  "Somitic_mesoderm",
  "Pharyngeal_mesoderm",
  "Cardiomyocytes",
  "Allantois",
  "ExE_mesoderm",
  "Mesenchyme",
  "Haematoendothelial_progenitors",
  "Endothelium",
  "Blood_progenitors_1",
  "Blood_progenitors_2",
  "Erythroid1",
  "Erythroid2",
  "Erythroid3",
  "NMP",
  "Rostral_neurectoderm",
  "Caudal_neurectoderm",
  "Neural_crest",
  "Forebrain_Midbrain_Hindbrain",
  "Spinal_cord",
  "Surface_ectoderm",
  "Visceral_endoderm",
  "ExE_endoderm",
  "ExE_ectoderm",
  "Parietal_endoderm"
)

opts$aggregated.celltypes <- c(
  # "Erythroid1" = "Erythroid",
  # "Erythroid2" = "Erythroid",
  # "Erythroid3" = "Erythroid",
  "Blood_progenitors_1" = "Blood_progenitors",
  "Blood_progenitors_2" = "Blood_progenitors",
  "Rostral_neurectoderm" = "Neurectoderm",
  "Caudal_neurectoderm" = "Neurectoderm",
  "Anterior_Primitive_Streak" = "Primitive_Streak"
)

opts$remove.ExE.celltypes <- F

opts$matrix <- "DeviationMatrix_Motif_cisbp"

# UMAP parameters
opts$umap.seed <- 42
opts$umap.neighbours <- 25
opts$umap.minDist <- 0.3

########################
## Load cell metadata ##
########################

sample_metadata <- fread(io$metadata) %>%
  .[pass_atacQC==TRUE & doublet_call==FALSE] %>%
  .[sample%in%opts$samples & celltype.predicted%in%opts$celltypes]

if (opts$remove.ExE.celltypes) {
  sample_metadata <- sample_metadata %>%
    .[!celltype.predicted%in%c("Visceral_endoderm","ExE_endoderm","ExE_ectoderm","Parietal_endoderm")]
}

sample_metadata %>%
  .[,aggregated_celltype:=stringr::str_replace_all(celltype.predicted,opts$aggregated.celltypes)]

opts$aggregated_celltype.colors <- opts$celltype.colors[names(opts$celltype.colors)%in%unique(sample_metadata$aggregated_celltype)]
stopifnot(unique(sample_metadata$aggregated_celltype)%in%names(opts$aggregated_celltype.colors))
sample_metadata %>%
  .[,aggregated_celltype:=factor(aggregated_celltype, levels=names(opts$celltype.colors))]

##################
## Subset ArchR ##
##################

ArchRProject.filt <- ArchRProject[sample_metadata$cell]

#############################
## Extract chromVAR matrix ##
#############################

chromvar.mtx <- getMatrixFromProject(ArchRProject.filt, useMatrix = opts$matrix) %>% assay(.,"z")

chromvar.mtx <- chromvar.mtx[apply(chromvar.mtx,1,var)>1,]

#########
## PCA ##
#########

pca.mtx <- irlba::prcomp_irlba(t(chromvar.mtx), n=50)$x
rownames(pca.mtx) <- colnames(chromvar.mtx)

##########
## UMAP ##
##########

umap.dt <- uwot::umap(pca.mtx, n_neighbors = opts$umap.neighbours, min_dist = opts$umap.minDist, metric = "cosine") %>%
  as.data.table(keep.rownames=F) %>%
  setnames(c("umap1","umap2")) %>%
  .[,cell:=rownames(pca.mtx)]

##########
## Plot ##
##########

to.plot <- umap.dt %>% merge(sample_metadata, by="cell")
  
p <- ggplot(to.plot, aes(x=umap1, y=umap2)) +
  # geom_point(aes(fill=celltype.predicted), size=2, shape=21, color="black", stroke=0.05) +
  geom_point(aes(fill=aggregated_celltype), size=2, shape=21, color="black", stroke=0.05) +
  scale_fill_manual(values=opts$celltype.colors) +
  guides(fill = guide_legend(override.aes = list(size=4))) +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank()
  )

# pdf(sprintf("%s/archr_umap_celltype_%s_LSIiter%s_nfeatures%s_neighb%s_mindist%s.pdf",io$outdir,opts$matrix,opts$lsi.iterations, opts$lsi.varFeatures, opts$umap.neighbours,opts$umap.minDist))
print(p)
# dev.off()
