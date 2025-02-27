
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
io$metadata <- paste0(io$basedir,"/sample_metadata.txt.gz")
# io$metadata <- paste0(io$basedir,"/results/atac/archR/celltype_assignment/sample_metadata_after_archR.txt.gz")
io$outdir <- paste0(io$basedir,"/results/atac/archR/dimensionality_reduction")

opts$samples <- c(
  "E7.5_rep1",
  "E7.5_rep2"
  # "E8.5_rep1",
  # "E8.5_rep2"
)

opts$celltypes = c(
  "Nascent_mesoderm",
  "Mixed_mesoderm",
  "Intermediate_mesoderm",
  "Caudal_Mesoderm",
  "Paraxial_mesoderm",
  "Somitic_mesoderm",
  "Pharyngeal_mesoderm",
  "Cardiomyocytes",
  "Allantois",
  "ExE_mesoderm"
  # "Mesenchyme"
)

# LSI parameters
# opts$lsi.iterations <- 1
# opts$lsi.cluster.resolution <- 2
opts$lsi.iterations = 2
opts$lsi.cluster.resolution = 2
opts$lsi.varFeatures <- 5000
opts$lsi.dims <- 15
opts$matrix <- "PeakMatrix"

# UMAP parameters
opts$umap.seed <- 42
opts$umap.neighbours <- 25
opts$umap.minDist <- 0.3

opts$batch.correction <- FALSE

########################
## Load cell metadata ##
########################

sample_metadata <- fread(io$metadata) %>%
  .[pass_atacQC==TRUE] %>%
  .[sample%in%opts$samples & celltype.mapped%in%opts$celltypes]

stopifnot(sample_metadata$archR_cell %in% rownames(ArchRProject))

##################
## Subset ArchR ##
##################

ArchRProject.filt <- ArchRProject[sample_metadata$archR_cell]
table(getCellColData(ArchRProject.filt,"Sample")[[1]])
table(getCellColData(ArchRProject.filt,"celltype.mapped")[[1]])

###########################
## Latent Semantic Index ##
###########################

# Iterative LSI: two iterations
ArchRProject.filt <- addIterativeLSI(ArchRProject.filt,
  useMatrix = opts$matrix, 
  name = "IterativeLSI", 
  firstSelection = "Top", # "Top" or "Var"
  dimsToUse = 1:opts$lsi.dims,
  depthCol = "nFrags",
  iterations = opts$lsi.iterations, 
  clusterParams = list(
    resolution = opts$lsi.cluster.resolution, 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = opts$lsi.varFeatures, 
  force = TRUE
)

# Correlation between latent dimensions and number of peaks
cor(getReducedDims(ArchRProject.filt, "IterativeLSI"),ArchRProject.filt$nFrags)[,1] %>% abs %>% sort(decreasing = T)

##########
## UMAP ##
##########

ArchRProject.filt <- addUMAP(ArchRProject.filt, 
  reducedDims = "IterativeLSI",
  dimsToUse = 2:opts$lsi.dims,
  name = "UMAP",
  metric = "cosine",
  nNeighbors = opts$umap.neighbours, 
  minDist = opts$umap.minDist, 
  seed = opts$umap.seed,
  force = TRUE
)
head(ArchRProject.filt@embeddings$UMAP$df)

##############
## Plot LSI ##
##############

to.plot <- getReducedDims(ArchRProject.filt,"IterativeLSI") %>%
  as.data.table(keep.rownames = T) %>%
  setnames("rn","archR_cell") %>% merge(sample_metadata,by="archR_cell")

p <- ggplot(to.plot, aes(x=LSI6, y=LSI)) +
  # geom_point(aes(fill=nFrags_atac), size=2, shape=21, color="black", stroke=0.05) +
  geom_point(aes(fill=celltype.mapped), size=2, shape=21, color="black", stroke=0.05) +
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

##############
## Plot UMAP ##
##############

to.plot <- getEmbedding(ArchRProject.filt,"UMAP") %>%
  as.data.table(keep.rownames = T) %>%
  setnames(c("archR_cell","umap1","umap2")) %>% merge(sample_metadata,by="archR_cell")

p1 <- ggplot(to.plot, aes(x=umap1, y=umap2)) +
  geom_point(aes(fill=celltype.mapped), size=2, shape=21, color="black", stroke=0.05) +
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
print(p1)
# dev.off()

p2 <- ggplot(to.plot, aes(x=umap1, y=umap2)) +
  geom_point(aes(color=log2(nFrags_atac)), size=1) +
  scale_color_gradient(low = "gray80", high = "red") +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank()
  )

# pdf(sprintf("%s/archr_umap_nFragsATAC_%s_LSIiter%s_nfeatures%s_neighb%s_mindist%s.pdf",io$outdir,opts$matrix,opts$lsi.iterations, opts$lsi.varFeatures, opts$umap.neighbours,opts$umap.minDist))
print(p2)
# dev.off()

p3 <- ggplot(to.plot, aes(x=umap1, y=umap2)) +
  geom_point(aes(fill=sample), size=1.5, shape=21, color="black", stroke=0.05) +
  guides(fill = guide_legend(override.aes = list(size=4))) +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    legend.position = "top",
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank()
  )

# pdf(sprintf("%s/archr_umap_sample_%s_LSIiter%s_nfeatures%s_neighb%s_mindist%s.pdf",io$outdir,opts$matrix,opts$lsi.iterations, opts$lsi.varFeatures, opts$umap.neighbours,opts$umap.minDist))
print(p3)
# dev.off()
