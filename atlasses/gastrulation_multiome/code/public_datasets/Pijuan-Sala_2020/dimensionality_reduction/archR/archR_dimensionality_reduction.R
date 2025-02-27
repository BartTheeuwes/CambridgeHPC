
########################
## Load ArchR project ##
########################

source("/Users/ricard/gastrulation_multiome_10x/public_datasets/Pijuan-Sala_2020/archR/load_archR_project.R")

names(ArchRProject@embeddings)

#####################
## Define settings ##
#####################

io$outdir <- paste0(io$basedir,"/results/dimensionality_reduction/archR")

# trial 1
# opts$lsi.iterations <- 2
# opts$lsi.cluster.resolution <- c(0.2)
# opts$lsi.varFeatures <- 25000
# opts$lsi.dims <- 30

# trial 2
# opts$lsi.iterations <- 3
# opts$lsi.cluster.resolution <- c(1, 2)
# opts$lsi.varFeatures <- 50000
# opts$lsi.dims <- 50

# trial 1
opts$lsi.iterations <- 1
opts$lsi.cluster.resolution <- 2
opts$lsi.varFeatures <- 50000
opts$lsi.dims <- 50

opts$umap.seed <- 42
opts$umap.neighbours <- 30
opts$umap.minDist <- 0.5

###########################
## Latent Semantic Index ##
###########################

# Iterative LSI: two iterations
ArchRProject <- addIterativeLSI(ArchRProject,
  useMatrix = "PeakMatrix", 
  name = "IterativeLSI", 
  iterations = opts$lsi.iterations, 
  clusterParams = list(
    resolution = opts$lsi.cluster.resolution, 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = opts$lsi.varFeatures, 
  dimsToUse = 1:opts$lsi.dims,
  force = TRUE
)

# Iterative LSI: four iterations and start from a lower intial clustering resolution
# ArchRProject <- addIterativeLSI(ArchRProject,
#   useMatrix = "PeakMatrix", 
#   name = "IterativeLSI2", 
#   iterations = 4, 
#   clusterParams = list( #See Seurat::FindClusters
#     resolution = c(0.1, 0.2, 0.4), 
#     sampleCells = 10000, 
#     n.start = 10
#   ), 
#   varFeatures = 15000, 
#   dimsToUse = 1:30
# )


##########
## UMAP ##
##########

ArchRProject <- addUMAP(ArchRProject, 
  reducedDims = "IterativeLSI", 
  name = "UMAP",
  metric = "cosine",
  nNeighbors = opts$umap.neighbours, 
  minDist = opts$umap.minDist, 
  seed = opts$umap.seed,
  force = TRUE
)
# head(ArchRProject@embeddings$UMAP$df)

###########
## t-SNE ##
###########

# ArchRProject <- addTSNE(ArchRProject,
#   reducedDims = "IterativeLSI",
#   name = "TSNE",
#   perplexity = 30
# )
# head(ArchRProject@embeddings$TSNE$df)

# plotEmbedding(ArchRProject, colorBy = "cellColData", name = "Sample", embedding = "UMAP")

##########
## Plot ##
##########

# plotEmbedding(ArchRProject, 
#   colorBy = "cellColData", 
#   name = "celltype",
#   embedding = "TSNE",
#   size = 0.1,
#   keepAxis = F
# )

to.plot <- getEmbedding(ArchRProject,"UMAP") %>%
  as.data.table(keep.rownames = T) %>%
  setnames(c("rn","umap1","umap2")) %>%
  merge(getCellColData(ArchRProject) %>% as.data.table(keep.rownames = T), by="rn")

p <- ggplot(to.plot, aes(x=umap1, y=umap2)) +
  # geom_point(aes(color=celltype), size=1.5) +
  # scale_color_manual(values=opts$celltype.colors) +
  # guides(color = guide_legend(override.aes = list(size=4))) +
  geom_point(aes(fill=celltype), size=2, shape=21, color="black", stroke=0.1) +
  scale_fill_manual(values=opts$celltype.colors) +
  guides(fill = guide_legend(override.aes = list(size=4))) +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    legend.position = "right",
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank()
  )

##########
## Save ##
##########


# save umap plot
pdf(sprintf("%s/archr_umap_celltype.pdf",io$outdir), width=8.5, height=7)
print(p)
dev.off()

# save options
options.to.save <- opts[c("lsi.iterations", "lsi.cluster.resolution", "lsi.varFeatures", "lsi.dims", "umap.seed", "umap.neighbours", "umap.minDist")]
saveRDS(options.to.save, sprintf("%s/hyperparameters.rds",io$outdir))