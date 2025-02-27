##########################################
## Dimensionality reduction with Seurat ##
##########################################

# Run term frequency inverse document frequency (TF-IDF) normalization
srat <- RunTFIDF(srat)

max(srat@assays$peaks@counts)
max(srat@assays$peaks@data)

# Feature selection
#   This can be a percentile specified as 'q' followed by the minimum percentile, for example 'q5' to set the top 95% most common features. 
#   Alternatively, this can be an integer specifying the minumum number of cells containing the feature for the feature to be included in the set of VariableFeatures. 
srat <- FindTopFeatures(srat, min.cutoff = 'q5')
length(srat@assays$peaks@var.features)

# Run SVD
srat <- RunSVD(
  object = srat,
  assay = 'peaks',
  reduction.key = 'LSI_',
  reduction.name = 'lsi'
)

# Run non-linear dimensionality reduction
srat <- RunTSNE(srat, reduction = 'lsi', dims = 1:15)
srat <- RunUMAP(srat, reduction = 'lsi', dims = 1:15)

# Plot
DimPlot(srat, reduction="lsi", dims=c(1,2))
DimPlot(srat, reduction="tsne", dims=c(1,2))
DimPlot(srat, reduction="umap", dims=c(1,2))

##############################################
## Dimensionality reduction with cis-topics ##
##############################################

cisTopicObject<- createcisTopicObject(m, project.name='CellLineMixture')

# Define number of topics
topic<- c(15, 20, 25)

# Run Latent Dirichlet Allocation
cisTopicObject<- cisTopic::runCGSModels(cisTopicObject, topic, burnin = 250, iterations = 500, nCores = 2)

# Run Latent Dirichlet Allocation with WarpLDA
# cisTopicObject <- runWarpLDAModels(cisTopicObject, topic, nCores=2, iterations = 500, seed = 42)

# Model selection
cisTopicObject <- selectModel(cisTopicObject)
# par(mfrow=c(3,3))
# cisTopicObject <- selectModel(cisTopicObject, type='maximum')
# cisTopicObject <- selectModel(cisTopicObject, type='perplexity')
# cisTopicObject <- selectModel(cisTopicObject, type='derivative')
# dev.off()

# cisTopicObject@selected.model <- cisTopicObject@models$`30`
# corrplot::corrplot(cor(t(cisTopicObject@selected.model$document_expects)))
# cis_topics.red <- modelMatSelection(cisTopicObject, 'cell', 'Probability')[1:5,1:5]
# cis_topics.red <- cisTopicObject@selected.model$document_expects[1:5,1:5]

# run PCA on topics
# cisTopicObject <- cisTopic::runPCA(cisTopicObject, target="cell")

# run t-SNE/UMAP on topics
cisTopicObject <- runtSNE(cisTopicObject, target='cell', seed=123, pca=TRUE, method='Z-score')
cisTopicObject <- runUmap(cisTopicObject, target='cell', seed=123, method='Z-score')

# cisTopicObject@dr$cell$Umap
plotFeatures(cisTopicObject, method='tSNE', target='cell', topic_contr="Z-score", topics=2, 
             colorBy=NULL, cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE)
# plotFeatures(cisTopicObject, method='Umap', target='cell', topic_contr="Z-score", topics=1)

# # Extract cell and region embeddings
# embeddings <- cisTopicObject@dr$cell$PCA$ind.coord
# colnames(embeddings) <- paste0("topic_",1:ncol(embeddings))
# loadings <- t( cisTopicObject@selected.model$topics )
# 
# # Create a DimReduc object for Seurat
# foo <- CreateDimReducObject(
#   embeddings = embeddings,
#   loadings = loadings,
#   key = "topic_",
#   assay = "peaks"
# )
# 
# srat@reductions$cis_topics <- foo
# DimPlot(srat, reduction="cis_topics", dims=c(1,2))

##########
## Save ##
##########

saveRDS(srat, io$seurat.output)
