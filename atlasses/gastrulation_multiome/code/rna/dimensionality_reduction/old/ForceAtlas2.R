library(igraph) 
library(ForceAtlas2) 

graph_k10 <- buildKNNGraph(sce_filt, k = 5, use.dimred = "PCA")
# graph_k10 <- buildSNNGraph(sce_filt, k = 10, use.dimred = "PCA", type = "rank")

layout <- layout.forceatlas2(graph_k10, iterations=1500)
# rownames(layout) <- colnames(sce_filt)
# reducedDim(sce_filt, "ForceAtlas2") <- layout

to.plot <- layout %>% as.data.table %>% 
  .[,cell:=colnames(sce_filt)] %>%
  merge(sample_metadata, by="cell")

p <- ggplot(to.plot, aes(x=V1, y=V2, fill=celltype.mapped)) +
  geom_point(size=3, shape=21, stroke=0.2) +
  scale_fill_manual(values=opts$celltype.colors) +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position="right",
    legend.title=element_blank()
  )

# pdf(paste0(io$outdir,"/umap_by_celltype.pdf"), width=6, height=4, useDingbats = F)
print(p)
# dev.off()