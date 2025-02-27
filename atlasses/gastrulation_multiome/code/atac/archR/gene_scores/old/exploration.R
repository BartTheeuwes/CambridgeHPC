
#################
## Exploration ##
#################

GeneScoreMatrix.se <- getMatrixFromProject(ArchRProject, useMatrix="GeneScoreMatrix")
GeneScoreMatrix.se <- getMatrixFromProject(ArchRProject, useMatrix="GeneScoreMatrix_nodistal")

# Pseudobulk
GeneScoreMatrix_pseudobulk.se <- getGroupSE(ArchRProject, groupBy = "celltype.mapped", useMatrix = "GeneScoreMatrix_nodistal")
dim(GeneScoreMatrix_pseudobulk.se)
rownames(GeneScoreMatrix_pseudobulk.se) <- rowData(GeneScoreMatrix_pseudobulk.se)$name

# Check genes
selected.genes <- rownames(GeneScoreMatrix_pseudobulk.se)[grep("Hbb",rownames(GeneScoreMatrix_pseudobulk.se))]
selected.genes <- rownames(GeneScoreMatrix_pseudobulk.se)[grep("Hoxc",rownames(GeneScoreMatrix_pseudobulk.se))]

to.plot <- assay(GeneScoreMatrix_pseudobulk.se[selected.genes,]) %>% as.matrix %>%
  as.data.table(keep.rownames = T) %>%
  setnames("rn","gene") %>%
  melt(id.vars="gene", variable.name="celltype") %>%
  .[,celltype:=factor(celltype,levels=opts$celltypes)]

ggbarplot(to.plot, x="celltype", y="value", fill="celltype") +
  facet_wrap(~gene) +
  scale_fill_manual(values=opts$celltype.colors) +
  labs(x="", y="Accessibility") +
  geom_hline(yintercept=0, colour="black", linetype="dashed") +
  guides(x = guide_axis(angle = 90)) +
  theme(
    # axis.text.x = element_text(size=rel(0.75)),
    axis.text.x = element_blank(),
    # axis.title.y = element_text(size=rel(0.75)),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    legend.title = element_blank()
  )
