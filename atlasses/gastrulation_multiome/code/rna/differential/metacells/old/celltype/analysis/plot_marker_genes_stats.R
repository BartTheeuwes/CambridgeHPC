################################################
## Plot number of marker genes per cell types ##
################################################

to.plot <- markers_genes.dt %>% .[,.N,by=c("celltype")]

p <- ggbarplot(to.plot, x="celltype", y="N", fill="celltype") +
  scale_fill_manual(values=opts$celltype.colors) +
  labs(x="", y="Number of marker genes") +
  theme(
    axis.text.y = element_text(size=rel(0.65)),
    axis.text.x = element_text(colour="black",size=rel(0.7), angle=90, hjust=1, vjust=0.5),
    axis.title = element_text(colour="black",size=rel(0.75)),
    axis.ticks.x = element_blank(),
    legend.position = "none"
)

pdf(file.path(io$outdir,"barplot_number_marker_genes.pdf"), width = 8, height = 4)
print(p)
dev.off()

##################################
## Plot gene marker exclusivity ##
##################################

to.plot <- markers_genes.dt %>%
  .[,.(Nx=.N),by="gene"] %>%
  .[,Nx:=factor(Nx)] %>%
  .[,.(Ny=.N),by="Nx"]

p <- ggbarplot(to.plot, x="Nx", y="Ny", fill="gray70") +
  labs(x="Number of different cell types per marker gene", y="") +
  theme(
    axis.text = element_text(size=rel(0.75)),
  )
pdf(file.path(io$outdir,"boxplot_exclusivity_per_gene.pdf"), width = 7, height = 5)
print(p)
dev.off()

################################################
## Plot gene marker exclusivity per cell type ##
################################################

to.plot <- markers_genes.dt %>% .[,N:=.N,by="gene"]

p <- ggboxplot(to.plot, x="celltype", y="N", fill="celltype", color="black") +
  scale_fill_manual(values=opts$celltype.colors) +
  labs(x="", y="Exclusivity of gene markers\n(the smaller the more exclusive)") +
  theme(
    axis.text.y = element_text(size=rel(0.75)),
    axis.title.y = element_text(size=rel(0.85)),
    axis.text.x = element_text(colour="black",size=rel(0.7), angle=90, hjust=1, vjust=0.5),
    legend.position = "none"
  )

pdf(file.path(io$outdir,"boxplot_exclusivity_per_celltype.pdf"), width = 9, height = 5)
print(p)
dev.off()


