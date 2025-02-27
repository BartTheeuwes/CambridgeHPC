library(ggpubr)

########################
## Load ArchR project ##
########################

source("/Users/ricard/gastrulation_multiome_10x/public_datasets/Pijuan-Sala_2020/archR/load_archR_project.R")

names(ArchRProject@embeddings)

#####################
## Define settings ##
#####################

io$outdir <- paste0(io$basedir,"/results/celltype_annotation")

##############################
## Create Gene Score matrix ##
#############################


# addGeneScoreMatrix(ArchRProject,
#   genes = getGenes(ArchRProject),
#   matrixName = "GeneScoreMatrix",
#   useTSS = FALSE,
#   extendTSS = FALSE
# )

###########################
## Get Gene Score matrix ##
###########################

# properties of the AchR model (https://www.archrproject.com/bookdown/calculating-gene-scores-in-archr.html): 
# - Accessibility within the entire gene body contributes to the gene score.
# - An exponential weighting function that accounts for the activity of putative distal regulatory elements in a distance-dependent fashion.
# - Imposed gene boundaries that minimizes the contribution of unrelated regulatory elements to the gene score.

# returns a SummarizedExperiment class
gene.score.matrix <- getMatrixFromProject(ArchRProject, useMatrix = "GeneScoreMatrix")
assay(gene.score.matrix)[1:3,1:3]
head(colData(gene.score.matrix))
head(rowData(gene.score.matrix))

hist(as.numeric(assay(gene.score.matrix[1:1000,1:5000])), breaks=30)

########################################
## Load marker genes based on the RNA ##
########################################

rna.marker_genes <- fread("/Users/ricard/data/gastrulation10x/results/marker_genes/E8.5/marker_genes.txt.gz") %>%
  .[gene%in%rowData(gene.score.matrix)[["name"]]]

# rna.marker_genes[,length(unique(gene)),by="celltype"]

################
## Parse data ##
################

gene.score.matrix <- gene.score.matrix %>%
  .[rowData(.)[["name"]] %in% unique(rna.marker_genes$gene),]


##################################################
## Calculate cell type affinities for each cell ##
##################################################

dt <- head(colnames(gene.score.matrix),n=100) %>% map(function(i) {
  data.table(value=assay(gene.score.matrix[,i])[,1]) %>%
  .[,gene:=rowData(gene.score.matrix)[["name"]]] %>%
    merge(rna.marker_genes[,c("celltype","gene")], by="gene", allow.cartesian=TRUE) %>%
    .[,.(score=mean(value)),by="celltype"] %>%
    .[,cell:=gsub("E8.5_gastrulation#","",i)]
}) %>% rbindlist %>% 
  merge(sample_metadata[,c("cell","celltype")] %>% setnames("celltype","true_celltype"), by="cell")

##########
## Plot ##
##########

for (i in unique(dt$cell)) {
  
  to.plot <- dt[cell==i] %>%
    # .[,celltype:=stringr::str_replace_all(celltype,"_"," ")] %>%
    .[,celltype:=factor(celltype,levels=names(opts$rna.celltype.colors))]
  
  p <- ggbarplot(to.plot, x="celltype", y="score", fill="celltype") +
    ggtitle(unique(to.plot$true_celltype)) +
    scale_fill_manual(values=opts$rna.celltype.colors) +
    labs(x="", y="Cell type affinity") +
    theme(
      plot.title = element_text(hjust = 0.5, size=rel(0.9)),
      axis.text.y = element_text(size=rel(0.75)),
      axis.text.x = element_text(colour="black",size=rel(0.7), angle=90, hjust=1, vjust=0.5),
      axis.ticks.x = element_blank(),
      legend.position = "none"
    )
  
  pdf(sprintf("%s/%s_affinity.pdf",io$outdir,i), width = 9, height = 5)
  print(p)
  dev.off()
}










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
