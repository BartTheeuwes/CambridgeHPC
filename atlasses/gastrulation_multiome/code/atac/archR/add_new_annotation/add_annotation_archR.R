library(GenomicRanges)

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

# io$metadata <- paste0(io$basedir,"/sample_metadata.txt.gz")
# io$outdir <- paste0(io$basedir,"/results/atac/archR/gene_scores")

#####################
## Load annotation ##
#####################

anno.gr <- fread("/Users/ricard/data/mm10_regulation/blood/Tal1_Chipseq/intersected_with_atac/D340004_Scl_intersected_with_atac.bed") %>%
  setnames(c("chr","start","end")) %>%
  makeGRangesFromDataFrame()

##################
## Add to ArchR ##
##################

# This function will add total counts of scATAC cells in provided features into ArchRProject.
foo <- addFeatureCounts(
  ArchRProject,
  features = anno.gr,
  name = "Tal1_ChIP_D340004_Scl",
  addRatio = FALSE
)

# This function for each sample will independently compute counts for each feature per cell in the provided ArchRProject or set of ArrowFiles.
addFeatureMatrix(
  input = ArchRProject,
  features = anno.gr,
  matrixName = "Tal1_ChIP_D340004_Scl",
  ceiling = 10^9,
  binarize = TRUE,
  verbose = TRUE
)

#############
## Explore ##
#############

sample_metadata <- fread(io$metadata) %>%
  .[pass_atacQC==TRUE & !is.na(celltype.mapped)]

acc.dt <- getCellColData(foo, c("Tal1_ChIP_D340004_SclCounts","Tal1_ChIP_D340004_SclRatio")) %>%
  as.data.table(keep.rownames = T) %>% setnames("rn","archR_cell") %>%
  merge(sample_metadata[,c("archR_cell","stage","celltype.mapped")], by="archR_cell")

to.plot <- acc.dt %>% 
  .[,N:=.N,by="celltype.mapped"] %>% 
  .[,sum(Tal1_ChIP_D340004_SclRatio>0.003)/unique(N),by="celltype.mapped"] %>%
  .[,celltype.mapped:=factor(celltype.mapped,levels=opts$celltypes)]

ggbarplot(to.plot, x="celltype.mapped", y="V1", fill="celltype.mapped", sort.val = "asc")   +
  scale_fill_manual(values=opts$celltype.colors) +
  labs(x="", y="Fraction of cells with Tal1 signal") +
  guides(x = guide_axis(angle = 90)) +
  theme(
    axis.text.x = element_text(color="black", size=rel(0.75)),
    # axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    # legend.title = element_blank(),
    legend.position = "none"
    # axis.title = element_blank(),
  )

# Correlate Tal1-ATAC signature with Tal1-RNA expression
cells <- fread(io$metadata) %>% .[pass_rnaQC==TRUE & !is.na(celltype.mapped),cell]
sce <- load_SingleCellExperiment(io$sce, normalise = TRUE, cells = cells)
colData(sce) <- sample_metadata %>% as.data.frame %>% tibble::column_to_rownames("cell") %>% .[colnames(sce),] %>% DataFrame()
rna.dt <- data.table(cell = colnames(sce), expr = logcounts(sce["Tal1",])[1,]) %>% merge(sample_metadata[,c("cell","archR_cell","stage","celltype.mapped")], by="cell")

to.plot <- merge(
  rna.dt[,c("archR_cell","expr")],
  acc.dt[,c("archR_cell", "Tal1_ChIP_D340004_SclRatio")],
  by="archR_cell"
) %>% merge(sample_metadata[,c("archR_cell","stage","celltype.mapped")])

to.plot[Tal1_ChIP_D340004_SclRatio>0.01,Tal1_ChIP_D340004_SclRatio:=0.01]

p <- ggscatter(to.plot, x="Tal1_ChIP_D340004_SclRatio", y="expr", size=0.5,
          add="reg.line", add.params = list(color="blue", fill="lightgray"), conf.int=TRUE) +
  facet_wrap(~celltype.mapped, scales="fixed") +
  labs(x="Tal1 accessibility signal", y="Tal1 RNA expression") +
  theme(
    axis.text.y = element_text(size=rel(0.75)),
    axis.text.x = element_blank()
  )

io$outdir <- "/Users/ricard/data/mm10_regulation/blood/Tal1_Chipseq/pdf"
pdf(sprintf("%s/Tal1_acc_vs_expr_per_cell.pdf",io$outdir), width=9, height=8)
print(p)
dev.off()

########################
## Save ARchR project ##
########################

# saveArchRProject(ArchRProject.filt)
