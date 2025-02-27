library(VennDiagram)

source(here::here("settings.R"))

#####################
## Define settings ##
#####################

# I/O
io$outdir <- file.path(io$basedir,"results/rna_atac/qc"); dir.create(io$outdir, showWarnings = F)

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(io$metadata)

##################################################
## Scatterplot of RNA QC stats vs ATAC QC stats ##
##################################################

cols <- c("cell", "sample","nCount_RNA", "nFeature_RNA", "mitochondrial_percent_RNA", "TSSEnrichment_atac","ReadsInTSS_atac", "PromoterRatio_atac", "nFrags_atac", "BlacklistRatio_atac","pass_rnaQC","pass_atacQC")

to.plot <- sample_metadata[,..cols]

ggscatter(to.plot, x="mitochondrial_percent_RNA", y="nFeature_RNA", size=0.5, color="pass_rnaQC",
          add="reg.line", add.params = list(color="blue", fill="lightgray"), conf.int=TRUE) +
  yscale("log10", .format = TRUE) + 
  # xscale("log10", .format = TRUE) +
  facet_wrap(~sample, scales="fixed")
  # xscale("log2", .format = TRUE) +
  # yscale("log2", .format = TRUE)
  # theme(
  #   axis.text = element_text(size = rel(0.75), color="black")
  # )
  
# pdf(paste0(io$outdir,"/qc_rna.pdf"), width=11, height=5, useDingbats = F)
# print(p)
# dev.off()

################################################
## Venn Diagram of pass_rnaQC and pass_atacQC ##
################################################

p <- venn.diagram(
  x = list("RNA"=sample_metadata[pass_rnaQC==T,cell],
           "ATAC"=sample_metadata[pass_atacQC==T,cell]
  ),
  filename=NULL,
  col="transparent", fill=c("#3CB54E","#6691CB"), alpha = 0.60, cex = 1.5, 
  fontfamily = "serif", fontface = "bold")

pdf(file.path(io$outdir,"venn_rnaQC_vs_atacQC.pdf"))
grid.draw(p)
dev.off()