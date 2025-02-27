here::i_am("atac/archR/chromvar/pseudobulk/run_chromvar_pseudobulk.R")

source(here::here("settings.R"))

suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(chromVAR))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--motif_annotation',             type="character",            help='Motif annotation')
p$add_argument('--indir',          type="character",                               help='Input directory')
p$add_argument('--outdir',          type="character",                               help='Output directory')

args <- p$parse_args(commandArgs(TRUE))


#####################
## Define settings ##
#####################

## START TEST ##
args$motif_annotation <- "CISBP"
args$indir <- file.path(io$basedir,"results/atac/archR/chromvar/pseudobulk")
args$outdir <- file.path(io$basedir,"results/atac/archR/chromvar/pseudobulk/archr_vs_chromvar_comparison")
## END TEST ##

dir.create(args$outdir, showWarnings=F)

###########################################
## Fetch chromVAR deviations: approach A ##
###########################################

# Load
io$chromvar.se.A <- file.path(args$indir,sprintf("chromVAR_deviations_%s_pseudobulk_chromvar.rds",args$motif_annotation))
chromvar.deviations.se_A <- readRDS(io$chromvar.se.A)

# Remove duplicated motifs
chromvar.deviations.se_A <- chromvar.deviations.se_A[!duplicated(rownames(chromvar.deviations.se_A)),]

# create data.table
chromvar.deviations.dt_A <- assay(chromvar.deviations.se_A,"z") %>% as.matrix %>% as.data.frame %>%
  as.data.table(keep.rownames = T) %>% setnames("rn","motif") %>%
  melt(id.vars="motif", variable.name="celltype") %>%
  .[,class:=factor("chromVAR")]

###########################################
## Fetch chromVAR deviations: approach B ##
###########################################

# Load
io$chromvar.se.B <- file.path(args$indir,sprintf("chromVAR_deviations_%s_pseudobulk_archr.rds",args$motif_annotation))
chromvar.deviations.se_B <- readRDS(io$chromvar.se.B)

# Remove duplicated motifs
chromvar.deviations.se_B <- chromvar.deviations.se_B[!duplicated(rownames(chromvar.deviations.se_B)),]

# create data.table
chromvar.deviations.dt_B <- assay(chromvar.deviations.se_B,"z") %>% as.matrix %>% as.data.frame %>%
  as.data.table(keep.rownames = T) %>% setnames("rn","motif") %>%
  melt(id.vars="motif", variable.name="celltype") %>%
  .[,class:=factor("ArchR")]

#############
## Combine ##
#############

motifs <- intersect(unique(chromvar.deviations.dt_A$motif),unique(chromvar.deviations.dt_B$motif))

chromvar.deviations.dt <- rbind(
  chromvar.deviations.dt_A[motif%in%motifs],
  chromvar.deviations.dt_B[motif%in%motifs]
)
length(unique(chromvar.deviations.dt$motif))

###################
## Rename motifs ##
###################

opts$motif_annotation <- args$motif_annotation
source(here::here("atac/archR/load_motif_annotation.R"))

mean(chromvar.deviations.dt$motif %in% motif2gene.dt$motif)

chromvar.deviations.dt <- chromvar.deviations.dt %>% 
  merge(motif2gene.dt[,c("motif","gene")],by="motif")

#############################################
## Boxplots of motif z-score per cell type ##
#############################################

opts$celltype.colors <- opts$celltype.colors[names(opts$celltype.colors)%in%unique(chromvar.deviations.dt$celltype)]

# i <- "FOXA2"
# TFs.to.plot <- unique(chromvar.deviations.dt$motif)# %>% head(n=3)
TFs.to.plot <- fread(io$rna.atlas.marker_TFs.up)[["gene"]] %>% unique %>% toupper
TFs.to.plot <- TFs.to.plot[TFs.to.plot%in%chromvar.deviations.dt$gene] %>% head(n=50)

for (i in TFs.to.plot) {
  
  to.plot <- chromvar.deviations.dt[gene==i] %>% 
    .[,celltype:=factor(celltype,levels=names(opts$celltype.colors))]
  
  p1 <- ggbarplot(to.plot, x="celltype", y="value", fill="celltype") +
    scale_fill_manual(values=opts$celltype.colors) +
    facet_wrap(~class, scales = "free_y") +
    geom_hline(yintercept=0, linetype="dashed") +
    guides(x = guide_axis(angle = 90)) +
    labs(x="", y=sprintf("%s chromVAR z-score",i)) +
    theme(
      legend.position = "none",
      # axis.text.x = element_text(color="black", angle=40, hjust=1, size=rel(0.75)),
      axis.text = element_text(color="black", size=rel(0.6)),
      axis.title = element_text(color="black", size=rel(0.8))
    )
  
  to.plot2 <- to.plot %>% dcast(motif+celltype~class, value.var="value")
  p2 <- ggscatter(to.plot2, x="chromVAR", y="ArchR", fill="celltype", shape=21, stroke=0.1, size=3.5, 
                  add="reg.line", add.params = list(color="black", fill="lightgray"), conf.int=TRUE) +
    stat_cor(method = "pearson", label.x.npc = "middle", label.y.npc = "bottom") +
    scale_fill_manual(values=opts$celltype.colors) +
    labs(x="", y="") +
    theme(
      legend.position = "none",
      axis.text = element_text(size=rel(0.7)),
      axis.title = element_text(size=rel(0.85))
    )
  
  p <- cowplot::plot_grid(plotlist=list(p1,p2), nrow=1, rel_widths = c(2/3,1/3))
  
  pdf(file.path(args$outdir,sprintf("%s_%s_chromvar_comparison_pseudobulk.pdf",i,args$motif_annotation)), width=11, height=5)
  print(p)
  dev.off()
}

# Completion token
file.create(file.path(args$outdir,"completed.txt"))
