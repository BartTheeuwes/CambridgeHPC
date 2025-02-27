library(chromVAR)

#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
} else {
  stop("Computer not recognised")
}

# Options
opts$motif_annotation <- "Motif_cisbp"

# I/O
io$chromvar.se.A <- sprintf("%s/results/atac/archR/chromvar/pseudobulk/chromVAR_deviations_summarized_experiment_%s_pseudobulk.rds",io$basedir,opts$motif_annotation)
io$chromvar.se.B <- sprintf("%s/processed/atac/archR_subset/pseudobulk/pseudobulk_DeviationMatrix_%s_summarized_experiment.rds",io$basedir,opts$motif_annotation)
io$outdir <- paste0(io$basedir,"/results/atac/archR/chromvar/pseudobulk/comparison")

###########################################
## Fetch chromVAR deviations: approach A ##
###########################################

# Load
chromvar.deviations.se_A <- readRDS(io$chromvar.se.A)

# Rename genes
# rownames(chromvar.deviations.se_A) <- rowData(chromvar.deviations.se_A)$name %>% toupper %>% stringr::str_split(.,"_") %>% map_chr(1)

# Remove duplicated motifs
chromvar.deviations.se_A <- chromvar.deviations.se_A[!duplicated(rownames(chromvar.deviations.se_A)),]

# create data.table
chromvar.deviations.dt_A <- assay(chromvar.deviations.se_A,"z") %>% as.matrix %>% as.data.frame %>%
  as.data.table(keep.rownames = T) %>% setnames("rn","motif") %>%
  melt(id.vars="motif", variable.name="celltype") %>%
  .[,class:=factor("pseudobulk")]

###########################################
## Fetch chromVAR deviations: approach B ##
###########################################

# Load
chromvar.deviations.se_B <- readRDS(io$chromvar.se.B)
chromvar.deviations.se_B <- chromvar.deviations.se_B[rowData(chromvar.deviations.se_B)$seqnames=="z",]

# Rename genes
rownames(chromvar.deviations.se_B) <- rowData(chromvar.deviations.se_B)$name %>% toupper %>% stringr::str_split(.,"_") %>% map_chr(1)

# Remove duplicated motifs
chromvar.deviations.se_B <- chromvar.deviations.se_B[!duplicated(rownames(chromvar.deviations.se_B)),]

# create data.table
chromvar.deviations.dt_B <- assay(chromvar.deviations.se_B) %>% as.matrix %>% as.data.frame %>%
  as.data.table(keep.rownames = T) %>% setnames("rn","motif") %>%
  melt(id.vars="motif", variable.name="celltype") %>%
  .[,class:=factor("single-cell")]

#############
## Combine ##
#############

motifs <- intersect(unique(chromvar.deviations.dt_A$motif),unique(chromvar.deviations.dt_B$motif))

chromvar.deviations.dt <- rbind(
  chromvar.deviations.dt_A[motif%in%motifs],
  chromvar.deviations.dt_B[motif%in%motifs]
)
length(unique(chromvar.deviations.dt$motif))

#############################################
## Boxplots of motif z-score per cell type ##
#############################################

opts$celltype.colors <- opts$celltype.colors[names(opts$celltype.colors)%in%unique(chromvar.deviations.dt$celltype)]

for (i in unique(chromvar.deviations.dt$motif)) {
  
  to.plot <- chromvar.deviations.dt[motif==i] %>% 
    .[,celltype:=factor(celltype,levels=names(opts$celltype.colors))]
  
  p1 <- ggbarplot(to.plot, x="celltype", y="value", fill="celltype") +
    scale_fill_manual(values=opts$celltype.colors) +
    facet_wrap(~class) +
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
  p2 <- ggscatter(to.plot2, x="single-cell", y="pseudobulk", fill="celltype", shape=21, stroke=0.1, size=3.5, 
                  add="reg.line", add.params = list(color="black", fill="lightgray"), conf.int=TRUE) +
    stat_cor(method = "pearson", label.x.npc = "middle", label.y.npc = "bottom") +
    scale_fill_manual(values=opts$celltype.colors) +
    labs(x="chromVAR z-score (single-cell)", y="chromVAR z-score (pseudobulk)") +
    theme(
      legend.position = "none",
      axis.text = element_text(size=rel(0.7)),
      axis.title = element_text(size=rel(0.85))
    )
  
  p <- cowplot::plot_grid(plotlist=list(p1,p2), nrow=1, rel_widths = c(2/3,1/3))
  
  pdf(sprintf("%s/%s_chromvar_singlecell_vs_pseudobulk.pdf",io$outdir,i), width=13, height=5)
  print(p)
  dev.off()
}
