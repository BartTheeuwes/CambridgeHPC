########################
## Load ArchR Project ##
########################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
} else {
  stop("Computer not recognised")
}

#####################
## Define settings ##
#####################

# I/O
# io$metadata <- paste0(io$basedir,"/processed/atac/archR/sample_metadata_after_archR.txt.gz")
io$outdir <- paste0(io$basedir,"/results/atac/archR/chromvar/comparison")

opts$motif_annotation <- "Motif_cisbp"

###########################################
## Fetch chromVAR deviations: approach A ##
###########################################

io$chromvar.se <- sprintf("%s/results/atac/archR/chromvar/deviations_summarized_experiment_%s.rds",io$basedir,opts$motif_annotation)

chromvar.deviations.se_A <- readRDS(io$chromvar.se)
chromvar.deviations.dt_A <- assay(chromvar.deviations.se_A,"z") %>% as.matrix %>% as.data.frame %>%
  as.data.table(keep.rownames = T) %>% setnames("rn","motif") %>%
  .[,motif:=stringr::str_split(motif,"_") %>% map_chr(1) %>% gsub("z:","",.)] %>%
  melt(id.vars="motif", variable.name="cell") %>%
  merge(sample_metadata[,c("cell","celltype.mapped")], by="cell") %>%
  setnames("celltype.mapped","celltype") %>%
  .[,class:=factor("A")]

###########################################
## Fetch chromVAR deviations: approach B ##
###########################################

io$chromvar.se <- sprintf("%s/results/atac/archR/chromvar/deviations_summarized_experiment_%s.rds",io$basedir,opts$motif_annotation)

chromvar.deviations.se_B <- readRDS(io$chromvar.se)
chromvar.deviations.dt_B <- assay(chromvar.deviations.se_B,"z") %>% as.matrix %>% as.data.frame %>%
  as.data.table(keep.rownames = T) %>% setnames("rn","motif") %>%
  .[,motif:=stringr::str_split(motif,"_") %>% map_chr(1) %>% gsub("z:","",.)] %>%
  melt(id.vars="motif", variable.name="cell") %>%
  merge(sample_metadata[,c("cell","celltype.mapped")], by="cell") %>%
  setnames("celltype.mapped","celltype") %>%
  .[,class:=factor("B")]

#############
## Combine ##
#############

chromvar.deviations.dt <- rbind(
  chromvar.deviations.dt_A,
  chromvar.deviations.dt_B
)

#############################################
## Boxplots of motif z-score per cell type ##
#############################################

opts$celltype.colors <- opts$celltype.colors[names(opts$celltype.colors)%in%unique(chromvar.deviations.dt$celltype)]

for (i in unique(chromvar.deviations.dt$motif)) {
  
  to.plot <- chromvar.deviations.dt[motif==i] %>% 
    .[,celltype:=factor(celltype,levels=names(opts$celltype.colors))]
  
  p <- ggboxplot(to.plot, x="celltype", y="value", fill="celltype", outlier.shape=NA) +
    scale_fill_manual(values=opts$celltype.colors) +
    facet_wrap(~class) +
    geom_hline(yintercept=0, linetype="dashed") +
    guides(x = guide_axis(angle = 90)) +
    labs(x="", y=sprintf("%s z-score",i)) +
    theme(
      legend.position = "none",
      # axis.text.x = element_text(color="black", angle=40, hjust=1, size=rel(0.75)),
      axis.text.x = element_text(color="black", size=rel(0.5)),
      axis.text.y = element_text(color="black", size=rel(0.8))
    )
  pdf(sprintf("%s/%s_chromvar_comparison.pdf",io$outdir,i), width=10, height=4)
  print(p)
  dev.off()
}
