library(ggpubr)

########################
## Load ArchR Project ##
########################

if (grepl("ricard",Sys.info()['nodename'])) {
  # source("/Users/ricard/gastrulation_multiome_10x/atac/archR/load_archR_project.R")
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  # source("/homes/ricard/gastrulation_multiome_10x/atac/archR/load_archR_project.R")
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
} else {
  stop("Computer not recognised")
}

#####################
## Define settings ##
#####################

# I/O
# io$metadata <- paste0(io$basedir,"/processed/atac/archR/sample_metadata_after_archR.txt.gz")
io$chromvar.dir <- paste0(io$basedir,"/results/atac/archR/chromvar/pseudobulk")
io$outdir <- paste0(io$basedir,"/results/atac/archR/chromvar/pseudobulk/pdf")

# Options
opts$motif.annotation <- "Motif_JASPAR2020_human"

######################################
## Load pseudobulk chromVAR results ##
######################################

chromvar.deviations.se <- readRDS(sprintf("%s/deviations_summarized_experiment_%s.rds",io$chromvar.dir,opts$motif.annotation))

#############################################
## Barplots of motif z-score per cell type ##
#############################################

chromvar.dt <- assay(chromvar.deviations.se,"z") %>% as.matrix %>% as.data.frame %>%
  as.data.table(keep.rownames = T) %>% setnames("rn","motif") %>%
  .[,motif:=stringr::str_split(motif,"_") %>% map_chr(1)] %>%
  melt(id.vars="motif", variable.name="celltype")

opts$celltype.colors <- opts$celltype.colors[names(opts$celltype.colors)%in%unique(chromvar.dt$celltype)]

max.value <- 5
min.value <- -4
chromvar.dt[value>=max.value,value:=max.value]
chromvar.dt[value<=min.value,value:=min.value]

for (i in unique(chromvar.dt$motif)) {
  to.plot <- chromvar.dt[motif==i] %>% .[,celltype:=factor(celltype,levels=names(opts$celltype.colors))]
  
  p <- ggbarplot(to.plot, x="celltype", y="value", fill="celltype") +
    coord_cartesian(ylim=c(min.value,max.value)) +
    scale_fill_manual(values=opts$celltype.colors) +
    geom_hline(yintercept=0, linetype="dashed") +
    labs(x="", y=sprintf("%s z-score",i)) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(color="black", angle=40, hjust=1, size=rel(0.75)),
      axis.text.y = element_text(color="black", size=rel(0.8))
    )
  pdf(sprintf("%s/%s_chromvar_pseudobulk.pdf",io$outdir,i), width=7, height=4)
  print(p)
  dev.off()
}

