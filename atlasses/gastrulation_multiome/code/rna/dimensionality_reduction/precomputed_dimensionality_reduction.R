source("/Users/ricard/gastrulation_multiome_10x/settings.R")

#####################
## I/O and options ##
#####################

io$outdir <- paste0(io$basedir,"/results/rna/dimensionality_reduction/pdf")

opts$aggregated.celltypes <- c(
  "Erythroid1" = "Erythroid",
  "Erythroid2" = "Erythroid",
  "Erythroid3" = "Erythroid",
  "Blood_progenitors_1" = "Blood_progenitors",
  "Blood_progenitors_2" = "Blood_progenitors",
  "Rostral_neurectoderm" = "Neurectoderm",
  "Caudal_neurectoderm" = "Neurectoderm",
  "Anterior_Primitive_Streak" = "Primitive_Streak"
)

###############
## Load data ##
###############

io$umap <- "/Users/ricard/data/gastrulation_multiome_10x/results/rna/dimensionality_reduction/umap.txt.gz"
umap.dt <- fread(io$umap)

to.plot <- sample_metadata %>%
  # Remove specific lineages
  # .[!celltype%in%c("ExE_ectoderm", "Parietal_endoderm","Visceral_endoderm","ExE_endoderm","PGC")] %>%
  # Aggregate celltypes
  .[,aggregated_celltype:=stringr::str_replace_all(celltype.mapped,opts$aggregated.celltypes)] %>%
  .[,aggregated_celltype:=factor(aggregated_celltype, levels=names(opts$celltype.colors))] %>%
  .[,aggregated_celltype:=stringr::str_replace_all(aggregated_celltype,"_"," ")]
  droplevels() %>%
  merge(umap.dt,by="cell")

################
## Parse data ##
################

  
names(opts$celltype.colors) <- names(opts$celltype.colors) %>% stringr::str_replace_all("_"," ")
opts$celltype.colors <- opts$celltype.colors[names(opts$celltype.colors) %in% unique(to.plot$aggregated_celltype)]

stopifnot(all(unique(to.plot$aggregated_celltype) %in% names(opts$celltype.colors)))
unique(to.plot$aggregated_celltype)[!unique(to.plot$aggregated_celltype) %in% names(opts$celltype.colors)]

##########
## Plot ##
##########

p <- ggplot(to.plot, aes(x=umapX, y=umapY)) +
  # geom_point(aes(colour=aggregated_celltype), size=0.1) +
  ggrastr::geom_point_rast(aes(colour=aggregated_celltype), size=0.05) +
  scale_color_manual(values=opts$celltype.colors) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    legend.position = "none",
    axis.text = element_blank(),
    axis.line = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank()
  )

pdf(paste0(io$outdir,"/umap.pdf"), width=4.5, height=4.5, useDingbats = F)
print(p)
dev.off()


# Plot each sample separately
for (i in unique(to.plot$sample)) {
  print(i)
  
  to.plot_i <- to.plot %>% copy %>%
    .[,alpha:=1.0] %>%
    .[sample!=i,c("celltype","alpha"):=list("None",0.25)]
  
  p <- ggplot() +
    ggrastr::geom_point_rast(aes(x=umapX, y=umapY), size=1.5, color="grey", alpha=0.25, data=to.plot_i[sample!=i]) +
    ggrastr::geom_point_rast(aes(x=umapX, y=umapY, fill=aggregated_celltype), size=2, shape=21, alpha=1.0, data=to.plot_i[sample==i]) +
    scale_fill_manual(values=opts$celltype.colors) +
    theme_classic() +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      legend.position="none"
    )
  
  pdf(sprintf("%s/umap_sample%s.pdf",io$outdir,i), width=5, height=5, useDingbats = F)
  print(p)
  dev.off()
}


