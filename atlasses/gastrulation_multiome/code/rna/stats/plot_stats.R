library(ggpubr)

#####################
## Define settings ##
#####################

source(here::here("settings.R"))
# io$metadata <- paste0(io$basedir,"/results/rna/mapping/sample_metadata_after_mapping.txt.gz")
io$outdir <- paste0(io$basedir,"/results/rna/general_stats")

sample_metadata <- fread(io$metadata) %>% 
  .[sample%in%opts$samples]

# sample_metadata <- fread("/Users/ricard/data/gastrulation_multiome_10x/processed/first_sample/sample_metadata.txt.gz")

###############################################
## Boxplots of general statistics per sample ##
###############################################

# to.plot <- sample_metadata %>% 
#   melt(id.vars=c("cell","sample","celltype"), measure.vars=c("nCount_RNA","nFeature_RNA"))

to.plot <- sample_metadata %>% 
  .[pass_QC==TRUE & !is.na(celltype.mapped)] %>% 
  .[,celltype.mapped:=factor(celltype.mapped, levels=names(opts$celltype.colors))]

p <- ggboxplot(to.plot, x="celltype.mapped", y="nCount_RNA", fill="celltype.mapped", outlier.shape=NA) +
  # yscale("log10", .format = TRUE) +
  labs(x="") +
  facet_wrap(~sample, scales="free_y", nrow=2) +
  coord_cartesian(ylim=c(0,3.5e4)) +
  scale_fill_manual(values=opts$celltype.colors) +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    axis.text.y = element_text(size=rel(0.75)),
    axis.text.x = element_text(colour="black",size=rel(0.8), angle=90, hjust=1, vjust=0.5),
    # axis.ticks.x = element_blank()
  )

pdf(paste0(io$outdir,"/general_stats_per_celltype.pdf"), width=8, height=6, useDingbats = F)
print(p)
dev.off()

########################################
## Barplots number of cells per sample ##
########################################

to.plot <- sample_metadata[,.N,by="sample"]

p <- ggbarplot(to.plot, x = "sample", y = "N", fill="gray70") +
  labs(x="", y="Number of cells (after QC)") +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    # axis.text.x = element_blank(),
    axis.text.x = element_text(colour="black",size=rel(0.8), angle=40, hjust=1),
    # axis.ticks.x = element_blank()
  )

pdf(paste0(io$outdir,"/N_per_sample.pdf"), width=8, height=6, useDingbats = F)
print(p)
dev.off()


##################################################
## Boxplots of general statistics per cell type ##
##################################################

to.plot <- sample_metadata %>% 
  melt(id.vars=c("cell","celltype.mapped"), measure.vars=c("nCount_RNA","nFeature_RNA"))

p <- ggboxplot(to.plot, x = "celltype.mapped", y = "value", fill="celltype.mapped", outlier.shape=NA) +
  yscale("log10", .format = TRUE) +
  labs(x="", y="") +
  scale_fill_manual(values=opts$celltype.colors) +
  facet_wrap(~variable, scales="free_y") +
  guides(fill = guide_legend(override.aes = list(size=0.25), ncol=1)) +
  scale_size(guide = 'none') +
  theme(
    legend.position = "right",
    legend.text = element_text(size=rel(0.75)),
    legend.title = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

pdf(paste0(io$outdir,"/general_stats_per_celltype.pdf"), width=16, height=10, useDingbats = F)
print(p)
dev.off()

