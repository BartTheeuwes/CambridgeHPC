#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
} else {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
}

io$metadata <- paste0(io$basedir,"/results/rna/doublets/sample_metadata_after_doublets.txt.gz")
io$outdir <- paste0(io$basedir,"/results/rna/celltype_proportions")

# opts$aggregate.celltypes <- c(
#   "Erythroid1" = "Erythroid",
#   "Erythroid2" = "Erythroid",
#   "Erythroid3" = "Erythroid",
#   "Blood_progenitors_1" = "Blood_progenitors",
#   "Blood_progenitors_2" = "Blood_progenitors",
#   "Rostral_neurectoderm" = "Neurectoderm",
#   "Caudal_neurectoderm" = "Neurectoderm",
#   "Anterior_Primitive_Streak" = "Primitive_Streak"
# )

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(io$metadata) %>%
  # .[pass_rnaQC==TRUE & doublet_call==FALSE & !is.na(celltype.mapped)]
  .[pass_QC==TRUE & hybrid_call==FALSE & !is.na(celltype.mapped)]

################################################
## Calculate cell type proportions per sample ##
################################################

to.plot <- sample_metadata %>%
  .[,N:=.N,by="sample"] %>%
  # .[,celltype.mapped:=stringr::str_replace_all(celltype.mapped,opts$aggregate.celltypes)] %>%
  .[,.(N=.N, celltype_proportion=.N/unique(N)),by=c("sample","stage","celltype.mapped")] %>%
  droplevels() %>% setorder(sample)  %>% .[,sample:=as.factor(sample)]

######################
## Stacked barplots ##
######################

p <- ggplot(to.plot, aes(x=sample, y=celltype_proportion)) +
  geom_bar(aes(fill=celltype.mapped), stat="identity", color="black") +
  facet_wrap(~stage, scales = "free_x", nrow=1) +
  scale_fill_manual(values=opts$celltype.colors) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank()
  )

pdf(sprintf("%s/celltype_proportions_stacked_barplots.pdf",io$outdir), width=7, height=5)
print(p)
dev.off()

#########################
## Horizontal barplots ##
#########################

# Rename "_" to " " in cell types
# to.plot[,celltype.mapped:=stringr::str_replace_all(celltype.mapped,"_"," ")]
# names(opts$celltype.colors) <- names(opts$celltype.colors) %>% stringr::str_replace_all("_"," ")

# Define colours and cell type order
opts$celltype.colors <- opts$celltype.colors[names(opts$celltype.colors) %in% unique(to.plot$celltype.mapped)]
to.plot[,celltype.mapped:=factor(celltype.mapped, levels=names(opts$celltype.colors))]

for (i in unique(to.plot$stage)) {
  p <- ggplot(to.plot[stage==i], aes(x=celltype.mapped, y=N)) +
    geom_bar(aes(fill=celltype.mapped), stat="identity", color="black") +
    scale_fill_manual(values=opts$celltype.colors) +
    facet_wrap(~sample, nrow=1, scales="fixed") +
    coord_flip() +
    labs(y="Number of cells") +
    theme_bw() +
    theme(
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_text(color="black", size=rel(1.2)),
      axis.title.x = element_text(color="black", size=rel(1)),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size=rel(1), color="black"),
      axis.text.x = element_text(size=rel(1), color="black")
    )
  
  pdf(sprintf("%s/celltype_proportions_%s_horizontal_barplots.pdf",io$outdir,i), width=10, height=5)
  print(p)
  dev.off()
}
