#####################
## Define settings ##
#####################

# load default setings
source(here::here("settings.R"))
source(here::here("utils.R"))

# Options
opts$celltypes = c(
  "Epiblast",
  "Primitive_Streak",
  "Caudal_epiblast",
  "PGC",
  "Anterior_Primitive_Streak",
  "Notochord",
  "Def._endoderm",
  "Gut",
  "Nascent_mesoderm",
  "Mixed_mesoderm",
  "Intermediate_mesoderm",
  "Caudal_Mesoderm",
  "Paraxial_mesoderm",
  "Somitic_mesoderm",
  "Pharyngeal_mesoderm",
  "Cardiomyocytes",
  "Allantois",
  "ExE_mesoderm",
  "Mesenchyme",
  "Haematoendothelial_progenitors",
  "Endothelium",
  "Blood_progenitors_1",
  "Blood_progenitors_2",
  "Erythroid1",
  "Erythroid2",
  "Erythroid3",
  "NMP",
  "Rostral_neurectoderm",
  "Caudal_neurectoderm",
  "Neural_crest",
  "Forebrain_Midbrain_Hindbrain",
  "Spinal_cord",
  "Surface_ectoderm",
  "Visceral_endoderm",
  "ExE_endoderm",
  "ExE_ectoderm",
  "Parietal_endoderm"
)

opts$motif_annotation <- "Motif_cisbp_lenient" # Motif_JASPAR2020_human

# I/O
# io$archR.pseudobulk.deviations.se <- sprintf("%s/pseudobulk/pseudobulk_DeviationMatrix_%s_summarized_experiment.rds",io$archR.directory,opts$motif_annotation)
# io$archR.pseudobulk.deviations.se <- sprintf("%s/results/atac/archR/chromvar/pseudobulk/chromVAR_deviations_summarized_experiment_%s_pseudobulk_correlated_peaks_archr.rds",io$basedir,opts$motif_annotation)
io$outdir <- paste0(io$basedir,"/results/rna_atac/rna_vs_chromvar/pseudobulk/per_celltype")

################################################
## Load pseudobulk RNA and chromVAR estimates ##
################################################

source(here::here("rna_atac/load_rna_atac_pseudobulk.R"))

#######################################################
## Load RNA vs chromVAR correlation results per gene ##
#######################################################

cor_rna_vs_chromvar_per_gene.dt <- fread(paste0(io$basedir,"/results/rna_atac/rna_vs_chromvar/pseudobulk/per_gene/cor_rna_vs_chromRNA_pseudobulk.txt.gz")) %>%
  .[,cor_sign:=as.factor(c("Repressor","Activator")[(r>0)+1])]

###########
## Merge ##
###########

rna_chromvar.dt <- merge(
  rna_tf_pseudobulk.dt,
  atac_chromvar_pseudobulk.dt,
  by = c("celltype","gene")
)

length(unique(rna_chromvar.dt$gene))

###############################
## Scatterplot per cell type ##
###############################

opts$max.chromvar <- 12
opts$max.expr <- 12

i <- "Neural_crest"
for (i in unique(rna_chromvar.dt$celltype)) {
  
  to.plot <- rna_chromvar.dt[celltype==i]  %>%
    merge(cor_rna_vs_chromvar_per_gene.dt[,c("gene","cor_sign")], by="gene") %>%
    .[chromvar_zscore<0,chromvar_zscore:=0] %>%
    .[chromvar_zscore>=opts$max.chromvar,chromvar_zscore:=opts$max.chromvar] %>%
    .[expr>=opts$max.expr,expr:=opts$max.expr]
  
  to.plot.text <- to.plot[expr>5 & chromvar_zscore>3.5] 
  to.plot.dots <- to.plot[!gene%in%to.plot.text$gene] 
  
  to.plot %>% .[,dot_size:=minmax.normalisation(expr)*minmax.normalisation(chromvar_zscore)]
  
  p <- ggplot(to.plot, aes(x=chromvar_zscore, y=expr)) +
    geom_point(aes(size=dot_size, fill=cor_sign), shape=21) + 
    # geom_point(size=0.5, data=to.plot.text) +
    ggrepel::geom_text_repel(data=to.plot.text, aes(label=gene), size=4) +
    # geom_text(aes(label=gene), size=3, data=to.plot.text) +
    scale_size_continuous(range = c(0.1,7)) +
    # scale_fill_gradient(low = "gray80", high = "darkgreen") +
    coord_cartesian(
      xlim = c(0,opts$max.chromvar+0.1),
      ylim = c(min(rna_chromvar.dt$expr),opts$max.expr+0.1)
      ) +
    labs(x="Motif accessibility (chromVAR+, z-score)", y="Gene expression") +
    guides(size="none", fill = guide_legend(override.aes = list(size=3))) +
    scale_fill_brewer(palette="Dark2") +
    theme_classic() +
    theme(
      legend.position = "none",
      legend.title = element_blank(),
      axis.text = element_text(size=rel(0.75), color="black")
    )
  
  pdf(sprintf("%s/scatterplots/%s.pdf",io$outdir,i), width = 6, height = 5)
  print(p)
  dev.off()
} 


#################################################
## Quantify number of active TFs per cell type ##
#################################################

to.plot <- rna_chromvar.dt %>%
  .[,sum(expr>5 & chromvar_zscore>15),by="celltype"] %>%
  # .[,sum(chromvar_zscore>25),by="celltype"] %>%
  setorder(-V1) %>% .[,cellype:=factor(celltype,levels=celltype)]


p <- ggbarplot(to.plot, x="celltype", y="V1", fill="celltype", sort.val = "asc") +
  scale_fill_manual(values=opts$celltype.colors) +
  labs(x="", y="Number of active TFs") +
  guides(x = guide_axis(angle = 90)) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(color="black", size=rel(0.8)),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank()
)

pdf(sprintf("%s/number_active_TFs.pdf",io$outdir), width = 9, height = 5)
print(p)
dev.off()


to.plot <- rna_chromvar.dt %>%
  merge(cor_rna_vs_chromvar_per_gene.dt[,c("gene","cor_sign")], by="gene") %>%
  .[,sum(expr>5 & chromvar_zscore>15),by=c("celltype","cor_sign")]# %>%
  # .[,celltype:=factor(celltype,levels=opts$celltypes)]# %>%
  # dcast(celltype~cor_sign, value.var="V1") %>%
  # .[,ratio_activators:=Activator/(Activator+Repressor)]

celltype.order <- to.plot %>%
  dcast(celltype~cor_sign, value.var="V1") %>%
  .[,ratio_activators:=Activator/(Activator+Repressor)] %>%
  setorder(ratio_activators) %>% .$celltype
to.plot[,celltype:=factor(celltype,levels=celltype.order)]

p <- ggbarplot(to.plot, x="celltype", y="V1", fill="cor_sign") +
  scale_fill_brewer(palette="Dark2") +
  # scale_fill_manual(values=opts$celltype.colors) +
  labs(x="", y="Number of active TFs") +
  guides(x = guide_axis(angle = 90)) +
  theme_classic() +
  theme(
    legend.position = "top",
    axis.text.x = element_text(color="black", size=rel(0.8)),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank()
  )
