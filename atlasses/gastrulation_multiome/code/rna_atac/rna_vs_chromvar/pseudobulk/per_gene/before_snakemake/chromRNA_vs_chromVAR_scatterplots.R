#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
  source("/Users/ricard/gastrulation_multiome_10x/utils.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
  source("/homes/ricard/gastrulation_multiome_10x/utils.R")
} else {
  stop("Computer not recognised")
}

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

opts$motif_annotation <- "Motif_cisbp"

# I/O
io$outdir <- paste0(io$basedir,"/results/rna_atac/rna_vs_chromvar/pseudobulk/per_gene/chromRNA_vs_chromVAR"); dir.create(io$outdir)

#######################################
## Load pseudobulk RNA and ATAC data ##
#######################################

io$archR.pseudobulk.deviations.se <- sprintf("%s/results/atac/archR/chromvar/pseudobulk/chromVAR_deviations_%s_pseudobulk_archr.rds",io$basedir,opts$motif_annotation)
if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/rna_atac/load_rna_atac_pseudobulk.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/rna_atac/load_rna_atac_pseudobulk.R")
} else {
  stop("Computer not recognised")
}
atac.chromVAR.se <- atac.chromvar.se
atac_chromVAR.dt <- chromvar.dt %>% setnames("chromvar_zscore","chromVAR")


io$archR.pseudobulk.deviations.se <- sprintf("%s/results/atac/archR/chromvar/pseudobulk/chromRNA_deviations_%s_pseudobulk_archr.rds",io$basedir,opts$motif_annotation)
if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/rna_atac/load_rna_atac_pseudobulk.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/rna_atac/load_rna_atac_pseudobulk.R")
} else {
  stop("Computer not recognised")
}
atac_chromRNA.se <- atac.chromvar.se
atac_chromRNA.dt <- chromvar.dt %>% setnames("chromvar_zscore","chromRNA")

rm(list=c("atac.dt","atac.peakMatrix.se","atac.chromvar.se","chromvar.dt"))

#############################################
## Load pre-computed correlation estimates ##
#############################################

cor_rna_vs_chromVAR.dt <- fread(paste0(io$basedir,"/results/rna_atac/rna_vs_chromvar/pseudobulk/per_gene/cor_rna_vs_chromVAR_pseudobulk.txt.gz"))
cor_rna_vs_chromRNA.dt <- fread(paste0(io$basedir,"/results/rna_atac/rna_vs_chromvar/pseudobulk/per_gene/cor_rna_vs_chromRNA_pseudobulk.txt.gz"))

###########
## Merge ##
###########

rna_chrom.dt <- merge(
  rna_tf.dt,
  atac_chromVAR.dt,
  by = c("celltype","gene")
) %>% merge(atac_chromRNA.dt, by = c("celltype","gene"))


cor.dt <- merge(cor_rna_vs_chromvar.dt[,c("gene","r")], cor_rna_vs_chromrna.dt[,c("gene","r")], by=c("gene")) %>% melt(id.vars="gene")


######################################
## Scatter plot of individual genes ##
######################################

genes.to.plot <- unique(cor.dt$gene)
# genes.to.plot <- cor.dt[sig==T & abs(r)>0.25,gene]

genes.to.plot <- c("FOXA2","FOXC2","FOXB1")

for (i in genes.to.plot) {
  
  outfile <- sprintf("%s/%s_scatterplot_chromRNA_vs_chromVAR_pseudobulk.pdf",io$outdir,i)
  
  to.plot <- rna_chrom.dt[gene==i] %>%
    melt(id.vars=c("celltype","gene","expr"), value.name="motif_accessibility") %>%
    .[,dot_size:=minmax.normalisation(abs(expr))] %>% 
    setorder(-dot_size)

  celltypes_to_label <- rna_tf.dt[gene==i] %>% setorder(-expr) %>% head(n=5) %>% .$celltype
  to.plot.text <- to.plot[celltype%in%celltypes_to_label]
  
  p <- ggscatter(to.plot, x="motif_accessibility", y="expr", fill="celltype", size="dot_size", shape=21, 
                  add="reg.line", add.params = list(color="black", fill="lightgray"), conf.int=TRUE) +
    stat_cor(method = "pearson", label.x.npc = "middle", label.y.npc = "bottom") +
    ggrepel::geom_text_repel(data=to.plot.text, aes(x=motif_accessibility, y=expr, label=celltype), size=3, max.overlaps=Inf) +
    scale_fill_manual(values=opts$celltype.colors) +
    scale_size_continuous(range = c(1,5)) +
    facet_wrap(~variable, nrow=1, scales="free_x") +
    labs(y="RNA expression", x=sprintf("Motif accessibility")) +
    guides(fill=F, size=F) +
    theme(
      axis.text = element_text(size=rel(0.7)),
      axis.title = element_text(size=rel(0.8))
    )
  
  pdf(outfile, width=8, height=3)
  print(p)
  dev.off()
}

