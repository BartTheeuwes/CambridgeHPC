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

io$outdir <- paste0(io$basedir,"/results/rna_atac/rna_vs_chromvar/pseudobulk/per_gene/cross_system_comparison"); dir.create(io$outdir, showWarnings = F)

#############################################
## Load pre-computed correlation estimates ##
#############################################

cor_gastrulation.dt <- fread(paste0(io$basedir,"/results/rna_atac/rna_vs_chromvar/pseudobulk/per_gene/cor_rna_vs_chromvar_correlated_peaks_pseudobulk.txt.gz"))
cor_pbmc.dt <- fread("/Users/ricard/data/multiome_10x_public_data/PBMC/results/rna_atac/rna_vs_chromvar/pseudobulk/per_gene/cor_rna_vs_chromvar_correlated_peaks_pseudobulk.txt.gz")
cor_e18_brain.dt <- fread("/Users/ricard/data/multiome_10x_public_data/e18_mouse_brain/results/rna_atac/rna_vs_chromvar/pseudobulk/per_gene/cor_rna_vs_chromvar_correlated_peaks_pseudobulk.txt.gz")

# (TO-DO Rename genes)

# Merge
cor.dt <- rbind(
  cor_gastrulation.dt[,c("gene","r","padj_fdr","sig")] %>% .[,class:="mouse_organogenesis"],
  # cor_pbmc.dt[,c("gene","r","padj_fdr","sig")] %>% .[,class:="human_pbmc"]
  cor_e18_brain.dt[,c("gene","r","padj_fdr","sig")] %>% .[,class:="mouse_e18_brain"]
) %>% dcast(gene~class, value.var = c("r","padj_fdr","sig"), suffixes=c("_gastrulation","_brain"))


##################
## Calculations ##
##################

tmp <- cor.dt %>% 
  .[abs(r_mouse_e18_brain)>=0.25 & abs(r_mouse_organogenesis)>=0.25]

nrow(tmp[sign(r_mouse_e18_brain)==sign(r_mouse_organogenesis)]) / nrow(tmp)

##########
## Plot ##
##########

to.plot <- cor.dt[!is.na(r_mouse_e18_brain) & !is.na(r_mouse_organogenesis)]
to.plot[,dot_size:=minmax.normalisation(abs(r_mouse_e18_brain)*abs(r_mouse_e18_brain))]

p <- ggscatter(to.plot, x="r_mouse_organogenesis", y="r_mouse_e18_brain", fill="gray70", size="dot_size", alpha=0.75, shape=21) +
  # stat_cor(method = "pearson", label.x.npc = "middle", label.y.npc = "bottom") +
  coord_cartesian(ylim=c(-1,1), xlim=c(-1,1)) +
  scale_size_continuous(range = c(0.5,3)) +
  geom_vline(xintercept=0, linetype="dashed", size=0.5) +
  geom_hline(yintercept=0, linetype="dashed", size=0.5) +
  ggrepel::geom_text_repel(aes(label=gene), size=3.5, max.overlaps=Inf, data=to.plot[r_mouse_organogenesis>0.75 & r_mouse_e18_brain>0.75]) +
  ggrepel::geom_text_repel(aes(label=gene), size=3.5, max.overlaps=Inf, data=to.plot[r_mouse_organogenesis<(-0.65) & r_mouse_e18_brain<(-0.65)]) +
  ggrepel::geom_text_repel(aes(label=gene), size=3.5, max.overlaps=Inf, data=to.plot[r_mouse_organogenesis>0.65 & r_mouse_e18_brain<(-0.65)]) +
  ggrepel::geom_text_repel(aes(label=gene), size=3.5, max.overlaps=Inf, data=to.plot[r_mouse_organogenesis<(-0.65) & r_mouse_e18_brain>0.65]) +
  labs(x="RNA vs chromVAR correlation\n(mouse E7.5-E8.5 organogenesis)", y="RNA vs chromVAR correlation\n(mouse E18.5 brain)") +
  guides(size=F) +
  theme(
    axis.text = element_text(size=rel(0.7))
  )
    
pdf(sprintf("%s/rna_vs_chromvar_mouse_e18_brain_vs_mouse_organogenesis.pdf",io$outdir), width=7, height=6)
print(p)
dev.off()